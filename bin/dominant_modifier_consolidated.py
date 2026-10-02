#!/usr/bin/env python3
"""
Dominant-modifier candidate search over the pipeline's consolidated +
snpEff-annotated VCF -- same inputs, annotation scope, and family.json-
derived mild/severe sibling groupings as recessive_modifier_consolidated.py
(shared code lives in bin/_modifier_common.py), but a single variant is
sufficient on its own under a dominant model: no second hit, no compound-
het trans confirmation, and no distinction between het/homo zygosity for
whether a variant qualifies (only in what it's reported as).

Gene-level severe exclusion is stricter here than the recessive script's:
since one variant is already sufficient, a gene is excluded from mild's
candidates if severe carries *any* variant there at all (passing the same
annotation-scope + protein-coding filter), not just a homo/compound-het
pattern -- there's no weaker "qualifying pattern" concept to require for a
model where a single hit is already enough.

Run against a manifest of per-sample consolidated+snpEff VCF paths:
    dominant_modifier_consolidated.py \\
        -family_json family.json \\
        -vcf_manifest vcf_manifest.tsv \\
        -output results.tsv
"""

import argparse
from collections import defaultdict
from datetime import datetime

import pandas as pd

from _modifier_common import (
    GnomadNFEAnnotator,
    build_severe_index,
    in_severe,
    load_family_groups,
    load_vcf_manifest,
    logger,
    parse_consolidated_vcf,
    timing,
)


def evaluate_mild_sample(records, small_keys, sv_positions, sv_merge_dist):
    """
    Apply the dominant model to one mild sample's variant set: any variant
    (het or homo) not in severe is sufficient on its own.

    Returns {(var_id, gene): {"reason":..., "gene":..., "rec":...}}, keyed
    by (var_id, gene) so a variant spanning multiple genes gets a separate
    passing entry per gene instead of collapsing onto just one.
    """
    surviving = {
        vid: rec for vid, rec in records.items()
        if not in_severe(rec, small_keys, sv_positions, sv_merge_dist)
    }

    passing = {}
    for vid, rec in surviving.items():
        for gene in rec["genes"]:
            passing[(vid, gene)] = {"reason": f"dominant_{rec['zygosity']}", "gene": gene, "rec": rec}
    return passing


def severe_qualifying_genes(severe_records_list):
    """
    Genes where >=1 severe sample carries ANY variant at all (already
    filtered to protein-coding + allowed annotation terms during parsing).
    Stricter than the recessive script's version by necessity -- under a
    dominant model there's no weaker bar than "has a variant here" to
    define "severe already has a qualifying hit in this gene", so any
    variant in severe (even a different one than mild's) rules the gene
    out entirely.
    """
    genes = set()
    for records in severe_records_list:
        for rec in records.values():
            genes.update(rec["genes"])
    return genes


@timing
def run(discordant_groups, concordant_severe, manifest, sv_merge_dist, gnomad_annotator=None, min_occurrence=1):
    logger.info(f"\n{'='*60}")
    logger.info(f"Discordant mild/severe groups: {len(discordant_groups)}")
    logger.info(f"Concordant-severe background samples: {len(concordant_severe)}")
    logger.info(f"{'='*60}\n")

    parsed_cache = {}

    def get_records(sample):
        if sample not in parsed_cache:
            if sample not in manifest:
                logger.warning(f"{sample}: no VCF in manifest, skipping")
                parsed_cache[sample] = {}
            else:
                parsed_cache[sample] = parse_consolidated_vcf(manifest[sample])
        return parsed_cache[sample]

    concordant_records = [get_records(s) for s in concordant_severe]

    # Phase 1: per-mild-sample evaluation
    gene_to_samples = defaultdict(set)        # gene -> set of mild sampleIDs with ANY qualifying variant in it
    all_passing = {}                          # (var_id, gene) -> {reason, gene, rec, samples:set}

    for mild, severes in discordant_groups:
        mild_records = get_records(mild)
        severe_records = [get_records(s) for s in severes] + concordant_records
        small_keys, sv_positions = build_severe_index(severe_records)

        passing = evaluate_mild_sample(mild_records, small_keys, sv_positions, sv_merge_dist)

        excluded_genes = severe_qualifying_genes(severe_records)
        if excluded_genes:
            before = len(passing)
            passing = {k: v for k, v in passing.items() if v["gene"] not in excluded_genes}
            dropped = before - len(passing)
            if dropped:
                logger.info(f"{mild}: dropped {dropped} candidate(s) in {len(excluded_genes)} gene(s) "
                            f"that also have a variant in severe")

        logger.info(f"{mild} (vs {', '.join(severes)} + {len(concordant_severe)} concordant-severe): "
                    f"{len(passing)} candidate variants")

        for key, info in passing.items():
            gene_to_samples[info["gene"]].add(mild)
            if key not in all_passing:
                all_passing[key] = {"reason": info["reason"], "gene": info["gene"], "rec": info["rec"], "samples": set()}
            all_passing[key]["samples"].add(mild)

    logger.info(f"\nCandidate variants across all mild samples (pre multi-sample filter): {len(all_passing)}")

    # Phase 2: require recurrence across >=min_occurrence independent mild
    # samples showing a qualifying variant in the same gene (allelic
    # recurrence of the exact same variant isn't required -- locus
    # heterogeneity, different private variants in the same gene across
    # different mild sibs, is expected under a dominant model same as it
    # is for the recessive script's compound-het path).
    final_pass = {}
    for key, info in all_passing.items():
        if len(gene_to_samples[info["gene"]]) >= min_occurrence:
            final_pass[key] = info

    logger.info(f"After multi-sample recurrence filter (>={min_occurrence} independent mild samples): {len(final_pass)}")

    if not final_pass:
        logger.warning("No variants passed filtering - returning empty results table")
        empty_variants = pd.DataFrame(columns=[
            "#CHROM", "POS", "ID", "REF", "ALT", "gene", "zygosity", "reason",
            "evidence_caller", "callers_detected", "GT", "PS", "samples", "N_samples", "gnomad_AF_NFE"
        ])
        empty_genes = pd.DataFrame(columns=["gene", "N_sibs_with_variant", "samples"])
        return empty_variants, empty_genes

    var_id_to_af = {}
    if gnomad_annotator is not None:
        logger.info("Annotating passing variants with gnomAD AF_nfe...")
        gnomad_annotator.load()
        batch = [(key, info["rec"]["chrom"], info["rec"]["pos"], info["rec"]["ref"], info["rec"]["alt"])
                 for key, info in final_pass.items()]
        results = gnomad_annotator.annotate_batch([(c, p, r, a) for _, c, p, r, a in batch])
        for key, chrom, pos, ref, alt in batch:
            var_id_to_af[key] = results.get(f"{chrom}:{pos}:{ref}>{alt}")

    rows = []
    gene_samples = defaultdict(set)
    for key, info in final_pass.items():
        vid, gene = key
        rec = info["rec"]
        rows.append({
            "#CHROM": rec["chrom"],
            "POS": rec["pos"],
            "ID": vid,
            "REF": rec["ref"],
            "ALT": rec["alt"],
            "gene": gene,
            "zygosity": rec["zygosity"],
            "reason": info["reason"],
            "evidence_caller": rec["caller"],
            "callers_detected": ",".join(rec["callers_all"]),
            "GT": rec["gt"],
            "PS": rec["ps"],
            "samples": ",".join(sorted(info["samples"])),
            "N_samples": len(info["samples"]),
            "gnomad_AF_NFE": var_id_to_af.get(key),
        })
        gene_samples[gene].update(info["samples"])

    final_df = pd.DataFrame(rows).sort_values(["#CHROM", "POS"]).reset_index(drop=True)

    N_homo = (final_df["zygosity"] == "homo").sum()
    N_het = (final_df["zygosity"] == "het").sum()
    logger.info(f"\n{'='*60}")
    logger.info(f"FINAL RESULTS")
    logger.info(f"{'='*60}")
    logger.info(f"Total variants: {len(final_df)} (homo: {N_homo}, het: {N_het})")
    logger.info(f"Unique genes: {final_df['gene'].nunique()}")
    logger.info(f"{'='*60}\n")

    gene_summary_rows = [
        {"gene": gene, "N_sibs_with_variant": len(samples), "samples": ",".join(sorted(samples))}
        for gene, samples in gene_samples.items()
    ]
    gene_summary_df = pd.DataFrame(gene_summary_rows).sort_values(
        ["N_sibs_with_variant", "gene"], ascending=[False, True]
    ).reset_index(drop=True)

    return final_df, gene_summary_df


def get_args():
    parser = argparse.ArgumentParser(
        description="Dominant-modifier search over the consolidated + snpEff VCF: a single variant is sufficient, no compound-het needed."
    )
    parser.add_argument("-family_json", type=str, required=True,
                         help="Path to family.json: {child: {father, mother, phenotype: mild|severe}}")
    parser.add_argument("-vcf_manifest", type=str, required=True,
                         help="TSV: sampleID<TAB>path to that sample's consolidated+snpEff VCF")
    parser.add_argument("-sv_merge_dist", type=int, default=500,
                         help="Max bp between SV start positions (same chrom+SVTYPE) to treat as the same site across samples")
    parser.add_argument("-gnomad_chm13_vcf", type=str, default=None,
                         help="Optional gnomAD VCF already lifted onto CHM13 coordinates, for AF_nfe annotation of the final pass")
    parser.add_argument("-output", type=str, default=None)
    parser.add_argument("-min_occurrence", type=int, default=1,
                         help="Minimum number of independent mild samples a gene must appear in to pass the multi-sample recurrence filter")
    return parser.parse_args()


def main():
    args = get_args()

    logger.info("Starting dominant-modifier search over consolidated+snpEff VCFs")

    discordant_groups, concordant_severe = load_family_groups(args.family_json)
    manifest = load_vcf_manifest(args.vcf_manifest)

    gnomad_annotator = None
    if args.gnomad_chm13_vcf:
        gnomad_annotator = GnomadNFEAnnotator(args.gnomad_chm13_vcf)

    res, gene_summary = run(discordant_groups, concordant_severe, manifest, args.sv_merge_dist, gnomad_annotator, args.min_occurrence)

    output_file = args.output or f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_dominant_modifier_consolidated.tsv"
    res.to_csv(output_file, index=False, sep="\t", header=True)
    logger.info(f"\nResults saved to: {output_file}")

    gene_summary_file = output_file.rsplit(".", 1)[0] + ".gene_summary.tsv"
    gene_summary.to_csv(gene_summary_file, index=False, sep="\t", header=True)
    logger.info(f"Gene summary saved to: {gene_summary_file}")


if __name__ == "__main__":
    main()
