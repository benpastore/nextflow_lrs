#!/usr/bin/env python3
"""
Recessive-modifier candidate search over the pipeline's consolidated +
snpEff-annotated VCF (one per sample: CALLERS/NUM_CALLERS + ANN, with each
detecting caller's own GT:PS kept in its own sample column -- see
bin/consolidate_variants.py and modules/snpeff/main.nf). Shared parsing,
annotation-scope, and family.json grouping logic lives in
bin/_modifier_common.py (also used by dominant_modifier_consolidated.py) --
this file holds only what's specific to the recessive model.

Unlike bin/recessive_modifier_ont.py (which works from a hand-written list
of discordant mild/severe VCF pairs and only has a co-occurrence heuristic
for compound-het), this script:
  - derives mild/severe sibling groupings directly from family.json
    ({child: {father, mother, phenotype: "mild"|"severe"}})
  - requires putative compound-het pairs to be CONFIRMED in trans (opposite
    haplotypes), using each variant's own calling caller's phase evidence
    (GT + PS for Clair3/Sniffles' read-based local phase blocks; GT alone
    for dipcall/hapdiff, which are phased genome-wide by construction since
    they come directly from a haplotype-resolved assembly)
  - restricts to a narrower annotation scope: CDS-affecting terms, plus
    exon/intron/upstream/downstream -- UTR is deliberately excluded (unlike
    recessive_modifier_ont.py's broader whitelist)
  - restricts further to protein-coding genes only (SnpEff ANN Transcript_
    BioType == 'protein_coding')

Model: recessive only (homozygous, or compound-het confirmed trans) -- for
a dominant-model search (single variant sufficient, no trans confirmation
needed), see dominant_modifier_consolidated.py.

Run against a manifest of per-sample consolidated+snpEff VCF paths:
    recessive_modifier_consolidated.py \\
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

# Read-based callers: haplotype comparisons only valid within the same PS
# (phase-block) value. Assembly-based callers: GT is phased genome-wide by
# construction (haplotype-resolved assembly), no PS/block concept needed.
LOCAL_PHASE_CALLERS = {'clair3_longphase', 'clair3', 'sniffles'}
ASSEMBLY_CALLERS = {'dipcall', 'hapdiff'}


def is_phased(gt):
    return '|' in gt


def alt_hap_index(gt):
    """For a phased het GT ('1|0'/'0|1'), which haplotype side (0/1) carries ALT."""
    alleles = gt.split('|')
    return 0 if alleles[0] != '0' else 1


def evaluate_mild_sample(records, small_keys, sv_positions, sv_merge_dist):
    """
    Apply the recessive model to one mild sample's variant set.

    Returns {(var_id, gene): {"reason":..., "gene":..., "rec":...}} for
    variants that pass: homozygous (not in severe), or heterozygous with a
    confirmed-trans compound-het partner in the same gene (not in severe).
    Keyed by (var_id, gene) rather than var_id alone so a variant spanning
    multiple genes (e.g. upstream/downstream of two neighboring genes, or a
    genuinely overlapping gene pair) gets a separate passing entry per gene
    instead of silently collapsing onto just one.
    """
    surviving = {
        vid: rec for vid, rec in records.items()
        if not in_severe(rec, small_keys, sv_positions, sv_merge_dist)
    }

    passing = {}

    for vid, rec in surviving.items():
        if rec["zygosity"] == "homo":
            for gene in rec["genes"]:
                passing[(vid, gene)] = {"reason": "homozygous", "gene": gene, "rec": rec}

    het_by_gene = defaultdict(list)
    for vid, rec in surviving.items():
        if rec["zygosity"] != "het":
            continue
        for gene in rec["genes"]:
            het_by_gene[gene].append((vid, rec))

    for gene, items in het_by_gene.items():
        local = [
            (vid, rec) for vid, rec in items
            if rec["caller"] in LOCAL_PHASE_CALLERS and is_phased(rec["gt"]) and rec["ps"] not in ("", ".")
        ]
        asm = [
            (vid, rec) for vid, rec in items
            if rec["caller"] in ASSEMBLY_CALLERS and is_phased(rec["gt"])
        ]

        by_ps = defaultdict(list)
        for vid, rec in local:
            by_ps[rec["ps"]].append((vid, rec))
        for ps, group in by_ps.items():
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    vid1, rec1 = group[i]
                    vid2, rec2 = group[j]
                    if alt_hap_index(rec1["gt"]) != alt_hap_index(rec2["gt"]):
                        passing.setdefault((vid1, gene), {
                            "reason": f"compound_het_trans:PS={ps},partner={vid2}", "gene": gene, "rec": rec1
                        })
                        passing.setdefault((vid2, gene), {
                            "reason": f"compound_het_trans:PS={ps},partner={vid1}", "gene": gene, "rec": rec2
                        })

        for i in range(len(asm)):
            for j in range(i + 1, len(asm)):
                vid1, rec1 = asm[i]
                vid2, rec2 = asm[j]
                if alt_hap_index(rec1["gt"]) != alt_hap_index(rec2["gt"]):
                    passing.setdefault((vid1, gene), {
                        "reason": f"compound_het_trans:assembly,partner={vid2}", "gene": gene, "rec": rec1
                    })
                    passing.setdefault((vid2, gene), {
                        "reason": f"compound_het_trans:assembly,partner={vid1}", "gene": gene, "rec": rec2
                    })

    return passing


def severe_qualifying_genes(severe_records_list, sv_merge_dist):
    """
    Genes where >=1 severe sample itself carries the recessive model's own
    qualifying pattern -- homozygous, or a compound-het pair confirmed in
    trans -- evaluated with no exclusion (severe isn't being filtered
    against anything here, just characterized on its own terms; passing
    small_keys=set() and sv_positions={} to evaluate_mild_sample means its
    in_severe() check is trivially False for every variant, so nothing is
    dropped before the homo/compound-het logic runs).

    Used to drop an entire gene from a mild sample's candidates when severe
    itself is also biallelically disrupted there -- even via a *different*
    specific variant than mild's. A true recessive model: if the severe
    sibling already carries a qualifying hit in this gene, that gene can't
    be what's differentiating mild from severe, so exact-variant matching
    alone (in_severe/build_severe_index above) isn't enough -- this closes
    that gap at the gene level.
    """
    genes = set()
    for records in severe_records_list:
        passing = evaluate_mild_sample(records, set(), {}, sv_merge_dist)
        genes.update(info["gene"] for info in passing.values())
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
    variant_sample_map = defaultdict(set)     # (var_id, gene) -> set of mild sampleIDs where it passed
    gene_sample_map = defaultdict(set)        # (mild_sample, gene) already tracked implicitly via variant_sample_map
    gene_to_samples = defaultdict(set)        # gene -> set of mild sampleIDs with ANY qualifying pattern in it
    all_passing = {}                          # (var_id, gene) -> {reason, gene, rec, samples:set}

    for mild, severes in discordant_groups:
        mild_records = get_records(mild)
        severe_records = [get_records(s) for s in severes] + concordant_records
        small_keys, sv_positions = build_severe_index(severe_records)

        passing = evaluate_mild_sample(mild_records, small_keys, sv_positions, sv_merge_dist)

        excluded_genes = severe_qualifying_genes(severe_records, sv_merge_dist)
        if excluded_genes:
            before = len(passing)
            passing = {k: v for k, v in passing.items() if v["gene"] not in excluded_genes}
            dropped = before - len(passing)
            if dropped:
                logger.info(f"{mild}: dropped {dropped} candidate(s) in {len(excluded_genes)} gene(s) "
                            f"also qualifying (homo/compound-het-trans) in severe")

        logger.info(f"{mild} (vs {', '.join(severes)} + {len(concordant_severe)} concordant-severe): "
                    f"{len(passing)} candidate variants")

        for key, info in passing.items():
            variant_sample_map[key].add(mild)
            gene_to_samples[info["gene"]].add(mild)
            if key not in all_passing:
                all_passing[key] = {"reason": info["reason"], "gene": info["gene"], "rec": info["rec"], "samples": set()}
            all_passing[key]["samples"].add(mild)

    logger.info(f"\nCandidate variants across all mild samples (pre multi-sample filter): {len(all_passing)}")

    # Phase 2: require recurrence across >=min_occurrence independent mild
    # samples -- either the same variant (homozygous), or (for compound-het)
    # >=min_occurrence samples each showing a qualifying trans pattern in
    # the same gene. Both branches use the same threshold now -- previously
    # the homozygous branch hardcoded >=2 regardless of -min_occurrence,
    # which meant a cohort with a single mild sample could never report any
    # homozygous candidate no matter how real the finding was.
    final_pass = {}
    for key, info in all_passing.items():
        if info["reason"] == "homozygous":
            if len(info["samples"]) >= min_occurrence:
                final_pass[key] = info
        else:
            if len(gene_to_samples[info["gene"]]) >= min_occurrence:
                final_pass[key] = info

    logger.info(f"After multi-sample recurrence filter (>={min_occurrence} independent mild samples): {len(final_pass)}")

    if not final_pass:
        logger.warning("No variants passed filtering - returning empty results table")
        empty_variants = pd.DataFrame(columns=[
            "#CHROM", "POS", "ID", "REF", "ALT", "gene", "zygosity", "reason",
            "evidence_caller", "callers_detected", "GT", "PS", "samples", "N_samples", "gnomad_AF_NFE"
        ])
        empty_genes = pd.DataFrame(columns=["gene", "N_sibs_with_variant", "N_homo_sibs", "N_compound_het_sibs", "samples"])
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
    gene_homo_samples = defaultdict(set)
    gene_comphet_samples = defaultdict(set)
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
        if info["reason"] == "homozygous":
            gene_homo_samples[gene].update(info["samples"])
        else:
            gene_comphet_samples[gene].update(info["samples"])

    final_df = pd.DataFrame(rows).sort_values(["#CHROM", "POS"]).reset_index(drop=True)

    N_homo = (final_df["zygosity"] == "homo").sum()
    N_het = (final_df["zygosity"] == "het").sum()
    logger.info(f"\n{'='*60}")
    logger.info(f"FINAL RESULTS")
    logger.info(f"{'='*60}")
    logger.info(f"Total variants: {len(final_df)} (homo: {N_homo}, het/compound-het-trans: {N_het})")
    logger.info(f"Unique genes: {final_df['gene'].nunique()}")
    logger.info(f"{'='*60}\n")

    # Per-gene summary: how many independent mild sibs carry a qualifying
    # variant (homo and/or compound-het-trans) in each gene -- a sib
    # counts once per gene even if it has multiple qualifying variants
    # there, and once even if it qualifies via both homo and compound-het.
    gene_summary_rows = []
    all_genes = set(gene_homo_samples) | set(gene_comphet_samples)
    for gene in all_genes:
        homo_samples = gene_homo_samples[gene]
        comphet_samples = gene_comphet_samples[gene]
        all_samples = homo_samples | comphet_samples
        gene_summary_rows.append({
            "gene": gene,
            "N_sibs_with_variant": len(all_samples),
            "N_homo_sibs": len(homo_samples),
            "N_compound_het_sibs": len(comphet_samples),
            "samples": ",".join(sorted(all_samples)),
        })
    gene_summary_df = pd.DataFrame(gene_summary_rows).sort_values(
        ["N_sibs_with_variant", "gene"], ascending=[False, True]
    ).reset_index(drop=True)

    return final_df, gene_summary_df


def get_args():
    parser = argparse.ArgumentParser(
        description="Recessive-modifier search over the consolidated + snpEff VCF, with phase-confirmed compound-het."
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

    logger.info("Starting recessive-modifier search over consolidated+snpEff VCFs")

    discordant_groups, concordant_severe = load_family_groups(args.family_json)
    manifest = load_vcf_manifest(args.vcf_manifest)

    gnomad_annotator = None
    if args.gnomad_chm13_vcf:
        gnomad_annotator = GnomadNFEAnnotator(args.gnomad_chm13_vcf)

    res, gene_summary = run(discordant_groups, concordant_severe, manifest, args.sv_merge_dist, gnomad_annotator, args.min_occurrence)

    output_file = args.output or f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_recessive_modifier_consolidated.tsv"
    res.to_csv(output_file, index=False, sep="\t", header=True)
    logger.info(f"\nResults saved to: {output_file}")

    gene_summary_file = output_file.rsplit(".", 1)[0] + ".gene_summary.tsv"
    gene_summary.to_csv(gene_summary_file, index=False, sep="\t", header=True)
    logger.info(f"Gene summary saved to: {gene_summary_file}")


if __name__ == "__main__":
    main()
