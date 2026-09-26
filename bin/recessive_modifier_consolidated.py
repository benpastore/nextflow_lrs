#!/usr/bin/env python3
"""
Recessive-modifier candidate search over the pipeline's consolidated +
snpEff-annotated VCF (one per sample: CALLERS/NUM_CALLERS + ANN, with each
detecting caller's own GT:PS kept in its own sample column -- see
bin/consolidate_variants.py and modules/snpeff/main.nf).

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

Model: recessive only (homozygous, or compound-het confirmed trans) -- for
a dominant-model search, see recessive_modifier_ont.py.

Run against a manifest of per-sample consolidated+snpEff VCF paths:
    recessive_modifier_consolidated.py \\
        -family_json family.json \\
        -vcf_manifest vcf_manifest.tsv \\
        -output results.tsv
"""

import argparse
import bisect
import gzip
import json
import logging
import os
import pickle
import time
import datetime
from collections import defaultdict

import numpy as np
import pandas as pd

try:
    import colorlog
except ImportError:
    colorlog = None

timestamp = datetime.datetime.now().strftime("%y%m%d_%H%M%S")

# Configure logging -- colorlog is a nice-to-have (not every environment
# this runs in has it, e.g. the rnaseq container); fall back to plain
# logging rather than failing outright.
level = "INFO"
logger = logging.getLogger("color_logger")

if not logger.handlers:
    if colorlog is not None:
        handler = colorlog.StreamHandler()
        handler.setFormatter(colorlog.ColoredFormatter(
            "%(asctime)s - %(log_color)s%(levelname)s:%(reset)s %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            log_colors={
                'DEBUG': 'blue',
                'INFO': 'green',
                'WARNING': 'yellow',
                'ERROR': 'red',
                'CRITICAL': 'bold_red',
            }
        ))
    else:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(
            "%(asctime)s - %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        ))
    logger.addHandler(handler)
    logger.setLevel(getattr(logging, level, logging.INFO))


def timing(f):
    """Decorator for timing functions."""
    def wrap(*args, **kwargs):
        time1 = time.time()
        ret = f(*args, **kwargs)
        time2 = time.time()
        logger.info(f"{f.__name__} function took {time2 - time1:.3f} s")
        return ret
    return wrap


def opener(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


class GnomadNFEAnnotator:
    """
    Lazy-loading annotator for gnomAD non-Finnish European (AF_nfe) allele
    frequencies, sourced from a gnomAD VCF already lifted onto the T2T-CHM13
    assembly (e.g. Ensembl's GCA_009914755.4 rapid-release gnomAD VCF) --
    since our variants are called against CHM13/hs1, no liftOver is needed
    as long as the resource is already in that coordinate system.
    """

    def __init__(self, gnomad_vcf_path=None):
        self.gnomad_vcf_path = gnomad_vcf_path
        self.af_dict = None

    @timing
    def load(self):
        """Load (or build+cache) the CHROM:POS:REF>ALT -> AF_nfe lookup."""
        if self.af_dict is not None:
            return  # Already loaded

        if not self.gnomad_vcf_path:
            logger.warning("No gnomAD CHM13 VCF provided")
            return

        pickle_file = f"{self.gnomad_vcf_path}.pickle"

        if os.path.exists(pickle_file):
            logger.info(f"Loading cached gnomAD AF_nfe lookup from {pickle_file}...")
            with open(pickle_file, 'rb') as f:
                self.af_dict = pickle.load(f)
            logger.info(f"gnomAD AF_nfe lookup loaded: {len(self.af_dict):,} variants")
            return

        logger.info(f"Building gnomAD AF_nfe lookup from {self.gnomad_vcf_path}...")

        af_dict = {}
        with opener(self.gnomad_vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")
                chrom, pos, ref, alt, info = fields[0], fields[1], fields[3], fields[4], fields[7]

                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"

                af_nfe = None
                for entry in info.split(';'):
                    if entry.startswith('AF_nfe='):
                        af_nfe = entry.split('=', 1)[1]
                        break

                if af_nfe is None:
                    continue

                alts = alt.split(',')
                afs = af_nfe.split(',')
                if len(alts) != len(afs):
                    continue

                for a, af in zip(alts, afs):
                    try:
                        af_dict[f"{chrom}:{pos}:{ref}>{a}"] = float(af)
                    except ValueError:
                        continue

        self.af_dict = af_dict

        logger.info(f"Caching gnomAD AF_nfe lookup to {pickle_file} for faster future loading...")
        with open(pickle_file, 'wb') as f:
            pickle.dump(self.af_dict, f)

        logger.info(f"gnomAD AF_nfe lookup built: {len(self.af_dict):,} variants")

    def annotate_batch(self, variants):
        """variants: list of (chrom, pos, ref, alt) -> dict variant_id -> AF_nfe or None"""
        if self.af_dict is None:
            return {f"{c}:{p}:{r}>{a}": None for c, p, r, a in variants}

        results = {}
        for chrom, pos, ref, alt in variants:
            variant_id = f"{chrom}:{pos}:{ref}>{alt}"
            results[variant_id] = self.af_dict.get(variant_id, None)
        return results


# CDS-affecting terms -- SnpEff doesn't tag most coding changes with a
# literal "cds" term, so "cds" is interpreted as this whole consequence set.
CDS_TERMS = {
    'missense_variant',
    'synonymous_variant',
    'stop_gained',
    'stop_lost',
    'start_lost',
    'frameshift_variant',
    'inframe_insertion',
    'inframe_deletion',
    'protein_altering_variant',
    'splice_donor_variant',
    'splice_acceptor_variant',
    'splice_region_variant',
}
ALLOWED_ANNOTATIONS = CDS_TERMS | {
    'exon_variant',
    'intron_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
}

# Read-based callers: haplotype comparisons only valid within the same PS
# (phase-block) value. Assembly-based callers: GT is phased genome-wide by
# construction (haplotype-resolved assembly), no PS/block concept needed.
LOCAL_PHASE_CALLERS = {'clair3_longphase', 'clair3', 'sniffles'}
ASSEMBLY_CALLERS = {'dipcall', 'hapdiff'}
CALLER_PRIORITY = ['clair3_longphase', 'clair3', 'sniffles', 'dipcall', 'hapdiff']


def get_ann(info_str):
    """Extract and parse the ANN field from VCF INFO string."""
    for field in info_str.split(';'):
        if field.startswith('ANN='):
            ann = field.split('=', 1)[1]
            return [x.split("|") for x in ann.split(",")]
    return False


def get_info_field(info_str, key):
    prefix = f"{key}="
    for field in info_str.split(';'):
        if field.startswith(prefix):
            return field[len(prefix):]
    return None


def classify_gt(gt):
    """Return 'homo', 'het', or None (not called / ref-only) for a GT string."""
    if not gt or gt in ('.', './.', '.|.'):
        return None
    sep = '|' if '|' in gt else '/'
    alleles = gt.split(sep)
    if len(alleles) != 2 or '.' in alleles:
        return None
    if alleles[0] == alleles[1]:
        return 'homo' if alleles[0] != '0' else None
    return 'het'


def is_phased(gt):
    return '|' in gt


def alt_hap_index(gt):
    """For a phased het GT ('1|0'/'0|1'), which haplotype side (0/1) carries ALT."""
    alleles = gt.split('|')
    return 0 if alleles[0] != '0' else 1


@timing
def parse_consolidated_vcf(path):
    """
    Parse one sample's consolidated+snpEff VCF into {var_id: record}.

    record: chrom, pos, ref, alt, svtype, ann (filtered), genes (set),
    caller (the one whose GT/PS we're trusting), callers_all (list),
    gt, ps, zygosity.
    """
    results = {}
    caller_col_idx = {}
    kept = 0
    filtered_ann_type = 0
    filtered_no_caller = 0
    filtered_not_called = 0

    with opener(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                for i, name in enumerate(cols):
                    if i >= 9:
                        caller_col_idx[name] = i
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, pos, vid, ref, alt, _qual, _filt, info, fmt = fields[:9]

            ann = get_ann(info)
            if not ann:
                continue
            filtered_ann = [a for a in ann if len(a) > 1 and any(t in a[1] for t in ALLOWED_ANNOTATIONS)]
            if not filtered_ann:
                filtered_ann_type += 1
                continue

            callers_field = get_info_field(info, "CALLERS")
            callers_present = callers_field.split(",") if callers_field else []
            chosen = next((c for c in CALLER_PRIORITY if c in callers_present and c in caller_col_idx), None)
            if chosen is None:
                filtered_no_caller += 1
                continue

            sample_val = fields[caller_col_idx[chosen]]
            fmt_keys = fmt.split(":")
            sample_vals = sample_val.split(":")
            d = dict(zip(fmt_keys, sample_vals))
            gt = d.get("GT", "")
            ps = d.get("PS", "")

            zygosity = classify_gt(gt)
            if zygosity is None:
                filtered_not_called += 1
                continue

            svtype = get_info_field(info, "SVTYPE")
            genes = {a[3] for a in filtered_ann if len(a) > 3 and a[3]}

            results[vid] = {
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "svtype": svtype,
                # filtered_ann itself is never read past this point -- only
                # the gene set it collapses to below is used downstream, so
                # don't hold onto the (potentially large, one-entry-per-
                # transcript) parsed ANN block for every kept variant.
                "genes": genes,
                "caller": chosen,
                "callers_all": callers_present,
                "gt": gt,
                "ps": ps,
                "zygosity": zygosity,
            }
            kept += 1

    logger.debug(
        f"Parsed {os.path.basename(path)}: kept {kept}, "
        f"filtered {filtered_ann_type} (annotation type), "
        f"{filtered_no_caller} (no usable caller column), "
        f"{filtered_not_called} (not called by chosen caller)"
    )
    return results


def build_severe_index(severe_records_list):
    """severe_records_list: list of per-sample {var_id: record} dicts.
    Returns (small_variant_keys, sv_positions_by_chrom_type)."""
    small_keys = set()
    sv_positions = defaultdict(list)
    for records in severe_records_list:
        for rec in records.values():
            if rec["svtype"]:
                sv_positions[(rec["chrom"], rec["svtype"])].append(rec["pos"])
            else:
                small_keys.add((rec["chrom"], rec["pos"], rec["ref"], rec["alt"]))
    for k in sv_positions:
        sv_positions[k].sort()
    return small_keys, sv_positions


def in_severe(rec, small_keys, sv_positions, sv_merge_dist):
    if rec["svtype"]:
        positions = sv_positions.get((rec["chrom"], rec["svtype"]))
        if not positions:
            return False
        i = bisect.bisect_left(positions, rec["pos"])
        for j in (i - 1, i):
            if 0 <= j < len(positions) and abs(positions[j] - rec["pos"]) <= sv_merge_dist:
                return True
        return False
    return (rec["chrom"], rec["pos"], rec["ref"], rec["alt"]) in small_keys


def evaluate_mild_sample(records, small_keys, sv_positions, sv_merge_dist):
    """
    Apply the recessive model to one mild sample's variant set.

    Returns {var_id: {"reason":..., "gene":...}} for variants that pass:
    homozygous (not in severe), or heterozygous with a confirmed-trans
    compound-het partner in the same gene (not in severe).
    """
    surviving = {
        vid: rec for vid, rec in records.items()
        if not in_severe(rec, small_keys, sv_positions, sv_merge_dist)
    }

    passing = {}

    for vid, rec in surviving.items():
        if rec["zygosity"] == "homo":
            for gene in rec["genes"]:
                passing[vid] = {"reason": "homozygous", "gene": gene, "rec": rec}

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
                        passing.setdefault(vid1, {
                            "reason": f"compound_het_trans:PS={ps},partner={vid2}", "gene": gene, "rec": rec1
                        })
                        passing.setdefault(vid2, {
                            "reason": f"compound_het_trans:PS={ps},partner={vid1}", "gene": gene, "rec": rec2
                        })

        for i in range(len(asm)):
            for j in range(i + 1, len(asm)):
                vid1, rec1 = asm[i]
                vid2, rec2 = asm[j]
                if alt_hap_index(rec1["gt"]) != alt_hap_index(rec2["gt"]):
                    passing.setdefault(vid1, {
                        "reason": f"compound_het_trans:assembly,partner={vid2}", "gene": gene, "rec": rec1
                    })
                    passing.setdefault(vid2, {
                        "reason": f"compound_het_trans:assembly,partner={vid1}", "gene": gene, "rec": rec2
                    })

    return passing


def load_family_groups(family_json_path):
    """
    Parse family.json ({child: {father, mother, phenotype}}) into:
      discordant_groups: list of (mild_sampleID, [severe_sampleIDs]) -- one
                          entry per mild child, paired against every severe
                          child sharing its (father, mother)
      concordant_severe: list of severe sampleIDs whose family has no mild
                          child (background filtering cohort, same role as
                          recessive_modifier_ont.py's -concordant list) --
                          also where a child with no father/mother in
                          family.json lands if it's severe (see below)

    A child missing father and/or mother in family.json never participates
    in the (father, mother) grouping below -- without both, there's no
    genuine trio/sibling relationship to establish, and keying on a
    missing value (e.g. both absent -> key (None, None)) would silently
    treat unrelated no-parent children as siblings of each other. Instead:
      - severe with no parents -> added directly to concordant_severe
        (its role there, filtering background evidence, needs no sibling
        relationship, so this is a safe default, not a loss of signal)
      - mild with no parents -> skipped entirely (no severe sibling can be
        established, so no trio comparison is possible for it); logged as
        a warning since this silently drops that sample from the search
    """
    with open(family_json_path) as f:
        families = json.load(f)

    by_parents = defaultdict(lambda: {"mild": [], "severe": []})
    concordant_severe = []
    for child, info in families.items():
        phenotype = info.get("phenotype")
        if phenotype not in ("mild", "severe"):
            logger.warning(f"{child}: phenotype '{phenotype}' is neither 'mild' nor 'severe', skipping")
            continue

        father = info.get("father")
        mother = info.get("mother")
        if not father or not mother:
            if phenotype == "severe":
                logger.info(
                    f"{child}: no father/mother in family.json -- adding "
                    f"directly to concordant-severe background (bypassing "
                    f"trio/sibling grouping)"
                )
                concordant_severe.append(child)
            else:
                logger.warning(
                    f"{child}: mild phenotype with no father/mother in "
                    f"family.json -- cannot establish a severe sibling to "
                    f"compare against, skipping (no trio analysis possible)"
                )
            continue

        by_parents[(father, mother)][phenotype].append(child)

    discordant_groups = []
    for (father, mother), grp in by_parents.items():
        if grp["mild"] and grp["severe"]:
            for mild in grp["mild"]:
                discordant_groups.append((mild, list(grp["severe"])))
        elif grp["severe"]:
            concordant_severe.extend(grp["severe"])
        elif grp["mild"]:
            logger.warning(
                f"Family (father={father}, mother={mother}): mild child(ren) "
                f"{grp['mild']} have no severe sibling -- no within-family "
                f"severe comparison possible for them (concordant-severe "
                f"background cohort still applies)"
            )

    return discordant_groups, concordant_severe


def load_vcf_manifest(path):
    """sampleID<TAB>vcf_path, one per line -> {sampleID: vcf_path}"""
    manifest = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sample, vcf = line.split("\t")
            manifest[sample] = vcf
    return manifest


@timing
def run(discordant_groups, concordant_severe, manifest, sv_merge_dist, gnomad_annotator=None):
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
    variant_sample_map = defaultdict(set)     # var_id -> set of mild sampleIDs where it passed
    gene_sample_map = defaultdict(set)        # (mild_sample, gene) already tracked implicitly via variant_sample_map
    gene_to_samples = defaultdict(set)        # gene -> set of mild sampleIDs with ANY qualifying pattern in it
    all_passing = {}                          # var_id -> {reason, gene, rec, samples:set}

    for mild, severes in discordant_groups:
        mild_records = get_records(mild)
        severe_records = [get_records(s) for s in severes] + concordant_records
        small_keys, sv_positions = build_severe_index(severe_records)

        passing = evaluate_mild_sample(mild_records, small_keys, sv_positions, sv_merge_dist)
        logger.info(f"{mild} (vs {', '.join(severes)} + {len(concordant_severe)} concordant-severe): "
                    f"{len(passing)} candidate variants")

        for vid, info in passing.items():
            variant_sample_map[vid].add(mild)
            gene_to_samples[info["gene"]].add(mild)
            if vid not in all_passing:
                all_passing[vid] = {"reason": info["reason"], "gene": info["gene"], "rec": info["rec"], "samples": set()}
            all_passing[vid]["samples"].add(mild)

    logger.info(f"\nCandidate variants across all mild samples (pre multi-sample filter): {len(all_passing)}")

    # Phase 2: require recurrence across >=2 independent mild samples --
    # either the same variant, or (for compound-het) >=2 samples each
    # showing a qualifying trans pattern in the same gene.
    final_pass = {}
    for vid, info in all_passing.items():
        if info["reason"] == "homozygous":
            if len(info["samples"]) >= 2:
                final_pass[vid] = info
        else:
            if len(gene_to_samples[info["gene"]]) >= 2:
                final_pass[vid] = info

    logger.info(f"After multi-sample recurrence filter (>=2 independent mild samples): {len(final_pass)}")

    if not final_pass:
        logger.warning("No variants passed filtering - returning empty results table")
        return pd.DataFrame(columns=[
            "#CHROM", "POS", "ID", "REF", "ALT", "gene", "zygosity", "reason",
            "evidence_caller", "callers_detected", "GT", "PS", "samples", "N_samples", "gnomad_AF_NFE"
        ])

    var_id_to_af = {}
    if gnomad_annotator is not None:
        logger.info("Annotating passing variants with gnomAD AF_nfe...")
        gnomad_annotator.load()
        batch = [(vid, info["rec"]["chrom"], info["rec"]["pos"], info["rec"]["ref"], info["rec"]["alt"])
                 for vid, info in final_pass.items()]
        results = gnomad_annotator.annotate_batch([(c, p, r, a) for _, c, p, r, a in batch])
        for vid, chrom, pos, ref, alt in batch:
            var_id_to_af[vid] = results.get(f"{chrom}:{pos}:{ref}>{alt}")

    rows = []
    for vid, info in final_pass.items():
        rec = info["rec"]
        rows.append({
            "#CHROM": rec["chrom"],
            "POS": rec["pos"],
            "ID": vid,
            "REF": rec["ref"],
            "ALT": rec["alt"],
            "gene": info["gene"],
            "zygosity": rec["zygosity"],
            "reason": info["reason"],
            "evidence_caller": rec["caller"],
            "callers_detected": ",".join(rec["callers_all"]),
            "GT": rec["gt"],
            "PS": rec["ps"],
            "samples": ",".join(sorted(info["samples"])),
            "N_samples": len(info["samples"]),
            "gnomad_AF_NFE": var_id_to_af.get(vid),
        })

    final_df = pd.DataFrame(rows).sort_values(["#CHROM", "POS"]).reset_index(drop=True)

    N_homo = (final_df["zygosity"] == "homo").sum()
    N_het = (final_df["zygosity"] == "het").sum()
    logger.info(f"\n{'='*60}")
    logger.info(f"FINAL RESULTS")
    logger.info(f"{'='*60}")
    logger.info(f"Total variants: {len(final_df)} (homo: {N_homo}, het/compound-het-trans: {N_het})")
    logger.info(f"Unique genes: {final_df['gene'].nunique()}")
    logger.info(f"{'='*60}\n")

    return final_df


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
    return parser.parse_args()


def main():
    args = get_args()

    logger.info("Starting recessive-modifier search over consolidated+snpEff VCFs")

    discordant_groups, concordant_severe = load_family_groups(args.family_json)
    manifest = load_vcf_manifest(args.vcf_manifest)

    gnomad_annotator = None
    if args.gnomad_chm13_vcf:
        gnomad_annotator = GnomadNFEAnnotator(args.gnomad_chm13_vcf)

    res = run(discordant_groups, concordant_severe, manifest, args.sv_merge_dist, gnomad_annotator)

    output_file = args.output or f"{timestamp}_recessive_modifier_consolidated.tsv"
    res.to_csv(output_file, index=False, sep="\t", header=True)
    logger.info(f"\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
