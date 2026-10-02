#!/usr/bin/env python3
"""
Shared parsing/annotation/family-grouping code for the consolidated+snpEff
modifier-search scripts (recessive_modifier_consolidated.py,
dominant_modifier_consolidated.py). Not a standalone entry point.

Both scripts work from the same input (one consolidated+snpEff VCF per
sample -- CALLERS/NUM_CALLERS + ANN, with each detecting caller's own
GT:PS kept in its own sample column -- see bin/consolidate_variants.py and
modules/snpeff/main.nf) and the same family.json-derived mild/severe
sibling groupings. What differs between the two models is only how many
qualifying variants a gene needs in one mild sample (recessive: homozygous,
or 2 het confirmed in trans; dominant: 1 variant, any zygosity) -- that
logic lives in each model's own script, not here.
"""

import bisect
import gzip
import json
import logging
import os
import pickle
import time
from collections import defaultdict

try:
    import colorlog
except ImportError:
    colorlog = None

# Configure logging -- colorlog is a nice-to-have (not every environment
# this runs in has it, e.g. the rnaseq container); fall back to plain
# logging rather than failing outright. Shared logger name so recessive
# and dominant runs interleave sensibly if ever run in the same process.
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
    logger.setLevel(logging.INFO)


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

# ANN field index 7 (Transcript_BioType, per SnpEff's standard ANN spec:
# Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|
# Feature_ID|Transcript_BioType|...). Restrict to protein-coding genes --
# SnpEff sets this structurally from the presence of CDS lines in the
# build GTF, so it's populated even for a GTF (like a UCSC RefSeq one)
# with no explicit gene_biotype/transcript_biotype attribute.
PROTEIN_CODING_BIOTYPE = 'protein_coding'

# Caller preference when a site was called by more than one -- read-based
# callers (Clair3/Sniffles) first since a PS-linked phase block is more
# directly comparable across variants than an assembly-wide haplotype;
# used to pick whose GT/PS a given variant's zygosity/phase is trusted from.
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


@timing
def parse_consolidated_vcf(path):
    """
    Parse one sample's consolidated+snpEff VCF into {var_id: record}.

    record: chrom, pos, ref, alt, svtype, genes (set), caller (the one
    whose GT/PS we're trusting), callers_all (list), gt, ps, zygosity.
    """
    results = {}
    caller_col_idx = {}
    kept = 0
    filtered_ann_type = 0
    filtered_not_protein_coding = 0
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
            filtered_ann = [a for a in filtered_ann if len(a) > 7 and a[7] == PROTEIN_CODING_BIOTYPE]
            if not filtered_ann:
                filtered_not_protein_coding += 1
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

    # protein-coding count logged at INFO (not DEBUG, which this module
    # never enables) since a GTF with no usable biotype info would silently
    # zero out every variant here otherwise -- this is the one counter that
    # needs to be visible by default to catch that failure mode.
    logger.info(
        f"Parsed {os.path.basename(path)}: kept {kept}, "
        f"filtered {filtered_ann_type} (annotation type), "
        f"{filtered_not_protein_coding} (not protein-coding), "
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


def load_family_groups(family_json_path):
    """
    Parse family.json ({child: {father, mother, phenotype}}) into:
      discordant_groups: list of (mild_sampleID, [severe_sampleIDs]) -- one
                          entry per mild child, paired against every severe
                          child sharing its (father, mother)
      concordant_severe: list of severe sampleIDs whose family has no mild
                          child (background filtering cohort) -- also where
                          a child with no father/mother in family.json
                          lands if it's severe (see below)

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
