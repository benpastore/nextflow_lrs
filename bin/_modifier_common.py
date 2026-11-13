#!/usr/bin/env python3
"""
Shared helpers for the pipeline's Python steps (bin/): logging + timing, a
gz-aware opener, the gnomAD / 1000 Genomes AF lookup (GnomadNFEAnnotator,
used by annotate_af.py), snpEff ANN / INFO / GT parsing, the caller
priority used to pick whose GT a site is read from, and the family.json
mild/severe groupings + VCF manifests. Not a standalone entry point.
"""

import gzip
import json
import logging
import time
from collections import defaultdict

try:
    import colorlog
except ImportError:
    colorlog = None

# Configure logging -- colorlog is a nice-to-have (not every environment
# this runs in has it, e.g. the rnaseq container); fall back to plain
# logging rather than failing outright.
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

    `field` picks the INFO allele-frequency key, so the same class also reads
    the 1000 Genomes CHM13 sites file (field="AF", see
    modules/kg_af/main.nf); `label` just names the source in the logs.
    """

    def __init__(self, gnomad_vcf_path=None, field="AF_nfe", label="gnomAD"):
        self.gnomad_vcf_path = gnomad_vcf_path
        self.field = field
        self.label = label

    def load(self):
        """No-op, kept for callers: the lookup is targeted, done in annotate_batch()."""
        if not self.gnomad_vcf_path:
            logger.warning(f"No {self.label} CHM13 VCF provided")

    @timing
    def annotate_batch(self, variants):
        """
        variants: list of (chrom, pos, ref, alt) -> dict variant_id -> AF (self.field) or None

        Single streaming pass over the gnomAD VCF that only keeps records
        matching the given (already-filtered) candidates -- the full
        resource is hundreds of millions of sites, far too large to hold
        in memory as a dict, while the candidate set is tiny.
        """
        results = {f"{c}:{p}:{r}>{a}": None for c, p, r, a in variants}
        if not self.gnomad_vcf_path or not results:
            return results

        # cheap (chrom, pos) prefilter so INFO is only parsed for candidate sites
        key = f"{self.field}="
        wanted_sites = {(c, str(p)) for c, p, _, _ in variants}

        logger.info(f"Scanning {self.gnomad_vcf_path} for {len(results):,} candidate variants...")
        n_records = 0
        n_found = 0
        n_malformed = 0
        with opener(self.gnomad_vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.split("\t", 8)
                # blank / truncated lines (e.g. a trailing newline, or a
                # partially-written download) -- skip rather than crash
                if len(fields) < 8:
                    n_malformed += 1
                    continue
                n_records += 1

                chrom, pos = fields[0], fields[1]
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"
                if (chrom, pos) not in wanted_sites:
                    continue

                ref, alt, info = fields[3], fields[4], fields[7]

                af_nfe = None
                for entry in info.split(';'):
                    if entry.startswith(key):
                        af_nfe = entry.split('=', 1)[1]
                        break

                if af_nfe is None:
                    continue

                alts = alt.split(',')
                afs = af_nfe.split(',')
                if len(alts) != len(afs):
                    continue

                for a, af in zip(alts, afs):
                    variant_id = f"{chrom}:{pos}:{ref}>{a}"
                    if variant_id in results and results[variant_id] is None:
                        try:
                            results[variant_id] = float(af)
                            n_found += 1
                        except ValueError:
                            continue

        logger.info(f"{self.label} scan done: {n_records:,} records scanned, "
                    f"{self.field} found for {n_found:,}/{len(results):,} candidates")
        if n_malformed:
            logger.warning(f"Skipped {n_malformed:,} malformed (<8 column) line(s) in "
                           f"{self.gnomad_vcf_path} -- check the file isn't truncated/corrupt")
        if n_records and not n_found:
            logger.warning(f"No candidate matched {self.label} -- check that the VCF's INFO field is "
                           f"named exactly '{self.field}' and its contig names match (chr1 / 1)")
        return results


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
