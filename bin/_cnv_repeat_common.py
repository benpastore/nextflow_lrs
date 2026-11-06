#!/usr/bin/env python3
"""
Shared code for cnv_compare.py and repeat_compare.py: CHM13 gene model
(from the same GTF snpEff is built from), interval lookups, the
mild/severe family groupings (reused from _modifier_common.py), and the
small amount of statistics both need. Not a standalone entry point.

Written for Python 3.6 (the rnaseq container): no math.comb, no walrus,
no dataclasses.
"""

import bisect
import gzip
import math
import re
from collections import defaultdict

from _modifier_common import load_family_groups, logger  # noqa: F401 (re-exported)

ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def present(path):
    import os
    return bool(path) and os.path.basename(path) != "NO_FILE" and os.path.exists(path)


def load_manifest(path):
    """sampleID<TAB>file, one per line -> {sampleID: file}"""
    manifest = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                sample, fpath = line.split("\t")[:2]
                manifest[sample] = fpath
    return manifest


def cohort_groups(family_json, available):
    """load_family_groups() restricted to samples that have caller output.

    Returns (groups, background, phenotype):
      groups      [(mild, [severe siblings])], only groups where the mild child
                  and at least one severe sibling have output
      background  concordant-severe samples with output
      phenotype   {sample: "mild"/"severe"} for every sample used
    """
    discordant, concordant = load_family_groups(family_json)
    groups, phenotype = [], {}
    for mild, severes in discordant:
        if mild not in available:
            logger.warning(f"{mild}: mild child has no caller output, skipping its family comparison")
            continue
        sev = [s for s in severes if s in available]
        missing = [s for s in severes if s not in available]
        if missing:
            logger.warning(f"{mild}: severe sibling(s) {missing} have no caller output, dropped")
        if not sev:
            logger.warning(f"{mild}: no severe sibling with caller output, skipping")
            continue
        groups.append((mild, sev))
        phenotype[mild] = "mild"
        for s in sev:
            phenotype[s] = "severe"
    background = [s for s in concordant if s in available]
    for s in background:
        phenotype[s] = "severe"
    logger.info(f"{len(groups)} mild child(ren) with severe sibling(s); {len(background)} concordant-severe background sample(s)")
    return groups, background, phenotype


# ----------------------------------------------------------------- gene model

class Transcript(object):
    __slots__ = ("tid", "gene", "chrom", "strand", "start", "end", "exons", "cds_start", "cds_end")

    def __init__(self, tid, gene, chrom, strand):
        self.tid, self.gene, self.chrom, self.strand = tid, gene, chrom, strand
        self.start = self.end = None
        self.exons = []
        self.cds_start = self.cds_end = None

    @property
    def coding(self):
        return self.cds_start is not None


class GeneLocus(object):
    __slots__ = ("gene", "chrom", "strand", "start", "end", "coding", "n_transcripts")

    def __init__(self, gene, chrom, strand, start, end, coding, n_transcripts):
        self.gene, self.chrom, self.strand = gene, chrom, strand
        self.start, self.end, self.coding, self.n_transcripts = start, end, coding, n_transcripts

    @property
    def key(self):
        return (self.gene, self.chrom, self.start)


class IntervalIndex(object):
    """Overlap / nearest queries over objects with .chrom/.start/.end (0-based, half-open)."""

    def __init__(self, items):
        by_chrom = defaultdict(list)
        for it in items:
            by_chrom[it.chrom].append(it)
        self.items, self.starts, self.maxlen = {}, {}, {}
        for chrom, lst in by_chrom.items():
            lst.sort(key=lambda it: it.start)
            self.items[chrom] = lst
            self.starts[chrom] = [it.start for it in lst]
            self.maxlen[chrom] = max(it.end - it.start for it in lst)

    def overlapping(self, chrom, start, end):
        lst = self.items.get(chrom)
        if not lst:
            return []
        starts = self.starts[chrom]
        lo = bisect.bisect_left(starts, start - self.maxlen[chrom])
        hi = bisect.bisect_left(starts, end)
        return [it for it in lst[lo:hi] if it.end > start]

    def nearest(self, chrom, start, end):
        """(item, distance) of the closest item not overlapping, or (None, None)."""
        lst = self.items.get(chrom)
        if not lst:
            return None, None
        best, best_d = None, None
        for window in (10000, 100000, 1000000, 10000000, 10 ** 10):
            for it in self.overlapping(chrom, start - window, end + window):
                d = max(it.start - end, start - it.end, 0)
                if best_d is None or d < best_d:
                    best, best_d = it, d
            if best is not None:
                return best, best_d
        return None, None


def load_gtf(path, exclude_chroms=()):
    """GTF -> (transcripts, gene_loci). Only transcript/exon/CDS lines are used,
    so it works for UCSC ncbiRefSeq GTFs (no "gene" lines) as well as
    Ensembl/GENCODE-style ones. Transcripts of the same gene name on the same
    chromosome are merged into one GeneLocus when they overlap; disjoint
    copies (paralogous duplicates sharing a symbol) stay separate loci."""
    tx = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("transcript", "exon", "CDS"):
                continue
            chrom = fields[0]
            if chrom in exclude_chroms:
                continue
            attrs = dict(ATTR_RE.findall(fields[8]))
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            gene = attrs.get("gene_name") or attrs.get("gene_id") or tid
            start, end = int(fields[3]) - 1, int(fields[4])
            key = (tid, chrom)
            t = tx.get(key)
            if t is None:
                t = tx[key] = Transcript(tid, gene, chrom, fields[6])
            if fields[2] == "exon":
                t.exons.append((start, end))
            elif fields[2] == "CDS":
                t.cds_start = start if t.cds_start is None else min(t.cds_start, start)
                t.cds_end = end if t.cds_end is None else max(t.cds_end, end)
            t.start = start if t.start is None else min(t.start, start)
            t.end = end if t.end is None else max(t.end, end)

    transcripts = [t for t in tx.values() if t.start is not None]
    for t in transcripts:
        t.exons.sort()

    by_gene = defaultdict(list)
    for t in transcripts:
        by_gene[(t.gene, t.chrom)].append(t)
    loci = []
    for (gene, chrom), ts in by_gene.items():
        ts.sort(key=lambda t: t.start)
        cluster = [ts[0]]
        cluster_end = ts[0].end
        for t in ts[1:] + [None]:
            if t is not None and t.start < cluster_end:
                cluster.append(t)
                cluster_end = max(cluster_end, t.end)
                continue
            loci.append(GeneLocus(gene, chrom, cluster[0].strand, min(c.start for c in cluster), cluster_end,
                                  any(c.coding for c in cluster), len(cluster)))
            if t is not None:
                cluster, cluster_end = [t], t.end
    logger.info(f"GTF {path}: {len(transcripts):,} transcripts, {len(loci):,} gene loci "
                f"({sum(1 for g in loci if g.coding):,} protein-coding)")
    return transcripts, loci


# ordered most -> least functionally direct; used to pick one feature per locus
FEATURE_RANK = ["CDS", "5UTR", "3UTR", "ncRNA_exon", "promoter", "intron"]


def _feature_in_transcript(t, start, end, promoter_bp):
    if start < t.end and end > t.start:
        found = []
        for es, ee in t.exons:
            if es >= end:
                break
            if ee <= start:
                continue
            if not t.coding:
                found.append("ncRNA_exon")
                continue
            os_, oe = max(es, start), min(ee, end)
            if os_ < t.cds_end and oe > t.cds_start:
                found.append("CDS")
            if os_ < t.cds_start:
                found.append("5UTR" if t.strand == "+" else "3UTR")
            if oe > t.cds_end:
                found.append("3UTR" if t.strand == "+" else "5UTR")
        if found:
            return min(found, key=FEATURE_RANK.index)
        return "intron"
    if t.strand == "-":
        prom_s, prom_e = t.end, t.end + promoter_bp
    else:
        prom_s, prom_e = t.start - promoter_bp, t.start
    if start < prom_e and end > prom_s:
        return "promoter"
    return None


class FeatureAnnotator(object):
    """Gene + feature class (CDS/5UTR/3UTR/ncRNA_exon/promoter/intron/intergenic)
    for an interval, choosing the most direct feature over all transcripts
    (protein-coding transcripts win ties)."""

    def __init__(self, transcripts, loci, promoter_bp=2000):
        self.tx_index = IntervalIndex(transcripts)
        self.gene_index = IntervalIndex(loci)
        self.promoter_bp = promoter_bp

    def annotate(self, chrom, start, end):
        per_gene = {}
        for t in self.tx_index.overlapping(chrom, start - self.promoter_bp, end + self.promoter_bp):
            feat = _feature_in_transcript(t, start, end, self.promoter_bp)
            if feat is None:
                continue
            rank = (FEATURE_RANK.index(feat), 0 if t.coding else 1)
            if t.gene not in per_gene or rank < per_gene[t.gene][0]:
                per_gene[t.gene] = (rank, feat, t.coding)
        if per_gene:
            ordered = sorted(per_gene.items(), key=lambda kv: kv[1][0])
            gene, (rank, feat, coding) = ordered[0]
            return {
                "gene": gene, "feature": feat, "gene_coding": "yes" if coding else "no",
                "all_genes": ";".join(f"{g}:{v[1]}" for g, v in ordered),
                "nearest_gene": gene, "nearest_gene_distance": 0,
            }
        near, dist = self.gene_index.nearest(chrom, start, end)
        return {
            "gene": "", "feature": "intergenic", "gene_coding": "",
            "all_genes": "",
            "nearest_gene": near.gene if near else "",
            "nearest_gene_distance": dist if near else "",
        }


# ----------------------------------------------------------------- statistics

def _log_comb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def sign_test_p(n_pos, n_neg):
    """Two-sided exact binomial (p=0.5) test of n_pos vs n_neg; "" if both 0.
    Summed in log space -- 2**n overflows a float once n > ~1000 (genome-wide
    repeat event counts get there easily)."""
    n = n_pos + n_neg
    if n == 0:
        return ""
    k = min(n_pos, n_neg)
    log_terms = [_log_comb(n, i) - n * math.log(2) for i in range(k + 1)]
    top = max(log_terms)
    log_tail = top + math.log(sum(math.exp(t - top) for t in log_terms))
    return min(1.0, 2 * math.exp(log_tail))


def mann_whitney_p(x, y):
    """Two-sided Mann-Whitney U, normal approximation with tie and continuity
    correction. "" when either group has < 2 values. With the small groups
    of a sibling cohort this is a screening statistic, not a final test."""
    n1, n2 = len(x), len(y)
    if n1 < 2 or n2 < 2:
        return ""
    combined = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    n = n1 + n2
    ranks = [0.0] * n
    tie_term = 0.0
    i = 0
    while i < n:
        j = i
        while j + 1 < n and combined[j + 1][0] == combined[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[k] = avg
        t = j - i + 1
        tie_term += t ** 3 - t
        i = j + 1
    r1 = sum(r for r, (_, g) in zip(ranks, combined) if g == 0)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    var = n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1)))
    if var <= 0:
        return 1.0
    diff = u1 - mu
    z = (abs(diff) - 0.5) / math.sqrt(var) if diff != 0 else 0.0
    return min(1.0, math.erfc(max(z, 0.0) / math.sqrt(2)))


def bh_qvalues(pvalues):
    """Benjamini-Hochberg q-values; "" entries stay "" and are not counted."""
    idx = [i for i, p in enumerate(pvalues) if p != ""]
    m = len(idx)
    q = [""] * len(pvalues)
    if not m:
        return q
    order = sorted(idx, key=lambda i: pvalues[i], reverse=True)
    prev = 1.0
    for rank, i in enumerate(order):
        r = m - rank
        prev = min(prev, pvalues[i] * m / r)
        q[i] = prev
    return q


def fmt(v, digits=4):
    if v == "" or v is None:
        return ""
    if isinstance(v, float):
        if v != 0 and (abs(v) < 10 ** -digits):
            return f"{v:.2e}"
        return f"{v:.{digits}g}" if abs(v) >= 1000 else f"{round(v, digits)}"
    return str(v)
