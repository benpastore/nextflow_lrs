#!/usr/bin/env python3
"""
Repeat-expansion comparison of mild vs severe siblings from each child's
Straglr genotypes (modules/straglr -> *.straglr.bed: one row per locus,
up to --max_num_clusters alleles with size / copy_number / support).

Loci are matched across children by position (within -merge_dist bp) AND
canonical motif (smallest rotation of the motif or its reverse complement,
reduced to its primitive period -- so CAG, AGC, CTG and CAGCAG are one
motif). Each child's value at a locus is its LONGEST allele, in repeat
copies. A child with no Straglr call at a locus is taken as the reference
length by default (genome-wide Straglr only reports loci with an
insertion; with a -loci catalog an absence means no coverage instead --
pass -absent_as_missing then).

Gene association: every locus is annotated from the GTF with the single
most direct feature over all transcripts (CDS > 5'UTR > 3'UTR > ncRNA exon
> promoter [-promoter_bp upstream of a TSS] > intron), plus every gene it
touches (all_genes) and, for intergenic loci, the nearest gene and its
distance. Known disease loci (e.g. a STRchive BED) can be named with
-known_loci_bed (columns: chrom start end name ...).

Two loci values "differ" when |delta| >= -min_copy_diff copies AND
|delta| / smaller >= -min_frac_diff. "Expanded" = longest allele that much
above the reference.

Outputs (<prefix>.repeat_*.tsv):
  loci           cohort catalog: locus, annotation, every child's longest allele
  mild_vs_severe one row per (mild, severe sibling, locus) that differ
  locus_summary  per locus:
                   n_mild_longer / n_mild_shorter  mild children longer / shorter
                       than EVERY severe sibling; sign_test_p on the two
                   n_mild_specific_expanded  mild children expanded where none
                       of their severe siblings are -- >= 2 = recurrent
                       mild-specific expansion, the main "pattern in the mild
                       group" signal
                   mild vs severe group medians / expanded counts, and a
                       Mann-Whitney p (all mild vs all severe incl. background)
                       with BH q -- a screen at sibling-cohort sizes
  gene_summary   locus_summary rolled up per gene (genic/promoter loci)
  pattern_summary  global patterns: within-family mild_longer vs mild_shorter
                   counts over all loci, by feature class, motif period and
                   motif, plus per-child expansion burden (mild vs siblings)
  sample_burden  per child: loci called, loci expanded, extra copies/bp,
                   expanded loci by feature class
"""
import argparse
import csv
import re
import statistics
from collections import defaultdict

from _cnv_repeat_common import (
    FEATURE_RANK, FeatureAnnotator, IntervalIndex, bh_qvalues, cohort_groups, fmt,
    load_gtf, load_manifest, logger, mann_whitney_p, opener, present, sign_test_p,
)

ALLELE_COL_RE = re.compile(r"allele(\d+):(size|copy_number|support)$")
MOTIF_COLS = ("repeat_unit", "motif", "target_repeat", "repeat")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
MISSING = ("", "-", ".", "NA", "nan", "None")
FEATURE_CLASSES = FEATURE_RANK + ["intergenic"]


def primitive(motif):
    n = len(motif)
    for p in range(1, n + 1):
        if n % p == 0 and motif[:p] * (n // p) == motif:
            return motif[:p]
    return motif


def canonical_motif(motif):
    m = primitive(motif.upper())
    if not m or any(c not in "ACGTN" for c in m):
        return m
    rc = m.translate(COMPLEMENT)[::-1]
    return min(s[i:] + s[:i] for s in (m, rc) for i in range(len(s)))


def period_class(period):
    if period <= 6:
        return str(period)
    return "7-12" if period <= 12 else "13+"


def _num(v):
    if v is None or v.strip() in MISSING:
        return None
    try:
        return float(v)
    except ValueError:
        return None


class Call(object):
    __slots__ = ("chrom", "start", "end", "motif", "period", "alleles", "max_cn", "max_size")

    def __init__(self, chrom, start, end, motif, alleles):
        self.chrom, self.start, self.end, self.motif = chrom, start, end, motif
        self.period = max(1, len(motif))
        self.alleles = sorted(alleles, key=lambda a: -a[1])  # (size, copies, support)
        self.max_size, self.max_cn = self.alleles[0][0], self.alleles[0][1]

    @property
    def allele_str(self):
        return "/".join(f"{a[1]:g}" for a in self.alleles)


def parse_straglr_bed(path, exclude_chroms):
    calls, header, n_bad = [], None, 0
    with opener(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header = line.lstrip("#").split("\t")
                continue
            fields = line.split("\t")
            if header is None:
                # straglr's documented layout when the header is missing
                header = ["chrom", "start", "end", "repeat_unit"]
                for i in range((len(fields) - 4) // 3):
                    header += [f"allele{i + 1}:size", f"allele{i + 1}:copy_number", f"allele{i + 1}:support"]
            col = {name: i for i, name in enumerate(header)}
            chrom = fields[0]
            if chrom in exclude_chroms:
                continue
            motif_i = next((col[c] for c in MOTIF_COLS if c in col), 3)
            alleles = defaultdict(dict)
            for name, i in col.items():
                m = ALLELE_COL_RE.match(name)
                if m and i < len(fields):
                    alleles[int(m.group(1))][m.group(2)] = _num(fields[i])
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                n_bad += 1
                continue
            motif = canonical_motif(fields[motif_i]) if motif_i < len(fields) else ""
            period = max(1, len(motif))
            parsed = []
            for _, a in sorted(alleles.items()):
                size, cn = a.get("size"), a.get("copy_number")
                if cn is None and size is not None:
                    cn = size / period
                if size is None and cn is not None:
                    size = cn * period
                if cn is not None:
                    parsed.append((size, cn, a.get("support")))
            if not motif or not parsed:
                n_bad += 1
                continue
            calls.append(Call(chrom, start, end, motif, parsed))
    return calls, n_bad


class Locus(object):
    __slots__ = ("lid", "chrom", "start", "end", "motif", "period", "ref_cn", "calls", "ann", "known")

    def __init__(self, chrom, motif, members):
        self.chrom, self.motif = chrom, motif
        self.period = max(1, len(motif))
        self.start = min(c.start for _, c in members)
        self.end = max(c.end for _, c in members)
        self.ref_cn = statistics.median([(c.end - c.start) / float(self.period) for _, c in members])
        self.calls = {}
        for sample, c in members:  # one call per sample: keep its longest
            if sample not in self.calls or c.max_cn > self.calls[sample].max_cn:
                self.calls[sample] = c
        self.lid = f"{chrom}:{self.start}-{self.end}:{motif}"
        self.ann, self.known = {}, ""


def build_loci(sample_calls, merge_dist):
    by_key = defaultdict(list)
    for sample, calls in sample_calls.items():
        for c in calls:
            by_key[(c.chrom, c.motif)].append((sample, c))
    loci = []
    for (chrom, motif), members in by_key.items():
        members.sort(key=lambda sc: sc[1].start)
        cluster, cluster_end = [], None
        for sc in members:
            if cluster and sc[1].start > cluster_end + merge_dist:
                loci.append(Locus(chrom, motif, cluster))
                cluster = []
            if not cluster:
                cluster_end = sc[1].end
            cluster.append(sc)
            cluster_end = max(cluster_end, sc[1].end)
        if cluster:
            loci.append(Locus(chrom, motif, cluster))
    loci.sort(key=lambda l: (l.chrom, l.start))
    return loci


class KnownLocus(object):
    __slots__ = ("chrom", "start", "end", "name")

    def __init__(self, chrom, start, end, name):
        self.chrom, self.start, self.end, self.name = chrom, start, end, name


def load_known(path):
    items = []
    with opener(path) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                items.append(KnownLocus(fields[0], int(fields[1]), int(fields[2]),
                                        fields[3] if len(fields) > 3 else f"{fields[0]}:{fields[1]}"))
            except (ValueError, IndexError):
                continue
    logger.info(f"{len(items)} known repeat loci from {path}")
    return IntervalIndex(items)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-family_json", required=True)
    p.add_argument("-manifest", required=True, help="sampleID<TAB>straglr.bed, one per line")
    p.add_argument("-gtf", default="", help="Gene model on the same reference (optional; no gene columns without it)")
    p.add_argument("-known_loci_bed", default="", help="Known/pathogenic repeat loci, chrom start end name (optional)")
    p.add_argument("-output_prefix", default="cohort")
    p.add_argument("-merge_dist", type=int, default=50, help="Max bp gap to treat calls as the same locus")
    p.add_argument("-min_copy_diff", type=float, default=3.0)
    p.add_argument("-min_frac_diff", type=float, default=0.1)
    p.add_argument("-promoter_bp", type=int, default=2000)
    p.add_argument("-min_motif_events", type=int, default=5,
                   help="Motifs with fewer within-family differences are pooled as 'other' in pattern_summary")
    p.add_argument("-absent_as_missing", action="store_true",
                   help="Treat a child with no call at a locus as missing, not reference length")
    p.add_argument("-exclude_chroms", default="chrM")
    args = p.parse_args()

    exclude = set(c for c in args.exclude_chroms.split(",") if c)
    manifest = load_manifest(args.manifest)
    groups, background, phenotype = cohort_groups(args.family_json, manifest)
    samples = sorted(phenotype)
    mild_samples = [m for m, _ in groups]
    severe_samples = [s for s in samples if phenotype[s] == "severe"]

    sample_calls = {}
    for s in samples:
        calls, n_bad = parse_straglr_bed(manifest[s], exclude)
        sample_calls[s] = calls
        logger.info(f"{s} ({phenotype[s]}): {len(calls):,} Straglr locus call(s)"
                    f"{f', {n_bad} unparseable row(s) skipped' if n_bad else ''}")
    loci = build_loci(sample_calls, args.merge_dist)
    logger.info(f"{len(loci):,} cohort loci after merging (merge_dist={args.merge_dist})")

    annotator = None
    if present(args.gtf):
        transcripts, gene_loci = load_gtf(args.gtf, exclude)
        annotator = FeatureAnnotator(transcripts, gene_loci, args.promoter_bp)
    known = load_known(args.known_loci_bed) if present(args.known_loci_bed) else None
    for l in loci:
        l.ann = annotator.annotate(l.chrom, l.start, l.end) if annotator else {
            "gene": "", "feature": "", "gene_coding": "", "all_genes": "", "nearest_gene": "", "nearest_gene_distance": ""}
        if known:
            l.known = ";".join(k.name for k in known.overlapping(l.chrom, l.start - args.merge_dist, l.end + args.merge_dist))

    def value(l, s):
        c = l.calls.get(s)
        if c is not None:
            return c.max_cn
        return None if args.absent_as_missing else l.ref_cn

    def differs(a, b):
        d = abs(a - b)
        return d >= args.min_copy_diff and d / max(min(a, b), 1.0) >= args.min_frac_diff

    def expanded(l, v):
        return v is not None and v > l.ref_cn and differs(v, l.ref_cn)

    prefix = args.output_prefix
    ann_cols = ["gene", "feature", "gene_coding", "all_genes", "nearest_gene", "nearest_gene_distance"]

    # ---- cohort catalog
    with open(f"{prefix}.repeat_loci.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["LOCUS_ID", "CHROM", "START", "END", "MOTIF", "PERIOD", "REF_COPIES"] + ann_cols
                   + ["known_locus"] + [f"{s}:{phenotype[s]}" for s in samples])
        for l in loci:
            w.writerow([l.lid, l.chrom, l.start, l.end, l.motif, l.period, fmt(l.ref_cn, 3)]
                       + [l.ann[c] for c in ann_cols] + [l.known]
                       + [l.calls[s].allele_str if s in l.calls else "" for s in samples])

    # ---- within-family comparisons
    events = []  # (locus, mild, direction) group-level, for pattern summary
    loc_stats = {}
    with open(f"{prefix}.repeat_mild_vs_severe.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["MILD_SAMPLE", "SEVERE_SAMPLE", "LOCUS_ID", "CHROM", "START", "END", "MOTIF", "REF_COPIES",
                    "MILD_ALLELES", "MILD_MAX_COPIES", "SEVERE_ALLELES", "SEVERE_MAX_COPIES",
                    "DELTA_COPIES", "DELTA_BP", "DIRECTION", "GROUP_DIRECTION", "MILD_EXPANDED", "SEVERE_EXPANDED"]
                   + ann_cols + ["known_locus"])
        for l in loci:
            st = loc_stats[l.lid] = defaultdict(int)
            st["longer"], st["specific"] = [], []
            for mild, sevs in groups:
                if mild not in l.calls and not any(s in l.calls for s in sevs):
                    continue
                m = value(l, mild)
                pairs = [(s, value(l, s)) for s in sevs]
                pairs = [(s, v) for s, v in pairs if v is not None]
                if m is None or not pairs:
                    continue
                if expanded(l, m) and not any(expanded(l, v) for _, v in pairs):
                    st["n_mild_specific_expanded"] += 1
                    st["specific"].append(mild)
                longer = [differs(m, v) and m > v for _, v in pairs]
                shorter = [differs(m, v) and m < v for _, v in pairs]
                if all(longer):
                    group_dir = "mild_longer"
                elif all(shorter):
                    group_dir = "mild_shorter"
                elif any(longer) or any(shorter):
                    group_dir = "mixed"
                else:
                    continue
                st["n_" + group_dir] += 1
                if group_dir == "mild_longer":
                    st["longer"].append(mild)
                if group_dir != "mixed":
                    events.append((l, mild, group_dir))
                for (s, v), is_long, is_short in zip(pairs, longer, shorter):
                    if not (is_long or is_short):
                        continue
                    mc, sc = l.calls.get(mild), l.calls.get(s)
                    w.writerow([mild, s, l.lid, l.chrom, l.start, l.end, l.motif, fmt(l.ref_cn, 3),
                                mc.allele_str if mc else "ref", fmt(m, 4), sc.allele_str if sc else "ref", fmt(v, 4),
                                fmt(m - v, 4), fmt((m - v) * l.period, 4), "mild_longer" if is_long else "mild_shorter",
                                group_dir, "yes" if expanded(l, m) else "no", "yes" if expanded(l, v) else "no"]
                               + [l.ann[c] for c in ann_cols] + [l.known])

    # ---- locus summary (every locus with a difference or an expansion)
    summary_rows, pvals = [], []
    for l in loci:
        st = loc_stats[l.lid]
        mv = [v for v in (value(l, s) for s in mild_samples) if v is not None]
        sv = [v for v in (value(l, s) for s in severe_samples) if v is not None]
        n_mild_exp = sum(1 for v in mv if expanded(l, v))
        n_sev_exp = sum(1 for v in sv if expanded(l, v))
        if not (st["n_mild_longer"] or st["n_mild_shorter"] or st["n_mixed"] or n_mild_exp or n_sev_exp):
            continue
        pvals.append(mann_whitney_p(mv, sv))
        summary_rows.append((l, st, mv, sv, n_mild_exp, n_sev_exp))
    qvals = bh_qvalues(pvals)
    sum_cols = ["LOCUS_ID", "CHROM", "START", "END", "MOTIF", "PERIOD", "REF_COPIES"] + ann_cols + [
        "known_locus", "n_mild_tested", "n_mild_longer", "n_mild_shorter", "n_mixed", "sign_test_p",
        "n_mild_specific_expanded", "mild_median_copies", "severe_median_copies", "median_diff_copies",
        "n_mild_expanded", "n_mild", "n_severe_expanded", "n_severe", "n_background_expanded",
        "mann_whitney_p", "mann_whitney_q", "mild_longer_samples", "mild_specific_samples"]
    out_rows = []
    for (l, st, mv, sv, n_mild_exp, n_sev_exp), pv, qv in zip(summary_rows, pvals, qvals):
        m_med = statistics.median(mv) if mv else ""
        s_med = statistics.median(sv) if sv else ""
        out_rows.append([
            l.lid, l.chrom, l.start, l.end, l.motif, l.period, fmt(l.ref_cn, 3)] + [l.ann[c] for c in ann_cols] + [
            l.known, len(groups), st["n_mild_longer"], st["n_mild_shorter"], st["n_mixed"],
            fmt(sign_test_p(st["n_mild_longer"], st["n_mild_shorter"])),
            st["n_mild_specific_expanded"], fmt(m_med, 4), fmt(s_med, 4),
            fmt(m_med - s_med, 4) if mv and sv else "",
            n_mild_exp, len(mv), n_sev_exp, len(sv),
            sum(1 for b in background if expanded(l, value(l, b))),
            fmt(pv), fmt(qv), ",".join(st["longer"]), ",".join(st["specific"])])
    i_spec, i_long, i_short = sum_cols.index("n_mild_specific_expanded"), sum_cols.index("n_mild_longer"), sum_cols.index("n_mild_shorter")
    out_rows.sort(key=lambda r: (-r[i_spec], -(r[i_long] - r[i_short]), r[1], r[2]))
    with open(f"{prefix}.repeat_locus_summary.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(sum_cols)
        w.writerows(out_rows)

    # ---- gene summary
    genes = {}
    for l, st, _mv, _sv, n_mild_exp, n_sev_exp in summary_rows:
        if not l.ann.get("gene"):
            continue
        g = genes.setdefault(l.ann["gene"], {"loci": [], "features": set(), "longer": 0, "shorter": 0,
                                             "specific": 0, "longer_samples": set(), "specific_samples": set(),
                                             "mild_exp": 0, "sev_exp": 0, "known": set()})
        g["loci"].append(f"{l.lid}({l.ann['feature']})")
        g["features"].add(l.ann["feature"])
        g["longer"] += st["n_mild_longer"]
        g["shorter"] += st["n_mild_shorter"]
        g["specific"] += st["n_mild_specific_expanded"]
        g["longer_samples"].update(st["longer"])
        g["specific_samples"].update(st["specific"])
        g["mild_exp"] += n_mild_exp
        g["sev_exp"] += n_sev_exp
        if l.known:
            g["known"].add(l.known)
    gene_rows = []
    for name, g in genes.items():
        feats = sorted(g["features"], key=lambda f: FEATURE_CLASSES.index(f) if f in FEATURE_CLASSES else 99)
        gene_rows.append([name, len(g["loci"]), ",".join(feats), g["longer"], g["shorter"],
                          fmt(sign_test_p(g["longer"], g["shorter"])), len(g["longer_samples"]),
                          g["specific"], len(g["specific_samples"]), g["mild_exp"], g["sev_exp"],
                          ";".join(sorted(g["known"])), ",".join(sorted(g["longer_samples"])),
                          ",".join(sorted(g["specific_samples"])), ";".join(g["loci"][:20])])
    gene_rows.sort(key=lambda r: (-r[8], -r[6], -(r[3] - r[4]), r[0]))
    with open(f"{prefix}.repeat_gene_summary.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["GENE", "n_loci", "features", "n_mild_longer_events", "n_mild_shorter_events", "sign_test_p",
                    "n_mild_children_longer", "n_mild_specific_expanded_events", "n_mild_children_specific_expanded",
                    "n_mild_expanded", "n_severe_expanded", "known_locus", "mild_longer_samples",
                    "mild_specific_samples", "loci"])
        w.writerows(gene_rows)

    # ---- per-sample burden
    burden = {}
    for s in samples:
        b = burden[s] = defaultdict(float)
        for l in loci:
            v = value(l, s)
            if s in l.calls:
                b["n_called"] += 1
            if expanded(l, v):
                b["n_expanded"] += 1
                b["expanded_copies"] += v - l.ref_cn
                b["expanded_bp"] += (v - l.ref_cn) * l.period
                b["exp_" + (l.ann.get("feature") or "unannotated")] += 1
    feat_cols = FEATURE_CLASSES if annotator else []
    role = {}
    for m, sevs in groups:
        role[m] = "mild"
        for s in sevs:
            role.setdefault(s, "severe_sibling")
    for s in background:
        role[s] = "severe_background"
    with open(f"{prefix}.repeat_sample_burden.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["SAMPLE", "PHENOTYPE", "ROLE", "n_loci_called", "n_expanded", "expanded_copies", "expanded_bp"]
                   + [f"n_expanded_{f}" for f in feat_cols])
        for s in samples:
            b = burden[s]
            w.writerow([s, phenotype[s], role.get(s, ""), int(b["n_called"]), int(b["n_expanded"]),
                        fmt(b["expanded_copies"], 4), fmt(b["expanded_bp"], 4)]
                       + [int(b["exp_" + f]) for f in feat_cols])

    # ---- global patterns
    def tally(rows_key):
        acc = defaultdict(lambda: [set(), 0, 0])
        for l, mild, d in events:
            for k in rows_key(l):
                a = acc[k]
                a[0].add(l.lid)
                a[1 if d == "mild_longer" else 2] += 1
        return acc

    pattern = []

    def add(kind, acc, order=None):
        keys = order if order is not None else sorted(acc, key=lambda k: -(acc[k][1] + acc[k][2]))
        for k in keys:
            if k not in acc:
                continue
            loci_set, n_long, n_short = acc[k]
            n = n_long + n_short
            pattern.append([kind, k, len(loci_set), n_long, n_short, fmt(n_long / float(n), 3) if n else "",
                            fmt(sign_test_p(n_long, n_short))])

    add("all_loci", tally(lambda l: ["all"]))
    if annotator:
        add("feature", tally(lambda l: [l.ann["feature"]]), FEATURE_CLASSES)
        add("genic", tally(lambda l: ["coding_gene" if l.ann["gene_coding"] == "yes" else
                                      "noncoding_gene" if l.ann["gene"] else "intergenic"]))
    add("motif_period", tally(lambda l: [period_class(l.period)]),
        [str(i) for i in range(1, 7)] + ["7-12", "13+"])
    motif_acc = tally(lambda l: [l.motif])
    pooled = defaultdict(lambda: [set(), 0, 0])
    for k, a in list(motif_acc.items()):
        if a[1] + a[2] < args.min_motif_events:
            pooled["other"][0].update(a[0])
            pooled["other"][1] += a[1]
            pooled["other"][2] += a[2]
            del motif_acc[k]
    motif_acc.update(pooled)
    add("motif", motif_acc)
    if known:
        add("known_locus", tally(lambda l: [l.known] if l.known else []))

    metrics = ["n_expanded", "expanded_copies", "expanded_bp"] + [f"exp_{f}" for f in feat_cols]
    for metric in metrics:
        hi = lo = 0
        for mild, sevs in groups:
            sib = statistics.mean([burden[s][metric] for s in sevs])
            if burden[mild][metric] > sib:
                hi += 1
            elif burden[mild][metric] < sib:
                lo += 1
        n = hi + lo
        pattern.append(["burden", metric.replace("exp_", "n_expanded_"), "", hi, lo,
                        fmt(hi / float(n), 3) if n else "", fmt(sign_test_p(hi, lo))])

    with open(f"{prefix}.repeat_pattern_summary.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["CATEGORY_TYPE", "CATEGORY", "n_loci", "n_mild_longer", "n_mild_shorter",
                    "frac_mild_longer", "sign_test_p"])
        w.writerows(pattern)

    n_spec = sum(1 for r in out_rows if r[i_spec] >= 2)
    logger.info(f"{len(out_rows):,} loci differ or are expanded; {len(events):,} within-family mild_longer/"
                f"mild_shorter event(s); {n_spec} locus/loci with recurrent (>= 2) mild-specific expansion")


if __name__ == "__main__":
    main()
