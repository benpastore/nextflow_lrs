#!/usr/bin/env python3
"""
Copy-number comparison of mild vs severe siblings, gene by gene, from each
child's Spectre CNV calls (modules/spectre -> *.spectre.vcf.gz).

1. Per child, every gene overlapped by a CNV call gets that call's copy
   number; any gene with no call is taken as CN 2 (Spectre only reports
   departures from diploid). When several calls overlap one gene, the one
   covering the most of the gene wins (n_cnvs records how many there were).
   CN comes from INFO/CN or FORMAT/CN when Spectre writes it, else from
   SVTYPE + GT (DEL het = 1, DEL hom = 0, DUP het = 3, DUP hom = 4) -- the
   cn_source column says which.

2. Within each family (family.json, same grouping as the modifier search:
   every mild child vs every severe sibling of the same father+mother),
   one row per (mild, severe sibling, gene) where the CN differs
   -> <prefix>.cnv_mild_vs_severe.tsv

3. Gene summary -> <prefix>.cnv_gene_summary.tsv. Per gene:
     n_mild_higher   mild children whose CN is ABOVE every severe sibling's
     n_mild_lower    mild children whose CN is BELOW every severe sibling's
     n_mild_mixed    differs from some severe siblings but not all
     n_pairs_mild_higher / n_pairs_mild_lower   same, counted per mild-severe pair
     sign_test_p     two-sided binomial of n_mild_higher vs n_mild_lower
     n_background_gain / _loss   concordant-severe children (no mild sibling)
                     with CN > 2 / < 2 at the gene -- a gain the mild child
                     has that's also common among severe children is a weak
                     modifier candidate
   Sorted by n_mild_higher, then fewest n_mild_lower.

4. Gene x sample copy-number matrix (every gene with CN != 2 in anyone)
   -> <prefix>.cnv_gene_cn.tsv
"""
import argparse
import csv
from collections import defaultdict

from _cnv_repeat_common import (
    IntervalIndex, cohort_groups, fmt, load_gtf, load_manifest, logger, opener, sign_test_p,
)

GT_CN = {("DEL", "het"): 1, ("DEL", "hom"): 0, ("DUP", "het"): 3, ("DUP", "hom"): 4}


class CNV(object):
    __slots__ = ("chrom", "start", "end", "svtype", "cn", "cn_source", "vid")

    def __init__(self, chrom, start, end, svtype, cn, cn_source, vid):
        self.chrom, self.start, self.end = chrom, start, end
        self.svtype, self.cn, self.cn_source, self.vid = svtype, cn, cn_source, vid

    @property
    def label(self):
        return f"{self.svtype}:{self.chrom}:{self.start}-{self.end}:CN{self.cn}"


def _int_cn(value):
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return None


def _zygosity(gt):
    alleles = [a for a in gt.replace("|", "/").split("/") if a not in (".", "")]
    if not alleles or all(a == "0" for a in alleles):
        return None
    return "hom" if len(set(alleles)) == 1 else "het"


def parse_spectre_vcf(path, min_size, keep_filtered, exclude_chroms):
    cnvs, sources, skipped = [], defaultdict(int), defaultdict(int)
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, vid, _ref, alt, _qual, filt, info_s = fields[:8]
            if chrom in exclude_chroms:
                skipped["excluded_chrom"] += 1
                continue
            if not keep_filtered and filt not in ("PASS", ".", ""):
                skipped["filtered"] += 1
                continue
            info = {}
            for kv in info_s.split(";"):
                k, _, v = kv.partition("=")
                info[k] = v
            fmt_d = dict(zip(fields[8].split(":"), fields[9].split(":"))) if len(fields) > 9 else {}

            svtype = (info.get("SVTYPE") or alt.strip("<>")).upper()
            if svtype.startswith("DUP"):
                svtype = "DUP"
            start = int(pos)
            try:
                end = int(info["END"])
            except (KeyError, ValueError):
                svlen = _int_cn(info.get("SVLEN"))
                end = start + abs(svlen) if svlen else start + 1

            cn, src = _int_cn(info.get("CN")), "INFO/CN"
            if cn is None:
                cn, src = _int_cn(fmt_d.get("CN")), "FORMAT/CN"
            if cn is None:
                cn, src = GT_CN.get((svtype, _zygosity(fmt_d.get("GT", "")))), "GT"
            if cn is None:
                skipped["no_copy_number"] += 1
                continue
            if svtype not in ("DEL", "DUP"):
                svtype = "DEL" if cn < 2 else "DUP"
            if cn == 2:
                skipped["cn_2"] += 1
                continue
            if end - start < min_size:
                skipped["below_min_size"] += 1
                continue
            sources[src] += 1
            cnvs.append(CNV(chrom, start, end, svtype, cn, src, vid))
    return cnvs, dict(sources), dict(skipped)


def gene_copy_numbers(cnvs, gene_index, min_gene_overlap):
    """{gene_key: (cn, cnv, overlap_frac, n_cnvs)} for genes hit by any CNV."""
    hits = defaultdict(list)
    for c in cnvs:
        for g in gene_index.overlapping(c.chrom, c.start, c.end):
            frac = (min(g.end, c.end) - max(g.start, c.start)) / float(max(1, g.end - g.start))
            if frac >= min_gene_overlap:
                hits[g.key].append((frac, abs(c.cn - 2), c))
    out = {}
    for key, lst in hits.items():
        frac, _, best = max(lst, key=lambda x: (x[0], x[1]))
        out[key] = (best.cn, best, frac, len(lst))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-family_json", required=True)
    p.add_argument("-manifest", required=True, help="sampleID<TAB>spectre.vcf.gz, one per line")
    p.add_argument("-gtf", required=True, help="Gene model on the same reference (e.g. params.snpeff_gtf)")
    p.add_argument("-output_prefix", default="cohort")
    p.add_argument("-min_size", type=int, default=0, help="Ignore CNVs shorter than this many bp")
    p.add_argument("-min_gene_overlap", type=float, default=0.0,
                   help="Fraction of a gene a CNV must cover to change its CN (0 = any overlap; 1 = whole gene)")
    p.add_argument("-include_noncoding", action="store_true", help="Also report non-protein-coding genes")
    p.add_argument("-keep_filtered", action="store_true", help="Keep CNVs whose FILTER isn't PASS/.")
    p.add_argument("-exclude_chroms", default="chrM", help="Comma-separated chromosomes to ignore")
    args = p.parse_args()

    exclude = set(c for c in args.exclude_chroms.split(",") if c)
    manifest = load_manifest(args.manifest)
    groups, background, phenotype = cohort_groups(args.family_json, manifest)
    if not groups:
        logger.warning("no mild child with a severe sibling has Spectre output -- tables will be empty")

    _, loci = load_gtf(args.gtf, exclude)
    if not args.include_noncoding:
        loci = [g for g in loci if g.coding]
    gene_index = IntervalIndex(loci)
    gene_by_key = {g.key: g for g in loci}

    samples = sorted(phenotype)
    gene_cn = {}
    for s in samples:
        cnvs, sources, skipped = parse_spectre_vcf(manifest[s], args.min_size, args.keep_filtered, exclude)
        gene_cn[s] = gene_copy_numbers(cnvs, gene_index, args.min_gene_overlap)
        logger.info(f"{s} ({phenotype[s]}): {len(cnvs)} CNV(s) [CN from {sources or 'none'}; "
                    f"skipped {skipped or 'none'}] -> {len(gene_cn[s])} gene(s) with CN != 2")

    def cn_of(sample, key):
        hit = gene_cn[sample].get(key)
        return hit[0] if hit else 2

    affected = set()
    for s in samples:
        affected.update(gene_cn[s])
    affected = sorted(affected, key=lambda k: (gene_by_key[k].chrom, gene_by_key[k].start))

    prefix = args.output_prefix
    # ---- per-pair table
    pair_cols = ["MILD_SAMPLE", "SEVERE_SAMPLE", "GENE", "CHROM", "GENE_START", "GENE_END",
                 "MILD_CN", "SEVERE_CN", "DELTA_CN", "DIRECTION", "GROUP_DIRECTION",
                 "MILD_CNV", "MILD_GENE_OVERLAP", "SEVERE_CNV", "SEVERE_GENE_OVERLAP",
                 "N_BACKGROUND_SAME_CN", "N_BACKGROUND"]
    summary = {}
    with open(f"{prefix}.cnv_mild_vs_severe.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(pair_cols)
        for key in affected:
            g = gene_by_key[key]
            stats = summary.setdefault(key, defaultdict(int))
            stats["higher_samples"] = []
            stats["lower_samples"] = []
            deltas = []
            bg_cns = [cn_of(b, key) for b in background]
            stats["bg_gain"] = sum(1 for c in bg_cns if c > 2)
            stats["bg_loss"] = sum(1 for c in bg_cns if c < 2)
            for mild, sevs in groups:
                m_cn = cn_of(mild, key)
                s_cns = [cn_of(s, key) for s in sevs]
                deltas.append(m_cn - sum(s_cns) / float(len(s_cns)))
                if m_cn > max(s_cns):
                    group_dir = "mild_higher"
                    stats["n_mild_higher"] += 1
                    stats["higher_samples"].append(mild)
                elif m_cn < min(s_cns):
                    group_dir = "mild_lower"
                    stats["n_mild_lower"] += 1
                    stats["lower_samples"].append(mild)
                elif any(c != m_cn for c in s_cns):
                    group_dir = "mixed"
                    stats["n_mild_mixed"] += 1
                else:
                    continue
                if m_cn > 2:
                    stats["n_mild_gain"] += 1
                elif m_cn < 2:
                    stats["n_mild_loss"] += 1
                m_hit = gene_cn[mild].get(key)
                for s, s_cn in zip(sevs, s_cns):
                    if s_cn == m_cn:
                        continue
                    stats["n_pairs_mild_higher" if m_cn > s_cn else "n_pairs_mild_lower"] += 1
                    s_hit = gene_cn[s].get(key)
                    w.writerow([
                        mild, s, g.gene, g.chrom, g.start, g.end,
                        m_cn, s_cn, m_cn - s_cn, "mild_higher" if m_cn > s_cn else "mild_lower", group_dir,
                        m_hit[1].label if m_hit else "", fmt(m_hit[2], 3) if m_hit else "",
                        s_hit[1].label if s_hit else "", fmt(s_hit[2], 3) if s_hit else "",
                        sum(1 for c in bg_cns if c == m_cn) if m_cn != 2 else "", len(background),
                    ])
            stats["mean_delta"] = sum(deltas) / len(deltas) if deltas else ""

    # ---- gene summary
    sum_cols = ["GENE", "CHROM", "GENE_START", "GENE_END",
                "n_mild_tested", "n_mild_higher", "n_mild_lower", "n_mild_mixed",
                "n_pairs_mild_higher", "n_pairs_mild_lower", "mean_delta_cn", "sign_test_p",
                "n_mild_gain", "n_mild_loss", "n_background_gain", "n_background_loss", "n_background",
                "mild_higher_samples", "mild_lower_samples"]
    rows = []
    for key in affected:
        st = summary[key]
        if not (st["n_mild_higher"] or st["n_mild_lower"] or st["n_mild_mixed"]):
            continue
        g = gene_by_key[key]
        rows.append([
            g.gene, g.chrom, g.start, g.end,
            len(groups), st["n_mild_higher"], st["n_mild_lower"], st["n_mild_mixed"],
            st["n_pairs_mild_higher"], st["n_pairs_mild_lower"], fmt(st["mean_delta"], 3),
            fmt(sign_test_p(st["n_mild_higher"], st["n_mild_lower"])),
            st["n_mild_gain"], st["n_mild_loss"], st["bg_gain"], st["bg_loss"], len(background),
            ",".join(st["higher_samples"]), ",".join(st["lower_samples"]),
        ])
    rows.sort(key=lambda r: (-r[5], r[6], -r[8], r[1], r[2]))
    with open(f"{prefix}.cnv_gene_summary.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(sum_cols)
        w.writerows(rows)

    # ---- gene x sample CN matrix
    with open(f"{prefix}.cnv_gene_cn.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["GENE", "CHROM", "GENE_START", "GENE_END"] + [f"{s}:{phenotype[s]}" for s in samples])
        for key in affected:
            g = gene_by_key[key]
            w.writerow([g.gene, g.chrom, g.start, g.end] + [cn_of(s, key) for s in samples])

    n_higher = sum(1 for r in rows if r[5])
    logger.info(f"{len(affected)} gene(s) with CN != 2 in any child; {len(rows)} differ within a family, "
                f"{n_higher} with mild CN above every severe sibling in >= 1 family")


if __name__ == "__main__":
    main()
