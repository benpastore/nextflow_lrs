#!/usr/bin/env python3
"""
hap-ibd output -> per mild x severe sibling pair, per parent: where the two
siblings carry the SAME parental copy (shared) or different ones
(not_shared), then IBD0/1/2.

Input is hap-ibd run on orient_family_phase.py's VCF, where haplotype 1 is
the paternal copy and haplotype 2 the maternal copy for every child. So:
  mild hap 1 IBD with severe hap 1  ->  paternal copy shared
  mild hap 2 IBD with severe hap 2  ->  maternal copy shared
  mild hap 1 IBD with severe hap 2 (or 2 with 1)  ->  cross-haplotype
      segment: shouldn't happen between full sibs unless the parents share
      that haplotype or orientation failed; counted in the summary, not used.
Everything between a chromosome's first and last marker that isn't shared
is not_shared. A not_shared stretch shorter than -bridge_cm with shared on
BOTH sides is relabelled shared (bridged=yes): a real one needs two
crossovers within a few Mb, which is rare, while hap-ibd segments are
routinely broken by genotype or switch errors.

Outputs (-output_prefix cohort), same layout the methylation comparison reads:
  cohort.sib_parent_segments.tsv  family mild severe parent chrom start end state bridged n_sites n_mismatch
  cohort.sib_ibd_segments.tsv     family mild severe chrom start end paternal maternal ibd
  cohort.sib_ibd_summary.tsv      genome fraction IBD0/1/2 per pair (expect ~25/50/25 for full sibs),
                                  mismatch rate inside shared segments (~ genotype error), cross-hap segments

n_mismatch = marker sites in the segment where the two sibs' alleles on that
copy differ (inside shared segments this approximates the error rate).
"""
import argparse
import bisect
import gzip
import json
import sys
from collections import defaultdict


def log(msg):
    sys.stderr.write(f"[hapibd_segments] {msg}\n")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_markers(path):
    """{chrom: (positions, {sample: [allele_hap1, allele_hap2] per site})}"""
    data = {}
    samples = []
    with opener(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            x = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = x[9:]
                continue
            chrom = x[0]
            if chrom not in data:
                data[chrom] = ([], {s: ([], []) for s in samples})
            pos, gts = data[chrom]
            pos.append(int(x[1]))
            for s, g in zip(samples, x[9:]):
                gts[s][0].append(g[0])
                gts[s][1].append(g[2])
    return data, samples


def read_ibd(path):
    """[(s1, h1, s2, h2, chrom, start, end)] from hap-ibd .ibd(.gz)."""
    rows = []
    with opener(path) as f:
        for line in f:
            x = line.split()
            if len(x) >= 7:
                rows.append((x[0], int(x[1]), x[2], int(x[3]), x[4], int(x[5]), int(x[6])))
    return rows


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def segments(shared, lo, hi, bridge_bp):
    """[[start, end, state, bridged]] tiling [lo, hi] from merged shared intervals."""
    segs = []
    cur = lo
    for s, e in shared:
        s, e = max(s, lo), min(e, hi)
        if e < s:
            continue
        if s > cur:
            segs.append([cur, s - 1, "not_shared", "no"])
        segs.append([s, e, "shared", "no"])
        cur = e + 1
    if cur <= hi:
        segs.append([cur, hi, "not_shared", "no"])
    # bridge short not_shared gaps flanked by shared
    for i in range(1, len(segs) - 1):
        a, b, c = segs[i - 1], segs[i], segs[i + 1]
        if b[2] == "not_shared" and a[2] == "shared" and c[2] == "shared" and b[1] - b[0] + 1 < bridge_bp:
            b[2], b[3] = "shared", "yes"
    # collapse neighbours with the same state
    out = []
    for s in segs:
        if out and out[-1][2] == s[2]:
            out[-1][1] = s[1]
            out[-1][3] = "yes" if "yes" in (out[-1][3], s[3]) else "no"
        else:
            out.append(list(s))
    return out


def count_sites(pos, a, b, start, end):
    i, j = bisect.bisect_left(pos, start), bisect.bisect_right(pos, end)
    return j - i, sum(1 for k in range(i, j) if a[k] != b[k])


def intersect(pat, mat):
    out, i, j = [], 0, 0
    while i < len(pat) and j < len(mat):
        lo, hi = max(pat[i][0], mat[j][0]), min(pat[i][1], mat[j][1])
        if lo <= hi:
            out.append((lo, hi, pat[i][2], mat[j][2]))
        if pat[i][1] < mat[j][1]:
            i += 1
        else:
            j += 1
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-family_json", required=True)
    p.add_argument("-families", required=True, nargs="+",
                   help="family IDs (<father>_<mother>); reads <fam>.oriented.vcf.gz + <fam>.hapibd.ibd.gz")
    p.add_argument("-bridge_cm", type=float, default=2.0,
                   help="not_shared gaps shorter than this between shared segments are bridged "
                        "(uniform 1 cM/Mb map, so cM = Mb)")
    p.add_argument("-output_prefix", default="cohort")
    args = p.parse_args()
    bridge_bp = int(args.bridge_cm * 1e6)

    with open(args.family_json) as f:
        fam = json.load(f)

    pre = args.output_prefix
    seg_out = open(f"{pre}.sib_parent_segments.tsv", "w")
    ibd_out = open(f"{pre}.sib_ibd_segments.tsv", "w")
    sum_out = open(f"{pre}.sib_ibd_summary.tsv", "w")
    seg_out.write("family\tmild\tsevere\tparent\tchrom\tstart\tend\tstate\tbridged\tn_sites\tn_mismatch\n")
    ibd_out.write("family\tmild\tsevere\tchrom\tstart\tend\tpaternal\tmaternal\tibd\n")
    sum_out.write("family\tmild\tsevere\tbp_scanned\tfrac_IBD0\tfrac_IBD1\tfrac_IBD2\t"
                  "paternal_shared_sites\tpaternal_shared_mismatch\tmaternal_shared_sites\tmaternal_shared_mismatch\t"
                  "n_cross_haplotype_segments\tn_bridged\n")

    for family in args.families:
        markers, samples = read_markers(f"{family}.oriented.vcf.gz")
        ibd = read_ibd(f"{family}.hapibd.ibd.gz")
        mild = [s for s in samples if fam.get(s, {}).get("phenotype") == "mild"]
        severe = [s for s in samples if fam.get(s, {}).get("phenotype") == "severe"]
        log(f"{family}: {len(ibd):,} hap-ibd segments; mild {mild}, severe {severe}")

        for m in mild:
            for sv in severe:
                shared = {"paternal": defaultdict(list), "maternal": defaultdict(list)}
                n_cross = 0
                for s1, h1, s2, h2, chrom, st, en in ibd:
                    if {s1, s2} != {m, sv}:
                        continue
                    hm, hs = (h1, h2) if s1 == m else (h2, h1)
                    if hm == hs:
                        shared["paternal" if hm == 1 else "maternal"][chrom].append((st, en))
                    else:
                        n_cross += 1
                bp = defaultdict(int)
                tot = {"paternal": [0, 0], "maternal": [0, 0]}
                n_bridged = 0
                for chrom in sorted(markers):
                    pos, gts = markers[chrom]
                    if not pos:
                        continue
                    per = {}
                    for parent, h in (("paternal", 0), ("maternal", 1)):
                        segs = segments(merge(shared[parent].get(chrom, [])), pos[0], pos[-1], bridge_bp)
                        per[parent] = segs
                        for s, e, state, bridged in segs:
                            n, mm = count_sites(pos, gts[m][h], gts[sv][h], s, e)
                            if state == "shared":
                                tot[parent][0] += n
                                tot[parent][1] += mm
                            n_bridged += bridged == "yes"
                            seg_out.write(f"{family}\t{m}\t{sv}\t{parent}\t{chrom}\t{s}\t{e}\t{state}\t{bridged}\t{n}\t{mm}\n")
                    for lo, hi, ps_, ms_ in intersect(per["paternal"], per["maternal"]):
                        k = (ps_ == "shared") + (ms_ == "shared")
                        bp[k] += hi - lo + 1
                        ibd_out.write(f"{family}\t{m}\t{sv}\t{chrom}\t{lo}\t{hi}\t{ps_}\t{ms_}\tIBD{k}\n")
                total = sum(bp.values()) or 1
                fr = [bp[k] / float(total) for k in (0, 1, 2)]
                rate = {k: (v[1] / float(v[0]) if v[0] else 0.0) for k, v in tot.items()}
                sum_out.write(f"{family}\t{m}\t{sv}\t{sum(bp.values())}\t{fr[0]:.3f}\t{fr[1]:.3f}\t{fr[2]:.3f}\t"
                              f"{tot['paternal'][0]}\t{rate['paternal']:.4f}\t{tot['maternal'][0]}\t{rate['maternal']:.4f}\t"
                              f"{n_cross}\t{n_bridged}\n")
                log(f"  {m} vs {sv}: IBD0 {fr[0]:.1%}, IBD1 {fr[1]:.1%}, IBD2 {fr[2]:.1%}; mismatch in shared "
                    f"pat {rate['paternal']:.2%} mat {rate['maternal']:.2%}; {n_cross} cross-haplotype, {n_bridged} bridged")
                if not (0.1 <= fr[2] <= 0.45 and 0.1 <= fr[0] <= 0.45):
                    log(f"  WARNING {m} vs {sv}: IBD fractions far from ~25/50/25% -- check orientation, "
                        f"hap-ibd parameters (segments broken by genotype errors?) or sample identities")
    for fh in (seg_out, ibd_out, sum_out):
        fh.close()


if __name__ == "__main__":
    main()
