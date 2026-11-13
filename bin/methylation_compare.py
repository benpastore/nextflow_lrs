#!/usr/bin/env python3
"""
Mild vs severe promoter methylation, per parental haplotype.

For every discordant (mild child, severe sibling) pair and every gene, the
mean CpG 5mC level of the gene's promoter(s) is compared:

  combined   all reads (HP1 + HP2 + unphased)
  paternal   only reads on the child's paternal haplotype
  maternal   only reads on the child's maternal haplotype

and each parental side is labelled with whether the two sibs inherited the
SAME haplotype from that parent there (bin/sib_haplotype_sharing.py
segments): a methylation difference on a not_shared haplotype can be
sequence-driven (cis); one on a shared haplotype is an epigenetic
difference on identical sequence.

Inputs
  -manifest           childID<TAB>hp1.bed.gz<TAB>hp2.bed.gz<TAB>ungrouped.bed.gz<TAB>phased.vcf.gz
                      (MODKIT_PILEUP --partition-tag HP output, and the
                      longphase_sv VCF whose phasing those HP tags came from)
  -parent_manifest    parentID<TAB>clair3.vcf.gz
  -family_json        {"child": {"father", "mother", "phenotype"}}
  -sib_parent_segments  cohort.sib_parent_segments.tsv
  -gtf                gene model on the reference (same as snpEff / cnv_compare)

HP1/HP2 -> paternal/maternal
  modkit's HP partition is longphase's haplotype label, which is only
  consistent WITHIN a phase block (PS). Each block's HP1 is assigned to a
  parent by voting over the child's phased het SNVs at sites where exactly
  one parent carries the ALT (that parent het or hom-alt, the other
  hom-ref): child 1|0 -> HP1 is that parent, 0|1 -> HP1 is the other.
  A block is assigned when it has >= -min_block_votes informative sites
  and >= -min_block_agree of them agree. CpGs outside an assigned block
  count only towards the combined level. Assumes the WhatsHap/longphase
  haplotag convention (HP:1 = reads carrying the first GT allele);
  cohort.methylation_imprinting_check.tsv lets you confirm it (e.g. SNRPN,
  PEG3, MEST promoters should be methylated on the MATERNAL side, H19 on
  the PATERNAL side).

Promoter = -upstream / +downstream bp around every transcript's TSS,
merged per gene. Level = sum(Nmod) / sum(Nvalid_cov) over the promoter's
CpGs ('m' rows; per-haplotype CpGs need >= -min_cov reads), reported when
>= -min_cpgs CpGs contribute.

Outputs (-output_prefix cohort)
  cohort.promoter_methylation.tsv          every pair x gene with data
  cohort.promoter_methylation.flagged.tsv  |delta| >= -min_delta on any side
  cohort.methylation_phase_blocks.tsv      per child: blocks assigned to a parent
  cohort.methylation_imprinting_check.tsv  known imprinted promoters, per side

Written for Python 3.6 (the rnaseq container).
"""
import argparse
import bisect
import csv
import json
from collections import defaultdict

from _cnv_repeat_common import load_gtf, logger, opener, present

# promoter methylated on this parent's allele (germline imprint)
IMPRINTED = {
    "SNRPN": "maternal", "PEG3": "maternal", "MEST": "maternal", "PLAGL1": "maternal",
    "KCNQ1OT1": "maternal", "PEG10": "maternal", "H19": "paternal", "MEG3": "paternal",
}
SIDES = ("paternal", "maternal")
OTHER = {"paternal": "maternal", "maternal": "paternal"}


# ------------------------------------------------------------------ promoters

class Promoter(object):
    __slots__ = ("idx", "gene", "chrom", "start", "end", "tss")

    def __init__(self, idx, gene, chrom, start, end, tss):
        self.idx, self.gene, self.chrom, self.start, self.end, self.tss = idx, gene, chrom, start, end, tss


def build_promoters(gtf, upstream, downstream, exclude_chroms):
    """One Promoter per (gene, chrom, merged window cluster). Windows are
    0-based half-open; tss = first TSS in the cluster (for sharing lookup)."""
    transcripts, _ = load_gtf(gtf, exclude_chroms)
    wins = defaultdict(list)
    for t in transcripts:
        if t.strand == "-":
            tss = t.end - 1
            s, e = tss - downstream + 1, tss + upstream + 1
        else:
            tss = t.start
            s, e = tss - upstream, tss + downstream
        wins[(t.gene, t.chrom)].append((max(0, s), e, tss))
    proms = []
    for (gene, chrom), ws in wins.items():
        ws.sort()
        cur = list(ws[0])
        for s, e, tss in ws[1:]:
            if s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                proms.append(Promoter(len(proms), gene, chrom, cur[0], cur[1], cur[2]))
                cur = [s, e, tss]
        proms.append(Promoter(len(proms), gene, chrom, cur[0], cur[1], cur[2]))
    logger.info(f"{len(proms):,} promoter windows over {len({p.gene for p in proms}):,} genes "
                f"(-{upstream}/+{downstream} bp around each TSS)")
    return proms


class PromoterIndex(object):
    """Point lookup: which promoters contain a position."""

    def __init__(self, proms):
        by_chrom = defaultdict(list)
        for p in proms:
            by_chrom[p.chrom].append(p)
        self.items, self.starts, self.maxlen = {}, {}, {}
        for chrom, lst in by_chrom.items():
            lst.sort(key=lambda p: p.start)
            self.items[chrom] = lst
            self.starts[chrom] = [p.start for p in lst]
            self.maxlen[chrom] = max(p.end - p.start for p in lst)

    def at(self, chrom, pos):
        lst = self.items.get(chrom)
        if not lst:
            return ()
        starts = self.starts[chrom]
        hi = bisect.bisect_right(starts, pos)
        lo = bisect.bisect_left(starts, pos - self.maxlen[chrom])
        return [p for p in lst[lo:hi] if p.end > pos]


# ------------------------------------------------------------------ phase blocks -> parent

def parse_gt(sample_field, fmt):
    keys = fmt.split(":")
    vals = sample_field.split(":")
    d = dict(zip(keys, vals))
    return d.get("GT", ""), d.get("PS", "")


def read_child_phased(path):
    """{(chrom,pos): (ref, alt, gt, ps)} for phased het biallelic SNVs, and
    {(chrom, ps): [min_pos, max_pos]} over every phased het record."""
    sites, spans = {}, {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10 or "," in c[4] or c[4].startswith("<"):
                continue
            gt, ps = parse_gt(c[9], c[8])
            if gt not in ("0|1", "1|0") or ps in ("", "."):
                continue
            chrom, pos = c[0], int(c[1])
            span = spans.get((chrom, ps))
            if span is None:
                spans[(chrom, ps)] = [pos, pos]
            else:
                span[0], span[1] = min(span[0], pos), max(span[1], pos)
            if len(c[3]) == 1 and len(c[4]) == 1:
                sites[(chrom, pos)] = (c[3], c[4], gt, ps)
    return sites, spans


def read_parent_classes(path, wanted):
    """{(chrom,pos): 'ref'/'alt'/'unknown'} for positions in `wanted`
    (absent from the VCF = hom-ref, as Clair3 only writes variant sites).
    'alt' = carries the child's ALT (het or hom); a different ALT, a
    non-PASS record or a missing GT = unknown."""
    out = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t", 10)
            key = (c[0], int(c[1]))
            site = wanted.get(key)
            if site is None:
                continue
            if c[6] not in ("PASS", "."):
                out[key] = "unknown"
                continue
            gt, _ = parse_gt(c[9].rstrip("\n"), c[8])
            alleles = gt.replace("|", "/").split("/")
            if "." in alleles or not gt:
                out[key] = "unknown"
                continue
            alts = c[4].split(",")
            if c[3] != site[0] or site[1] not in alts:
                out[key] = "unknown" if any(a != "0" for a in alleles) else "ref"
                continue
            idx = str(alts.index(site[1]) + 1)
            others = [a for a in alleles if a not in ("0", idx)]
            out[key] = "unknown" if others else ("alt" if idx in alleles else "ref")
    return out


class BlockMap(object):
    """Per chrom sorted (start, end, hp1_parent) for assigned blocks."""

    def __init__(self, blocks):
        by_chrom = defaultdict(list)
        for chrom, start, end, parent in blocks:
            by_chrom[chrom].append((start, end, parent))
        self.items, self.starts = {}, {}
        for chrom, lst in by_chrom.items():
            lst.sort()
            self.items[chrom] = lst
            self.starts[chrom] = [b[0] for b in lst]

    def hp1_parent(self, chrom, pos):
        lst = self.items.get(chrom)
        if not lst:
            return None
        i = bisect.bisect_right(self.starts[chrom], pos) - 1
        if i < 0:
            return None
        s, e, parent = lst[i]
        return parent if s <= pos <= e else None


def assign_blocks(child, phased_vcf, father_vcf, mother_vcf, min_votes, min_agree, block_rows):
    sites, spans = read_child_phased(phased_vcf)
    fa = read_parent_classes(father_vcf, sites)
    mo = read_parent_classes(mother_vcf, sites)
    votes = defaultdict(lambda: {"paternal": 0, "maternal": 0})
    for key, (ref, alt, gt, ps) in sites.items():
        f, m = fa.get(key, "ref"), mo.get(key, "ref")
        if f == "alt" and m == "ref":
            alt_parent = "paternal"
        elif m == "alt" and f == "ref":
            alt_parent = "maternal"
        else:
            continue
        hp1 = alt_parent if gt == "1|0" else OTHER[alt_parent]
        votes[(key[0], ps)][hp1] += 1

    assigned, bp_assigned, bp_total = [], 0, 0
    for (chrom, ps), (s, e) in spans.items():
        v = votes.get((chrom, ps), {"paternal": 0, "maternal": 0})
        n = v["paternal"] + v["maternal"]
        best = "paternal" if v["paternal"] >= v["maternal"] else "maternal"
        agree = v[best] / n if n else 0.0
        ok = n >= min_votes and agree >= min_agree
        bp_total += e - s + 1
        if ok:
            assigned.append((chrom, s, e, best))
            bp_assigned += e - s + 1
        block_rows.append({"sample": child, "chrom": chrom, "PS": ps, "start": s, "end": e,
                           "votes_hp1_paternal": v["paternal"], "votes_hp1_maternal": v["maternal"],
                           "hp1_parent": best if ok else "", "agreement": f"{agree:.3f}" if n else ""})
    logger.info(f"{child}: {len(assigned):,}/{len(spans):,} phase blocks assigned to a parent "
                f"({bp_assigned / 1e6:,.1f} / {bp_total / 1e6:,.1f} Mb)")
    return BlockMap(assigned)


# ------------------------------------------------------------------ methylation

def accumulate(child, beds, blocks, pindex, min_cov):
    """{promoter idx: {side: [nmod, nvalid, ncpg]}} for side in
    combined / paternal / maternal."""
    acc = defaultdict(lambda: {"combined": [0, 0, 0], "paternal": [0, 0, 0], "maternal": [0, 0, 0]})
    combined_cpg = defaultdict(set)
    n_rows = 0
    for hp, path in beds:
        if not present(path):
            continue
        with opener(path) as f:
            for line in f:
                c = line.split()
                if len(c) < 12 or c[3] != "m":
                    continue
                chrom, pos = c[0], int(c[1])
                proms = pindex.at(chrom, pos)
                if not proms:
                    continue
                nvalid, nmod = int(c[9]), int(c[11])
                if nvalid == 0:
                    continue
                n_rows += 1
                side = None
                if hp in ("1", "2"):
                    hp1 = blocks.hp1_parent(chrom, pos)
                    if hp1 is not None:
                        side = hp1 if hp == "1" else OTHER[hp1]
                for p in proms:
                    a = acc[p.idx]
                    a["combined"][0] += nmod
                    a["combined"][1] += nvalid
                    combined_cpg[p.idx].add(pos)
                    if side is not None and nvalid >= min_cov:
                        s = a[side]
                        s[0] += nmod
                        s[1] += nvalid
                        s[2] += 1
    for idx, cpgs in combined_cpg.items():
        acc[idx]["combined"][2] = len(cpgs)
    logger.info(f"{child}: {n_rows:,} promoter CpG rows from {sum(1 for _, p in beds if present(p))} bedMethyl file(s)")
    return acc


def level(a, min_cpgs):
    nmod, nvalid, ncpg = a
    return (nmod / nvalid, ncpg) if nvalid and ncpg >= min_cpgs else (None, ncpg)


# ------------------------------------------------------------------ sibling sharing

class Sharing(object):
    """Per (mild, severe, parent, chrom): sorted segments (as in haplotype_modifier.py)."""

    def __init__(self, path):
        self.segs = defaultdict(list)
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                self.segs[(row["mild"], row["severe"], row["parent"], row["chrom"])].append(
                    (int(row["start"]), int(row["end"]), row["state"] == "shared"))
        self.starts = {}
        for k, v in self.segs.items():
            v.sort()
            self.starts[k] = [s[0] for s in v]

    def status(self, mild, severe, parent, chrom, pos):
        k = (mild, severe, parent, chrom)
        if k not in self.segs:
            return "unknown"
        i = bisect.bisect_right(self.starts[k], pos) - 1
        if i < 0:
            return "unknown"
        lo, hi, st = self.segs[k][i]
        if not lo <= pos <= hi:
            return "unknown"
        return "shared" if st else "not_shared"


# ------------------------------------------------------------------ run

def fmt(v):
    return "" if v is None else f"{v:.3f}"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-family_json", required=True)
    p.add_argument("-manifest", required=True)
    p.add_argument("-parent_manifest", required=True)
    p.add_argument("-sib_parent_segments", required=True)
    p.add_argument("-gtf", required=True)
    p.add_argument("-upstream", type=int, default=2000)
    p.add_argument("-downstream", type=int, default=500)
    p.add_argument("-min_cov", type=int, default=3, help="Reads per CpG for the per-haplotype levels")
    p.add_argument("-min_cpgs", type=int, default=5, help="CpGs needed for a promoter level")
    p.add_argument("-min_delta", type=float, default=0.2, help="|mild - severe| flagged")
    p.add_argument("-min_block_votes", type=int, default=5)
    p.add_argument("-min_block_agree", type=float, default=0.8)
    p.add_argument("-exclude_chroms", default="chrM")
    p.add_argument("-output_prefix", default="cohort")
    args = p.parse_args()

    with open(args.family_json) as f:
        family = json.load(f)
    children = {}
    with open(args.manifest) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 5:
                children[c[0]] = {"beds": [("1", c[1]), ("2", c[2]), ("ungrouped", c[3])], "vcf": c[4]}
    parents = {}
    with open(args.parent_manifest) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 2:
                parents[c[0]] = c[1]

    pairs = []
    for mild, info in family.items():
        if info.get("phenotype") != "mild" or mild not in children:
            continue
        for sev, sinfo in family.items():
            if (sinfo.get("phenotype") == "severe" and sev in children
                    and sinfo.get("father") == info.get("father") and sinfo.get("mother") == info.get("mother")):
                pairs.append((mild, sev))
    logger.info(f"{len(pairs)} mild/severe pair(s): {pairs}")

    proms = build_promoters(args.gtf, args.upstream, args.downstream,
                            set(c for c in args.exclude_chroms.split(",") if c))
    pindex = PromoterIndex(proms)

    block_rows = []
    levels = {}   # child -> {prom idx: {side: (level, ncpg)}}
    for child in sorted({c for pr in pairs for c in pr}):
        fam = family[child]
        fv, mv = parents.get(fam.get("father")), parents.get(fam.get("mother"))
        if fv and mv:
            blocks = assign_blocks(child, children[child]["vcf"], fv, mv,
                                   args.min_block_votes, args.min_block_agree, block_rows)
        else:
            logger.warning(f"{child}: missing a parent VCF -- combined levels only")
            blocks = BlockMap([])
        acc = accumulate(child, children[child]["beds"], blocks, pindex, args.min_cov)
        levels[child] = {idx: {side: level(a[side], args.min_cpgs) for side in a} for idx, a in acc.items()}

    sharing = Sharing(args.sib_parent_segments)
    cols = ["MILD_SAMPLE", "SEVERE_SAMPLE", "GENE", "CHROM", "PROMOTER_START", "PROMOTER_END",
            "mild_meth", "severe_meth", "delta", "n_cpg_mild", "n_cpg_severe",
            "pat_mild", "pat_severe", "pat_delta", "pat_haplotype",
            "mat_mild", "mat_severe", "mat_delta", "mat_haplotype", "meth_flag"]
    rows, flagged = [], []
    for mild, sev in pairs:
        lm, ls = levels[mild], levels[sev]
        for idx in sorted(set(lm) & set(ls)):
            pr = proms[idx]
            (m, ncm), (s, ncs) = lm[idx]["combined"], ls[idx]["combined"]
            row = {"MILD_SAMPLE": mild, "SEVERE_SAMPLE": sev, "GENE": pr.gene, "CHROM": pr.chrom,
                   "PROMOTER_START": pr.start, "PROMOTER_END": pr.end,
                   "mild_meth": fmt(m), "severe_meth": fmt(s), "n_cpg_mild": ncm, "n_cpg_severe": ncs,
                   "delta": fmt(m - s) if m is not None and s is not None else ""}
            flags = []
            if m is not None and s is not None and abs(m - s) >= args.min_delta:
                flags.append(f"combined:mild_{'higher' if m > s else 'lower'}({m - s:+.2f})")
            any_side = False
            for side, pre in (("paternal", "pat"), ("maternal", "mat")):
                hm, hs = lm[idx][side][0], ls[idx][side][0]
                hap = sharing.status(mild, sev, side, pr.chrom, pr.tss)
                row.update({f"{pre}_mild": fmt(hm), f"{pre}_severe": fmt(hs), f"{pre}_haplotype": hap,
                            f"{pre}_delta": fmt(hm - hs) if hm is not None and hs is not None else ""})
                if hm is not None and hs is not None:
                    any_side = True
                    if abs(hm - hs) >= args.min_delta:
                        flags.append(f"{side}:mild_{'higher' if hm > hs else 'lower'}({hm - hs:+.2f},{hap})")
            if m is None and s is None and not any_side:
                continue
            row["meth_flag"] = ";".join(flags)
            rows.append(row)
            if flags:
                flagged.append(row)

    pre = args.output_prefix
    for path, data in ((f"{pre}.promoter_methylation.tsv", rows), (f"{pre}.promoter_methylation.flagged.tsv", flagged)):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(data)
        logger.info(f"wrote {len(data):,} row(s) -> {path}")

    with open(f"{pre}.methylation_phase_blocks.tsv", "w", newline="") as f:
        bcols = ["sample", "chrom", "PS", "start", "end", "votes_hp1_paternal", "votes_hp1_maternal",
                 "hp1_parent", "agreement"]
        w = csv.DictWriter(f, fieldnames=bcols, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(block_rows, key=lambda r: (r["sample"], r["chrom"], r["start"])))

    # imprinting sanity check: expected-methylated side should be high, other low
    with open(f"{pre}.methylation_imprinting_check.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "gene", "expected_methylated", "paternal_meth", "maternal_meth", "consistent"])
        for child in sorted(levels):
            for pr in proms:
                if pr.gene not in IMPRINTED or pr.idx not in levels[child]:
                    continue
                pm = levels[child][pr.idx]["paternal"][0]
                mm = levels[child][pr.idx]["maternal"][0]
                exp = IMPRINTED[pr.gene]
                ok = ""
                if pm is not None and mm is not None:
                    hi, lo = (pm, mm) if exp == "paternal" else (mm, pm)
                    ok = "yes" if hi > lo else "no"
                w.writerow([child, pr.gene, exp, fmt(pm), fmt(mm), ok])
    logger.info(f"imprinting check -> {pre}.methylation_imprinting_check.tsv "
                f"(consistent=no across most genes means HP1/HP2 -> parent is inverted)")


if __name__ == "__main__":
    main()
