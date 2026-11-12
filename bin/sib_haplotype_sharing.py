#!/usr/bin/env python3
"""
Sibling haplotype sharing (identity-by-descent) from family genotypes
alone: father + mother + their children.

Why: each child carries one paternal and one maternal haplotype. At any
locus, two siblings got either the same or the other haplotype from each
parent, so they share 0, 1 or 2 haplotypes (IBD0 / IBD1 / IBD2). A variant
on a haplotype the severe sibling also inherited can't, on its own, explain
why the mild sibling is mild. This script works out, along the genome,
which parental haplotypes each mild/severe sibling pair shares, so the
modifier search (bin/haplotype_modifier.py) can tell shared from
mild-specific haplotypes.

How (per parent, per mild x severe sibling pair):
  1. Informative sites: biallelic SNVs where that parent is het and the
     other parent is homozygous (ref or alt). There, each child's allele
     from the het parent can be read straight off its genotype:
        other parent 0/0: child 0/0 -> got REF, child 0/1 -> got ALT
        other parent 1/1: child 1/1 -> got ALT, child 0/1 -> got REF
  2. Same allele in both sibs = "match". Shared haplotype -> nearly all
     sites match; different haplotypes -> nearly all mismatch (the parent
     is het there, so its two haplotypes differ by definition).
  3. A two-state HMM (shared / not_shared) with genotype-error rate
     -error_rate and a distance-scaled switch probability (-switch_rate per
     bp, ~ the crossover rate) smooths the match/mismatch run into segments
     (Viterbi). Segments with < -min_segment_sites informative sites are
     merged into their neighbours, so short runs of genotyping error
     (e.g. in segmental duplications) don't show up as fake crossovers.
  Paternal and maternal segments are then intersected into IBD0/1/2.

No phasing is needed -- only genotypes. (Which parental haplotype each
variant is on is assigned separately, in bin/consolidate_variants.py, as
INFO/ORIGIN.)

Conventions (same as bin/sv_trio_origin.py):
  - a site missing from a sample's VCF is taken as hom-ref (Clair3 only
    writes variant sites); a record that is present but non-PASS, or whose
    GT is missing, is "unknown" and the site is skipped
  - chrX/chrY/chrM are skipped by default (-exclude_chroms): the father is
    hemizygous on X, so the two-state model above doesn't hold there

Inputs
  -family_json   {child: {father, mother, phenotype: mild|severe}}
  -vcf_manifest  sampleID<TAB>vcf for every parent (Clair3) and child
                 (*.trio.whatshapphase.vcf.gz)

Outputs (prefix -output_prefix)
  .sib_ibd_segments.tsv       mild x severe pair: IBD0/1/2 segments
  .sib_parent_segments.tsv    same, per parent, with informative-site counts (QC)
  .sib_ibd_summary.tsv        per pair: fraction of the scanned genome in IBD0/1/2
                              (expect roughly 25/50/25% for full siblings)
"""
import argparse
import gzip
import json
import math
import sys
from collections import Counter, defaultdict

REF, HET, HOM = 0, 1, 2
UNKNOWN = None


def log(msg):
    sys.stderr.write(f"[sib_haplotype_sharing] {msg}\n")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def gt_code(gt):
    """0/1/2 copies of ALT for a biallelic GT, or None if missing/odd."""
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or any(a not in ("0", "1") for a in alleles):
        return UNKNOWN
    return alleles.count("1")


def load_snvs(path, exclude_chroms):
    """{chrom: {pos: (alt, code)}} for biallelic SNV records.

    code is None for a record that exists but can't be trusted (non-PASS,
    missing GT) -- that is "unknown", not hom-ref."""
    out = defaultdict(dict)
    n = 0
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t", 10)
            if len(fields) < 10:
                continue
            chrom, pos, _id, ref, alt, _q, filt = fields[:7]
            if chrom in exclude_chroms or len(ref) != 1 or len(alt) != 1 or alt in (".", "*"):
                continue
            fmt = fields[8].split(":")
            vals = fields[9].split(":")
            d = dict(zip(fmt, vals))
            gt = d.get("GT", "")
            out[chrom][int(pos)] = (alt, gt_code(gt) if filt in ("PASS", ".") else UNKNOWN)
            n += 1
    log(f"  {path}: {n:,} biallelic SNV records")
    return out


def code_at(sample, chrom, pos, alt):
    """Genotype code of `sample` at a site whose ALT is `alt`. Absent -> REF
    (Clair3 writes variant sites only); different ALT at the same position,
    or an untrusted record -> None."""
    rec = sample.get(chrom, {}).get(pos)
    if rec is None:
        return REF
    if rec[0] != alt:
        return UNKNOWN
    return rec[1]


def transmitted_allele(child_code, other_parent_code):
    """Allele (0 REF / 1 ALT) a child got from the het parent, given the
    other parent is homozygous; None if Mendelian-inconsistent/unknown."""
    if child_code is UNKNOWN:
        return None
    if other_parent_code == REF:
        return {REF: 0, HET: 1}.get(child_code)
    if other_parent_code == HOM:
        return {HOM: 1, HET: 0}.get(child_code)
    return None


# ------------------------------------------------------------------ HMM

def viterbi(obs, err, rate):
    """obs: [(pos, matched)] sorted by pos. Returns a state per obs,
    1 = shared, 0 = not shared."""
    if not obs:
        return []
    lm = (math.log(err), math.log(1 - err))      # log P(match | not_shared, shared)
    lx = (math.log(1 - err), math.log(err))      # log P(mismatch | ...)
    e = lm if obs[0][1] else lx
    v0, v1 = math.log(0.5) + e[0], math.log(0.5) + e[1]
    back = bytearray(len(obs) - 1)
    prev = obs[0][0]
    for i in range(1, len(obs)):
        pos, m = obs[i]
        p = min(0.5, max(1e-12, (pos - prev) * rate))
        ls, lk = math.log(p), math.log(1 - p)
        e = lm if m else lx
        a0, b0 = v0 + lk, v1 + ls    # into not_shared: stay / switch from shared
        a1, b1 = v1 + lk, v0 + ls    # into shared:     stay / switch from not_shared
        bp0 = 0 if a0 >= b0 else 1
        bp1 = 1 if a1 >= b1 else 0
        v0, v1 = max(a0, b0) + e[0], max(a1, b1) + e[1]
        back[i - 1] = bp0 | (bp1 << 1)
        prev = pos
    s = 1 if v1 > v0 else 0
    states = [s]
    for bp in reversed(back):
        s = (bp >> s) & 1
        states.append(s)
    states.reverse()
    return states


def to_segments(obs, states):
    """[[start, end, state, n_sites, n_discordant]] runs of equal state."""
    segs = []
    for (pos, m), s in zip(obs, states):
        disc = 0 if m == bool(s) else 1
        if segs and segs[-1][2] == s:
            seg = segs[-1]
            seg[1] = pos
            seg[3] += 1
            seg[4] += disc
        else:
            segs.append([pos, pos, s, 1, disc])
    return segs


def merge_short(segs, min_sites):
    """Flip the shortest segment below min_sites into its neighbours'
    state, repeatedly, so isolated short runs don't become crossovers.
    A flipped segment's discordant count becomes its concordant count."""
    segs = [list(s) for s in segs]
    while len(segs) > 1:
        idx = min(range(len(segs)), key=lambda i: segs[i][3])
        if segs[idx][3] >= min_sites:
            break
        seg = segs[idx]
        seg[2] = 1 - seg[2]
        seg[4] = seg[3] - seg[4]
        merged = []
        for s in segs:
            if merged and merged[-1][2] == s[2]:
                m = merged[-1]
                m[1] = s[1]
                m[3] += s[3]
                m[4] += s[4]
            else:
                merged.append(s)
        segs = merged
    return segs


def intersect(pat, mat):
    """Two-pointer overlap of paternal and maternal segment lists ->
    [(start, end, pat_state, mat_state)]."""
    out = []
    i = j = 0
    while i < len(pat) and j < len(mat):
        lo = max(pat[i][0], mat[j][0])
        hi = min(pat[i][1], mat[j][1])
        if lo <= hi:
            out.append((lo, hi, pat[i][2], mat[j][2]))
        if pat[i][1] < mat[j][1]:
            i += 1
        else:
            j += 1
    return out


# ------------------------------------------------------------------ main

def chrom_key(c):
    s = c[3:] if c.startswith("chr") else c
    return (0, int(s), "") if s.isdigit() else (1, 0, s)


def pair_observations(het_parent, other_parent, sib_a, sib_b, chrom):
    obs = []
    for pos, (alt, code) in sorted(het_parent.get(chrom, {}).items()):
        if code != HET:
            continue
        other = code_at(other_parent, chrom, pos, alt)
        if other not in (REF, HOM):
            continue
        a = transmitted_allele(code_at(sib_a, chrom, pos, alt), other)
        b = transmitted_allele(code_at(sib_b, chrom, pos, alt), other)
        if a is None or b is None:
            continue
        obs.append((pos, a == b))
    return obs


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-family_json", required=True)
    p.add_argument("-vcf_manifest", required=True, help="sampleID<TAB>vcf, parents and children")
    p.add_argument("-error_rate", type=float, default=0.05,
                   help="Per-site genotype-error rate for the HMM (default 0.05)")
    p.add_argument("-switch_rate", type=float, default=1e-8,
                   help="Per-bp probability of a sharing switch (crossover in either sib), default 1e-8")
    p.add_argument("-min_segment_sites", type=int, default=50,
                   help="Segments with fewer informative sites get merged into their neighbours (default 50)")
    p.add_argument("-exclude_chroms", default="chrX,chrY,chrM,X,Y,MT")
    p.add_argument("-output_prefix", default="cohort")
    args = p.parse_args()

    exclude = {c for c in args.exclude_chroms.split(",") if c}
    with open(args.family_json) as f:
        families = json.load(f)
    manifest = {}
    with open(args.vcf_manifest) as f:
        for line in f:
            if line.strip():
                s, v = line.rstrip("\n").split("\t")[:2]
                manifest[s] = v

    by_parents = defaultdict(lambda: {"mild": [], "severe": [], "other": []})
    for child, info in families.items():
        fa, mo = info.get("father"), info.get("mother")
        if not fa or not mo:
            continue
        ph = info.get("phenotype")
        by_parents[(fa, mo)][ph if ph in ("mild", "severe") else "other"].append(child)

    pre = args.output_prefix
    seg_out = open(f"{pre}.sib_parent_segments.tsv", "w")
    ibd_out = open(f"{pre}.sib_ibd_segments.tsv", "w")
    sum_out = open(f"{pre}.sib_ibd_summary.tsv", "w")
    seg_out.write("family\tmild\tsevere\tparent\tchrom\tstart\tend\tstate\tn_sites\tn_discordant\n")
    ibd_out.write("family\tmild\tsevere\tchrom\tstart\tend\tpaternal\tmaternal\tibd\n")
    sum_out.write("family\tmild\tsevere\tbp_scanned\tfrac_IBD0\tfrac_IBD1\tfrac_IBD2\t"
                  "paternal_sites\tpaternal_discordance\tmaternal_sites\tmaternal_discordance\n")

    state_name = {1: "shared", 0: "not_shared"}
    for (fa, mo), grp in sorted(by_parents.items()):
        fam = f"{fa}+{mo}"
        if fa not in manifest or mo not in manifest:
            log(f"family {fam}: parent VCF(s) missing from manifest, skipping")
            continue
        children = [c for c in grp["mild"] + grp["severe"] + grp["other"] if c in manifest]
        if not children:
            continue
        log(f"family {fam}: children {', '.join(children)}")
        father = load_snvs(manifest[fa], exclude)
        mother = load_snvs(manifest[mo], exclude)
        kids = {c: load_snvs(manifest[c], exclude) for c in children}

        chroms = sorted(set(father) | set(mother), key=chrom_key)
        for mild in [c for c in grp["mild"] if c in kids]:
            for sev in [c for c in grp["severe"] if c in kids]:
                bp = Counter()
                site_tot = Counter()
                disc_tot = Counter()
                for chrom in chroms:
                    per_parent = {}
                    for name, het_p, other_p in (("paternal", father, mother), ("maternal", mother, father)):
                        obs = pair_observations(het_p, other_p, kids[mild], kids[sev], chrom)
                        segs = merge_short(to_segments(obs, viterbi(obs, args.error_rate, args.switch_rate)),
                                           args.min_segment_sites)
                        per_parent[name] = segs
                        for s in segs:
                            site_tot[name] += s[3]
                            disc_tot[name] += s[4]
                            seg_out.write(f"{fam}\t{mild}\t{sev}\t{name}\t{chrom}\t{s[0]}\t{s[1]}\t"
                                          f"{state_name[s[2]]}\t{s[3]}\t{s[4]}\n")
                    for lo, hi, ps_, ms_ in intersect(per_parent["paternal"], per_parent["maternal"]):
                        ibd = ps_ + ms_
                        bp[ibd] += hi - lo + 1
                        ibd_out.write(f"{fam}\t{mild}\t{sev}\t{chrom}\t{lo}\t{hi}\t"
                                      f"{state_name[ps_]}\t{state_name[ms_]}\tIBD{ibd}\n")
                tot = sum(bp.values()) or 1
                fr = [bp[k] / float(tot) for k in (0, 1, 2)]
                pd_ = disc_tot["paternal"] / float(site_tot["paternal"] or 1)
                md_ = disc_tot["maternal"] / float(site_tot["maternal"] or 1)
                sum_out.write(f"{fam}\t{mild}\t{sev}\t{sum(bp.values())}\t{fr[0]:.3f}\t{fr[1]:.3f}\t{fr[2]:.3f}\t"
                              f"{site_tot['paternal']}\t{pd_:.4f}\t{site_tot['maternal']}\t{md_:.4f}\n")
                log(f"  {mild} vs {sev}: IBD0 {fr[0]:.1%}, IBD1 {fr[1]:.1%}, IBD2 {fr[2]:.1%} of "
                    f"{sum(bp.values()):,} bp; within-segment discordance "
                    f"pat {pd_:.2%}, mat {md_:.2%}")
                if not (0.1 <= fr[2] <= 0.45 and 0.1 <= fr[0] <= 0.45):
                    log(f"  WARNING {mild} vs {sev}: IBD fractions far from the ~25/50/25% expected for "
                        f"full siblings -- check sample identities / family.json")
        del father, mother, kids

    for fh in (seg_out, ibd_out, sum_out):
        fh.close()


if __name__ == "__main__":
    main()
