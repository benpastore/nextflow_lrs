#!/usr/bin/env python3
"""
Prepare a WhatsHap family-phased VCF for hap-ibd: orient every child's phase
to paternal|maternal along the whole chromosome, keep the children only,
and keep only sites hap-ibd can use.

Why: hap-ibd reads each sample's GT order (a|b) as two chromosome-long
haplotypes and ignores phase blocks (PS). WhatsHap's orientation can flip
from one block to the next, which hap-ibd would see as a haplotype switch
at every block boundary. Each child block's first allele (HP1) is assigned
to a parent by voting over its Mendelian-informative het SNVs (one parent
carries the ALT, the other doesn't), using the parents' genotypes already
in the same family VCF; blocks whose HP1 is maternal are flipped. After
this, hap-ibd haplotype 1 = paternal and haplotype 2 = maternal for every
child, so its segments read directly as "same paternal copy" / "same
maternal copy".

Sites kept (all children together, so every output site is complete):
  biallelic SNV, FILTER PASS/., on an included chromosome, and for every
  child either homozygous (0/0 -> 0|0, 1/1 -> 1|1; missing = not called by
  Clair3 = 0|0) or a phased het in a block assigned to a parent, with
  GQ >= -min_gq when the VCF has GQ (ONT genotype errors break hap-ibd's
  identity-by-state seeds, so low-confidence hets are dropped).

    orient_family_phase.py -vcf FAM.family.whatshapphase.vcf.gz \\
        -family_json family.json -output FAM.oriented.vcf.gz \\
        -blocks FAM.orient_blocks.tsv
"""
import argparse
import gzip
import json
import sys
from collections import defaultdict


def log(msg):
    sys.stderr.write(f"[orient_family_phase] {msg}\n")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def gt_class(gt):
    """'ref' / 'het' / 'hom' for a biallelic GT; missing counts as ref
    (Clair3 only writes variant sites, bcftools merge fills them with ./.)."""
    a = [x for x in gt.replace("|", "/").split("/") if x not in ("", ".")]
    if not a or all(x == "0" for x in a):
        return "ref"
    return "hom" if all(x != "0" for x in a) else "het"


def alt_parent(f, m):
    """Which parent transmitted a child-het ALT, or None if uninformative."""
    if f == "ref" and m in ("het", "hom"):
        return "maternal"
    if m == "ref" and f in ("het", "hom"):
        return "paternal"
    if f == "hom" and m == "het":
        return "paternal"
    if m == "hom" and f == "het":
        return "maternal"
    return None


OTHER = {"paternal": "maternal", "maternal": "paternal"}


def sample_fields(fmt, val):
    return dict(zip(fmt.split(":"), val.split(":")))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-vcf", required=True, help="WHATSHAP_PHASE_FAMILY <family>.family.whatshapphase.vcf.gz")
    p.add_argument("-family_json", required=True)
    p.add_argument("-output", required=True, help="oriented, children-only VCF (.vcf.gz)")
    p.add_argument("-blocks", required=True, help="per child phase block -> HP1 parent TSV")
    p.add_argument("-min_votes", type=int, default=3)
    p.add_argument("-min_agree", type=float, default=0.8)
    p.add_argument("-min_gq", type=int, default=20)
    p.add_argument("-exclude_chroms", default="chrX,chrY,chrM")
    args = p.parse_args()
    exclude = {c for c in args.exclude_chroms.split(",") if c}

    with open(args.family_json) as f:
        fam = json.load(f)

    # pass 1: sample columns, which children belong here, block votes
    header, samples = [], []
    with opener(args.vcf) as f:
        for line in f:
            if line.startswith("##"):
                header.append(line)
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
    col = {s: i + 9 for i, s in enumerate(samples)}
    children = [c for c in samples if c in fam and fam[c].get("father") in col and fam[c].get("mother") in col]
    if not children:
        sys.exit(f"[orient_family_phase] no child in {args.vcf} has both parents in the VCF ({samples})")
    log(f"children: {children}; samples in VCF: {samples}")

    votes = {c: defaultdict(lambda: {"paternal": 0, "maternal": 0}) for c in children}
    spans = {c: {} for c in children}
    with opener(args.vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if x[0] in exclude or len(x[3]) != 1 or len(x[4]) != 1 or "," in x[4]:
                continue
            fmt = x[8]
            for c in children:
                d = sample_fields(fmt, x[col[c]])
                gt, ps = d.get("GT", ""), d.get("PS", "")
                if gt not in ("0|1", "1|0") or ps in ("", "."):
                    continue
                bk = (x[0], ps)
                pos = int(x[1])
                lo_hi = spans[c].get(bk)
                spans[c][bk] = (pos, pos) if lo_hi is None else (min(lo_hi[0], pos), max(lo_hi[1], pos))
                fa = gt_class(sample_fields(fmt, x[col[fam[c]["father"]]]).get("GT", ""))
                mo = gt_class(sample_fields(fmt, x[col[fam[c]["mother"]]]).get("GT", ""))
                o = alt_parent(fa, mo)
                if o:
                    votes[c][bk][o if gt[0] == "1" else OTHER[o]] += 1

    hp1 = {c: {} for c in children}
    with open(args.blocks, "w") as out:
        out.write("child\tchrom\tPS\tstart\tend\tvotes_hp1_paternal\tvotes_hp1_maternal\thp1_parent\n")
        for c in children:
            n_ok = 0
            for bk, (s, e) in sorted(spans[c].items(), key=lambda kv: (kv[0][0], kv[1][0])):
                v = votes[c].get(bk, {"paternal": 0, "maternal": 0})
                n = v["paternal"] + v["maternal"]
                best = "paternal" if v["paternal"] >= v["maternal"] else "maternal"
                ok = n >= args.min_votes and v[best] / float(n) >= args.min_agree
                if ok:
                    hp1[c][bk] = best
                    n_ok += 1
                out.write(f"{c}\t{bk[0]}\t{bk[1]}\t{s}\t{e}\t{v['paternal']}\t{v['maternal']}\t{best if ok else ''}\n")
            log(f"{c}: {n_ok:,}/{len(spans[c]):,} phase blocks assigned a parent")

    # pass 2: write oriented children-only VCF
    kept = dropped = 0
    with opener(args.vcf) as f, gzip.open(args.output, "wt", compresslevel=4) as out:
        for h in header:
            if not h.startswith("##FORMAT") and not h.startswith("##INFO"):
                out.write(h)
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, oriented paternal|maternal">\n')
        out.write("##orient_family_phase=haplotype 1 = paternal, haplotype 2 = maternal\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(children) + "\n")
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if (x[0] in exclude or len(x[3]) != 1 or len(x[4]) != 1 or "," in x[4]
                    or x[6] not in ("PASS", ".")):
                continue
            fmt = x[8]
            gts = []
            for c in children:
                d = sample_fields(fmt, x[col[c]])
                gt, ps, gq = d.get("GT", ""), d.get("PS", ""), d.get("GQ", "")
                cls = gt_class(gt)
                if cls == "ref":
                    gts.append("0|0")
                elif cls == "hom":
                    gts.append("1|1")
                else:
                    o = hp1[c].get((x[0], ps)) if gt in ("0|1", "1|0") else None
                    try:
                        low = gq not in ("", ".") and float(gq) < args.min_gq
                    except ValueError:
                        low = False
                    if o is None or low:
                        gts = None
                        break
                    gts.append(gt if o == "paternal" else gt[2] + "|" + gt[0])
            if gts is None or all(g == "0|0" for g in gts):
                dropped += 1
                continue
            kept += 1
            out.write("\t".join(x[:5]) + "\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
    log(f"{kept:,} sites written, {dropped:,} dropped (unoriented / low-GQ het, or no child carries ALT)")


if __name__ == "__main__":
    main()
