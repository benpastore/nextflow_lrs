#!/usr/bin/env python3
"""
Parent of origin for a child's structural variants (Sniffles), from two
independent lines of evidence:

1. Transmission -- is the SV also called in the father's / mother's own
   Sniffles VCF? Matched by chromosome, type (INS and DUP treated as
   compatible -- Sniffles often reports a tandem duplication as an insertion
   in one sample and a DUP in another), start within -merge_dist bp and, for
   sized SVs, |SVLEN| ratio >= -min_size_ratio (BND: same mate chromosome,
   mate position within -merge_dist).
     child het: father only -> paternal, mother only -> maternal,
                both -> both_parents (transmission can't tell which),
                neither -> not_in_parents (de novo candidate, or missed in a parent)
     child hom: both -> both_parents, one -> hom_one_parent, neither -> not_in_parents

2. Haplotype -- the child's Sniffles call carries its longphase haplotype and
   phase block (INFO/PHASE=hp,ps,...; trusted only when Sniffles' own QC
   says PASS). Inside each longphase block, the child's phased het SNVs at
   which the parents' Clair3 genotypes differ say which haplotype came from
   which parent (father 0/0 + mother carries ALT -> ALT is maternal, etc.).
   A block's HP1 is called paternal/maternal when >= -min_informative such
   SNVs vote and >= -min_agreement of them agree; the SV's haplotype then
   gives its origin. This resolves SVs transmission can't (present in both
   parents) and cross-checks the ones it can.

   The SNV VCF must be the phased VCF whose blocks match the BAM Sniffles ran
   on (in this pipeline: LONGPHASE round 1, *.longphase.vcf.gz). A child
   site absent from a parent's Clair3 VCF is taken as hom-ref (Clair3 only
   writes variant sites); -missing_parent_as_missing turns that off.

Outputs
  -output        one row per child SV: transmission, haplotype origin, final
                 origin + origin_source (transmission+haplotype / transmission /
                 haplotype / conflict)
  -block_output  per longphase block: informative SNVs and HP1's parent
"""
import argparse
import bisect
import gzip
import re
import sys
from collections import Counter, defaultdict

BND_MATE_RE = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")
COMPATIBLE = {"INS": "INS_DUP", "DUP": "INS_DUP"}


def log(msg):
    sys.stderr.write(f"[sv_trio_origin] {msg}\n")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def info_dict(info):
    d = {}
    for kv in info.split(";"):
        k, _, v = kv.partition("=")
        d[k] = v
    return d


def gt_class(gt):
    """'ref' / 'het' / 'hom' / None (missing)."""
    alleles = [a for a in gt.replace("|", "/").split("/") if a not in (".", "")]
    if not alleles:
        return None
    nonref = [a for a in alleles if a != "0"]
    if not nonref:
        return "ref"
    return "hom" if len(nonref) == len(alleles) and len(set(alleles)) == 1 else "het"


def sample_fields(fields):
    if len(fields) < 10:
        return {}
    return dict(zip(fields[8].split(":"), fields[9].split(":")))


# ------------------------------------------------------------------ SVs

class SV(object):
    __slots__ = ("chrom", "pos", "end", "svtype", "svlen", "vid", "filt", "gt", "hp", "ps", "mate")

    def __init__(self, fields):
        self.chrom, self.pos, self.vid, _ref, alt, _q, self.filt, info = fields[:8]
        self.pos = int(self.pos)
        d = info_dict(info)
        self.svtype = (d.get("SVTYPE") or alt.strip("<>")).upper()
        try:
            self.svlen = abs(int(float(d.get("SVLEN", "0"))))
        except ValueError:
            self.svlen = 0
        try:
            self.end = int(d.get("END", self.pos))
        except ValueError:
            self.end = self.pos
        self.gt = sample_fields(fields).get("GT", "")
        self.hp = self.ps = ""
        phase = d.get("PHASE", "").split(",")
        # PHASE=<hp>,<phaseset>,<nSupport>,<nTotal>,<qcHP>,<qcPhaseset>
        if len(phase) >= 6:
            if phase[4] == "PASS" and phase[0] in ("1", "2"):
                self.hp = phase[0]
            if phase[5] == "PASS" and phase[1] not in ("", "NULL"):
                self.ps = phase[1]
        self.mate = None
        if self.svtype == "BND":
            m = BND_MATE_RE.search(alt)
            if m:
                self.mate = (m.group(1), int(m.group(2)))

    @property
    def key_type(self):
        return COMPATIBLE.get(self.svtype, self.svtype)


def read_svs(path):
    out = []
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8:
                out.append(SV(fields))
    return out


class SVIndex(object):
    def __init__(self, svs):
        by = defaultdict(list)
        for s in svs:
            by[(s.chrom, s.key_type)].append(s)
        self.by = {}
        for k, lst in by.items():
            lst.sort(key=lambda s: s.pos)
            self.by[k] = ([s.pos for s in lst], lst)

    def match(self, sv, merge_dist, min_ratio):
        """Best-matching carried parental SV (GT with an ALT allele), or None."""
        entry = self.by.get((sv.chrom, sv.key_type))
        if not entry:
            return None
        starts, lst = entry
        lo = bisect.bisect_left(starts, sv.pos - merge_dist)
        hi = bisect.bisect_right(starts, sv.pos + merge_dist)
        best, best_d = None, None
        for p in lst[lo:hi]:
            if gt_class(p.gt) not in ("het", "hom"):
                continue
            if sv.svtype == "BND":
                if not (sv.mate and p.mate and sv.mate[0] == p.mate[0]
                        and abs(sv.mate[1] - p.mate[1]) <= merge_dist):
                    continue
            elif sv.svlen and p.svlen and min(sv.svlen, p.svlen) / float(max(sv.svlen, p.svlen)) < min_ratio:
                continue
            d = abs(p.pos - sv.pos)
            if best_d is None or d < best_d:
                best, best_d = p, d
        return best


# ------------------------------------------------------------------ SNV blocks

def read_child_phased_snvs(path):
    """{(chrom,pos,ref,alt): (hp1_is_alt, ps)} for phased het biallelic SNVs."""
    snvs = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            ref, alt = fields[3], fields[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            sf = sample_fields(fields)
            gt, ps = sf.get("GT", ""), sf.get("PS", "")
            if gt not in ("0|1", "1|0") or ps in ("", "."):
                continue
            snvs[(fields[0], fields[1], ref, alt)] = (gt[0] == "1", ps)
    return snvs


def read_parent_classes(path, wanted):
    out = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 10)
            if len(fields) < 10:
                continue
            key = (fields[0], fields[1], fields[3], fields[4])
            if key in wanted:
                out[key] = gt_class(sample_fields([x.rstrip("\n") for x in fields]).get("GT", ""))
    return out


def alt_origin(father, mother):
    """Which parent transmitted the child's ALT at a child-het site, or None."""
    f_alt, m_alt = father in ("het", "hom"), mother in ("het", "hom")
    if father == "ref" and m_alt:
        return "maternal"
    if mother == "ref" and f_alt:
        return "paternal"
    if father == "hom" and mother == "het":
        return "paternal"
    if mother == "hom" and father == "het":
        return "maternal"
    return None


def other(origin):
    return {"paternal": "maternal", "maternal": "paternal"}.get(origin)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-child", required=True)
    p.add_argument("-child_sv", required=True, help="Child Sniffles VCF (run with --phase)")
    p.add_argument("-father_sv", required=True)
    p.add_argument("-mother_sv", required=True)
    p.add_argument("-child_snv", required=True, help="Child phased SNV VCF whose blocks match Sniffles' BAM")
    p.add_argument("-father_snv", required=True, help="Father Clair3 VCF")
    p.add_argument("-mother_snv", required=True, help="Mother Clair3 VCF")
    p.add_argument("-merge_dist", type=int, default=500)
    p.add_argument("-min_size_ratio", type=float, default=0.5)
    p.add_argument("-min_informative", type=int, default=3)
    p.add_argument("-min_agreement", type=float, default=0.8)
    p.add_argument("-missing_parent_as_missing", action="store_true")
    p.add_argument("-output", required=True)
    p.add_argument("-block_output", required=True)
    args = p.parse_args()

    # ---- haplotype -> parent, per longphase block
    snvs = read_child_phased_snvs(args.child_snv)
    log(f"{len(snvs):,} phased het SNVs in the child")
    fcls = read_parent_classes(args.father_snv, snvs)
    mcls = read_parent_classes(args.mother_snv, snvs)
    missing = None if args.missing_parent_as_missing else "ref"
    votes = defaultdict(Counter)
    span = {}
    n_inf = 0
    for key, (hp1_alt, ps) in snvs.items():
        origin = alt_origin(fcls.get(key, missing), mcls.get(key, missing))
        if origin is None:
            continue
        n_inf += 1
        hp1 = origin if hp1_alt else other(origin)
        bk = (key[0], ps)
        votes[bk][hp1] += 1
        pos = int(key[1])
        lo, hi = span.get(bk, (pos, pos))
        span[bk] = (min(lo, pos), max(hi, pos))
    log(f"{n_inf:,} informative SNVs (parental genotypes differ) across {len(votes):,} block(s)")

    block_hp1 = {}
    with open(args.block_output, "w") as out:
        out.write("CHROM\tPS\tfirst_informative_pos\tlast_informative_pos\tn_informative\t"
                  "n_hp1_paternal\tn_hp1_maternal\tagreement\tHP1_origin\n")
        for (chrom, ps), c in sorted(votes.items(), key=lambda kv: (kv[0][0], int(kv[0][1]) if kv[0][1].isdigit() else 0)):
            n = c["paternal"] + c["maternal"]
            top, k = c.most_common(1)[0]
            agree = k / float(n)
            call = top if n >= args.min_informative and agree >= args.min_agreement else ""
            if call:
                block_hp1[(chrom, ps)] = call
            lo, hi = span[(chrom, ps)]
            out.write(f"{chrom}\t{ps}\t{lo}\t{hi}\t{n}\t{c['paternal']}\t{c['maternal']}\t{agree:.3f}\t{call}\n")
    log(f"{len(block_hp1):,} block(s) with HP1 assigned to a parent")

    # ---- SVs
    child = [s for s in read_svs(args.child_sv) if gt_class(s.gt) in ("het", "hom")]
    fidx = SVIndex(read_svs(args.father_sv))
    midx = SVIndex(read_svs(args.mother_sv))
    log(f"{len(child):,} child SV call(s) with a non-ref genotype")

    counts = Counter()
    ps_seen = ps_found = 0
    cols = ["CHILD", "CHROM", "POS", "END", "SVTYPE", "SVLEN", "ID", "FILTER", "CHILD_GT", "CHILD_HP", "CHILD_PS",
            "father_match", "father_gt", "mother_match", "mother_gt", "transmission",
            "haplotype_origin", "block_n_informative", "block_agreement", "origin", "origin_source"]
    with open(args.output, "w") as out:
        out.write("\t".join(cols) + "\n")
        for s in sorted(child, key=lambda s: (s.chrom, s.pos)):
            fm = fidx.match(s, args.merge_dist, args.min_size_ratio)
            mm = midx.match(s, args.merge_dist, args.min_size_ratio)
            zyg = gt_class(s.gt)
            if fm and mm:
                trans = "both_parents"
            elif zyg == "hom":
                trans = "hom_one_parent" if (fm or mm) else "not_in_parents"
            else:
                trans = "paternal" if fm else "maternal" if mm else "not_in_parents"

            hap, n_b, agr = "", "", ""
            if s.ps:
                ps_seen += 1
                bk = (s.chrom, s.ps)
                if bk in votes:
                    ps_found += 1
                    c = votes[bk]
                    n_b = c["paternal"] + c["maternal"]
                    agr = f"{c.most_common(1)[0][1] / float(n_b):.3f}"
                if zyg == "het" and s.hp and bk in block_hp1:
                    hap = block_hp1[bk] if s.hp == "1" else other(block_hp1[bk])

            t = trans if trans in ("paternal", "maternal") else ""
            if t and hap:
                origin, src = (t, "transmission+haplotype") if t == hap else ("", "conflict")
            elif t:
                origin, src = t, "transmission"
            elif hap:
                origin, src = hap, "haplotype"
            else:
                origin, src = "", ""
            counts[(trans, src)] += 1
            out.write("\t".join(str(x) for x in [
                args.child, s.chrom, s.pos, s.end, s.svtype, s.svlen, s.vid, s.filt, s.gt, s.hp, s.ps,
                fm.vid if fm else "", fm.gt if fm else "", mm.vid if mm else "", mm.gt if mm else "",
                trans, hap, n_b, agr, origin, src]) + "\n")

    if ps_seen:
        log(f"{ps_found:,}/{ps_seen:,} phased child SVs fall in a block seen in the SNV VCF"
            + ("" if ps_found else " -- block IDs don't match: is -child_snv the VCF of the BAM Sniffles ran on?"))
    by_trans = Counter()
    by_src = Counter()
    for (t, src), n in counts.items():
        by_trans[t] += n
        by_src[src or "unresolved"] += n
    log("transmission: " + ", ".join(f"{n:,} {t}" for t, n in by_trans.most_common()))
    log("origin source: " + ", ".join(f"{n:,} {s}" for s, n in by_src.most_common()))


if __name__ == "__main__":
    main()
