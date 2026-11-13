#!/usr/bin/env python3
"""
Mild-child variants on a parental copy the severe sibling did NOT inherit.

For each mild x severe sibling pair (family.json phenotypes), every variant
in the mild child's annotated master VCF is placed on a parental copy by
its INFO/ORIGIN (consolidate_variants.py: paternal / maternal / both), and
that copy is looked up in the hap-ibd sharing segments
(hapibd_segments.py, cohort.sib_parent_segments.tsv):

  paternal / maternal   reported if that copy is not_shared at the variant
  both (homozygous)     reported if either copy is not_shared
  unknown origin        reported only where BOTH copies are not_shared
                        (IBD0), since then it can't be on a shared copy

No genetic model and no filtering beyond that: every such variant is
written, with present_in_severe saying whether the severe sib carries the
same allele anyway (e.g. a parent homozygous for it, or recurrent), plus
snpEff (worst consequence / impact, genes) and population AF so you can
filter downstream. Variants on chromosomes the sharing step doesn't scan
(chrX/Y/M by default) are counted, not reported.

Outputs (-output_prefix cohort):
  cohort.mild_haplotype_variants.tsv.gz  one row per (mild, severe, variant)
  cohort.mild_haplotype_genes.tsv        per (mild, severe, gene): counts, absent-in-severe,
                                         HIGH/MODERATE, rare (min AF < -rare_af or no AF)
  cohort.mild_haplotype_funnel.tsv       per pair: variants in, placed, on a non-shared copy, reported
"""
import argparse
import bisect
import csv
import gzip
import json
import sys
from collections import defaultdict

CALLER_COLS = ("clair3_longphase", "sniffles")
IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}


def log(msg):
    sys.stderr.write(f"[mild_haplotype_variants] {msg}\n")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def info_get(info, key):
    for f in info.split(";"):
        if f.startswith(key + "="):
            return f[len(key) + 1:]
    return None


def zygosity(gt):
    a = [x for x in gt.replace("|", "/").split("/") if x not in ("", ".")]
    if not a or all(x == "0" for x in a):
        return None
    return "hom" if all(x != "0" for x in a) else "het"


def to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def read_vcf(path):
    """Yield one dict per called variant (GT from the first caller column that
    called it, in CALLER_COLS order, else the first sample column)."""
    cols = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            x = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                cols = {name: i for i, name in enumerate(x) if i >= 9}
                order = [cols[c] for c in CALLER_COLS if c in cols] or sorted(cols.values())
                continue
            if len(x) < 10:
                continue
            keys = x[8].split(":")
            gt = ps = ""
            for i in order:
                d = dict(zip(keys, x[i].split(":")))
                if zygosity(d.get("GT", "")):
                    gt, ps = d.get("GT", ""), d.get("PS", "")
                    break
            zyg = zygosity(gt)
            if not zyg:
                continue
            info = x[7]
            ann = info_get(info, "ANN")
            genes, worst, worst_eff = set(), "", ""
            if ann:
                for a in ann.split(","):
                    a = a.split("|")
                    if len(a) > 3:
                        if a[3]:
                            genes.add(a[3])
                        if IMPACT_RANK.get(a[2], -1) > IMPACT_RANK.get(worst, -1):
                            worst, worst_eff = a[2], a[1]
            yield {
                "chrom": x[0], "pos": int(x[1]), "id": x[2], "ref": x[3], "alt": x[4],
                "svtype": info_get(info, "SVTYPE") or "", "gt": gt, "ps": ps, "zyg": zyg,
                "origin": info_get(info, "ORIGIN") or "unknown",
                "origin_source": info_get(info, "ORIGIN_SOURCE") or "",
                "de_novo": "DE_NOVO" in info.split(";"),
                "callers": info_get(info, "CALLERS") or "",
                "genes": ",".join(sorted(genes)), "impact": worst, "effect": worst_eff,
                "gnomad": to_float(info_get(info, "gnomAD_AF_NFE")),
                "kg": to_float(info_get(info, "KG_AF_EUR")),
            }


class Sharing(object):
    def __init__(self, path):
        self.segs = defaultdict(list)
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                self.segs[(row["mild"], row["severe"], row["parent"], row["chrom"])].append(
                    (int(row["start"]), int(row["end"]), row["state"]))
        self.starts = {}
        for k, v in self.segs.items():
            v.sort()
            self.starts[k] = [s[0] for s in v]

    def state(self, mild, severe, parent, chrom, pos):
        k = (mild, severe, parent, chrom)
        if k not in self.segs:
            return "unscanned"
        i = bisect.bisect_right(self.starts[k], pos) - 1
        if i < 0:
            return "unscanned"
        lo, hi, st = self.segs[k][i]
        return st if lo <= pos <= hi else "unscanned"


class SevereIndex(object):
    """Exact small-variant keys (-> GT) and SV positions by (chrom, type)."""

    def __init__(self, path, sv_dist):
        self.small, self.sv, self.sv_dist = {}, defaultdict(list), sv_dist
        for v in read_vcf(path):
            if v["svtype"]:
                self.sv[(v["chrom"], v["svtype"])].append(v["pos"])
            else:
                self.small[(v["chrom"], v["pos"], v["ref"], v["alt"])] = v["gt"]
        for k in self.sv:
            self.sv[k].sort()

    def lookup(self, v):
        """(present, severe GT or '')"""
        if v["svtype"]:
            ps = self.sv.get((v["chrom"], v["svtype"]), [])
            i = bisect.bisect_left(ps, v["pos"])
            near = any(0 <= j < len(ps) and abs(ps[j] - v["pos"]) <= self.sv_dist for j in (i - 1, i))
            return near, ""
        gt = self.small.get((v["chrom"], v["pos"], v["ref"], v["alt"]))
        return gt is not None, gt or ""


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-family_json", required=True)
    p.add_argument("-vcf_manifest", required=True, help="childID<TAB>annotated master VCF")
    p.add_argument("-sib_parent_segments", required=True)
    p.add_argument("-sv_merge_dist", type=int, default=500)
    p.add_argument("-rare_af", type=float, default=0.01)
    p.add_argument("-output_prefix", default="cohort")
    args = p.parse_args()

    with open(args.family_json) as f:
        fam = json.load(f)
    manifest = {}
    with open(args.vcf_manifest) as f:
        for line in f:
            if line.strip():
                s, v = line.rstrip("\n").split("\t")[:2]
                manifest[s] = v
    pairs = []
    for m, mi in fam.items():
        if mi.get("phenotype") != "mild" or m not in manifest:
            continue
        for s, si in fam.items():
            if (si.get("phenotype") == "severe" and s in manifest
                    and (si.get("father"), si.get("mother")) == (mi.get("father"), mi.get("mother"))):
                pairs.append((m, s))
    log(f"pairs: {pairs}")
    sharing = Sharing(args.sib_parent_segments)

    cols = ["MILD", "SEVERE", "CHROM", "POS", "ID", "REF", "ALT", "SVTYPE", "GT", "zygosity",
            "parent_of_origin", "origin_source", "de_novo", "copy_status", "ibd",
            "present_in_severe", "severe_GT", "genes", "worst_effect", "worst_impact",
            "gnomAD_AF_NFE", "KG_AF_EUR", "callers"]
    pre = args.output_prefix
    out = gzip.open(f"{pre}.mild_haplotype_variants.tsv.gz", "wt", compresslevel=4)
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    funnel = []
    genes = defaultdict(lambda: {"n": 0, "absent": 0, "high": 0, "moderate": 0, "rare": 0, "min_af": None})
    severe_cache = {}

    for m, s in pairs:
        if s not in severe_cache:
            severe_cache = {s: SevereIndex(manifest[s], args.sv_merge_dist)}
        sev = severe_cache[s]
        n = {"variants": 0, "unscanned": 0, "origin_unresolved": 0, "on_shared_copy": 0, "reported": 0,
             "reported_absent_in_severe": 0}
        for v in read_vcf(manifest[m]):
            n["variants"] += 1
            st = {par: sharing.state(m, s, par, v["chrom"], v["pos"]) for par in ("paternal", "maternal")}
            if "unscanned" in st.values():
                n["unscanned"] += 1
                continue
            if v["origin"] in ("paternal", "maternal"):
                copies = [v["origin"]]
            elif v["origin"] == "both":
                copies = ["paternal", "maternal"]
            elif all(x == "not_shared" for x in st.values()):
                copies = ["paternal", "maternal"]       # IBD0: on a non-shared copy either way
            else:
                n["origin_unresolved"] += 1
                continue
            if not any(st[c] == "not_shared" for c in copies):
                n["on_shared_copy"] += 1
                continue
            present, sgt = sev.lookup(v)
            n["reported"] += 1
            n["reported_absent_in_severe"] += not present
            afs = [a for a in (v["gnomad"], v["kg"]) if a is not None]
            ibd = (st["paternal"] == "shared") + (st["maternal"] == "shared")
            w.writerow({
                "MILD": m, "SEVERE": s, "CHROM": v["chrom"], "POS": v["pos"], "ID": v["id"],
                "REF": v["ref"], "ALT": v["alt"], "SVTYPE": v["svtype"], "GT": v["gt"], "zygosity": v["zyg"],
                "parent_of_origin": v["origin"], "origin_source": v["origin_source"],
                "de_novo": "yes" if v["de_novo"] else "",
                "copy_status": ";".join(f"{c}:{st[c]}" for c in copies), "ibd": f"IBD{ibd}",
                "present_in_severe": "yes" if present else "no", "severe_GT": sgt,
                "genes": v["genes"], "worst_effect": v["effect"], "worst_impact": v["impact"],
                "gnomAD_AF_NFE": "" if v["gnomad"] is None else v["gnomad"],
                "KG_AF_EUR": "" if v["kg"] is None else v["kg"], "callers": v["callers"],
            })
            for g in (v["genes"].split(",") if v["genes"] else []):
                gs = genes[(m, s, g)]
                gs["n"] += 1
                gs["absent"] += not present
                gs["high"] += v["impact"] == "HIGH"
                gs["moderate"] += v["impact"] == "MODERATE"
                rare = (min(afs) < args.rare_af) if afs else True
                gs["rare"] += rare
                if afs and (gs["min_af"] is None or min(afs) < gs["min_af"]):
                    gs["min_af"] = min(afs)
        funnel.append(dict(MILD=m, SEVERE=s, **n))
        log(f"{m} vs {s}: " + ", ".join(f"{k}={v:,}" for k, v in n.items()))
    out.close()

    with open(f"{pre}.mild_haplotype_genes.tsv", "w", newline="") as f:
        gw = csv.writer(f, delimiter="\t")
        gw.writerow(["MILD", "SEVERE", "gene", "n_variants", "n_absent_in_severe", "n_HIGH", "n_MODERATE",
                     f"n_rare_AF_lt_{args.rare_af}", "min_AF"])
        for (m, s, g), d in sorted(genes.items(), key=lambda kv: (kv[0][0], kv[0][1], -kv[1]["high"],
                                                                   -kv[1]["moderate"], kv[0][2])):
            gw.writerow([m, s, g, d["n"], d["absent"], d["high"], d["moderate"], d["rare"],
                         "" if d["min_af"] is None else f"{d['min_af']:.6g}"])
    with open(f"{pre}.mild_haplotype_funnel.tsv", "w", newline="") as f:
        fcols = ["MILD", "SEVERE", "variants", "unscanned", "origin_unresolved", "on_shared_copy", "reported",
                 "reported_absent_in_severe"]
        fw = csv.DictWriter(f, fieldnames=fcols, delimiter="\t")
        fw.writeheader()
        fw.writerows(funnel)


if __name__ == "__main__":
    main()
