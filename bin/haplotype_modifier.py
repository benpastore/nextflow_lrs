#!/usr/bin/env python3
"""
Haplotype-aware modifier search: compare mild vs severe siblings by
haplotype FIRST, then run the recessive / dominant models.

Inputs already carry each variant's parental haplotype: INFO/ORIGIN in the
consolidated+snpEff VCF (paternal / maternal / both / unknown, assigned by
bin/consolidate_variants.py from parental genotypes, the child's WhatsHap
trio phase blocks, and bin/sv_trio_origin.py for SVs). This script adds
whether a severe sibling shares that haplotype (bin/sib_haplotype_sharing.py).

Step 1 -- mild vs severe haplotype comparison (cohort.haplotype_variants.tsv)
  For every in-scope variant of every mild child (annotation scope as before,
  bin/_modifier_common.py's parse_consolidated_vcf):
    haplotype_status  is the parental haplotype carrying it also carried by a
                        severe sib?  shared / mild_specific / unknown
    category          shared_haplotype  -- on a haplotype a severe sib also has
                      mild_specific     -- on a haplotype no severe sib has, and
                                           no severe sample carries the allele
                      shared_allele     -- haplotype not shared, but a severe
                                           sample carries the same allele anyway
                                           (e.g. a parent homozygous for it)
                      de_novo           -- in neither parent nor any severe sample
                      unknown_haplotype -- origin/sharing couldn't be resolved

Step 2 -- genetic models on those haplotypes
  dominant  (cohort.dominant_modifier.tsv)
      a variant whose category is mild_specific or de_novo. Variants on a
      haplotype shared with a severe sib are removed here -- this is the
      "remove haplotypes common to mild and severe" step. Then, gene level
      (as in bin/dominant_modifier_consolidated.py): drop the gene if any
      severe sib or concordant-severe sample carries ANY in-scope variant
      in it -- one hit is enough under this model, so severe already
      having one rules the gene out.
  recessive (cohort.recessive_modifier.tsv)
      the gene is biallelic in mild -- homozygous, or qualifying variants on
      BOTH the paternal and the maternal haplotype (or two hets in trans by
      read phasing, when origin is unknown) -- and
        - no severe sib (or concordant-severe sample) is itself biallelic
          in that gene, and
        - the mild child is not IBD2 with a severe sib across the gene
          (identical diplotype; guards against variants missed in severe).
      Shared-haplotype variants are NOT removed before this test: in a
      compound het, one hit may sit on the haplotype the severe sib shares
      and the other on the mild-specific haplotype --
          mild   = hapA(m1) + hapB(m2)    -> biallelic
          severe = hapA(m1) + hapC        -> monoallelic
      m1 is "shared", but the diplotype still differs, so both m1 and m2 are
      reported (reason says which side is shared vs mild_specific). Dropping
      shared-haplotype variants first would lose exactly these cases.

Output TSVs keep the schema the downstream steps (AlphaGenome,
AlphaMissense, variant report) already read, with extra columns added.
"""
import argparse
import bisect
import csv
from collections import Counter, defaultdict

from _modifier_common import (
    GnomadNFEAnnotator,
    annotate_final_pass,
    build_severe_index,
    in_severe,
    load_family_groups,
    load_vcf_manifest,
    logger,
    parse_consolidated_vcf,
)

LOCAL_PHASE_CALLERS = {"clair3_longphase", "clair3", "sniffles"}

BASE_COLS = ["#CHROM", "POS", "ID", "REF", "ALT", "gene", "zygosity", "reason",
             "evidence_caller", "callers_detected", "GT", "PS", "samples", "N_samples",
             "gnomad_AF_NFE", "kg_AF_EUR"]
HAP_COLS = ["parent_of_origin", "origin_source", "haplotype_status", "category", "ibd_vs_severe"]


# ------------------------------------------------------------------ sibling sharing

class Sharing:
    """Per (mild, severe, parent): sorted segments, for position lookup."""

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

    def shared(self, mild, severe, parent, chrom, pos):
        """True / False, or None outside any segment (unscanned chrom,
        breakpoint gap, or chromosome end)."""
        k = (mild, severe, parent, chrom)
        if k not in self.segs:
            return None
        i = bisect.bisect_right(self.starts[k], pos) - 1
        if i < 0:
            return None
        lo, hi, st = self.segs[k][i]
        return st if lo <= pos <= hi else None


def side_status(sharing, mild, severes, side, chrom, pos):
    """('shared', [sibs]) / ('mild_specific', []) / ('unknown', [])."""
    vals = [sharing.shared(mild, s, side, chrom, pos) for s in severes]
    with_ = [s for s, v in zip(severes, vals) if v is True]
    if with_:
        return "shared", with_
    if vals and all(v is False for v in vals):
        return "mild_specific", []
    return "unknown", []


def ibd_label(sharing, mild, severes, chrom, pos):
    parts = []
    for s in severes:
        p = sharing.shared(mild, s, "paternal", chrom, pos)
        m = sharing.shared(mild, s, "maternal", chrom, pos)
        n = None if p is None or m is None else int(p) + int(m)
        fmt = lambda x: "?" if x is None else ("shared" if x else "not_shared")
        parts.append(f"{s}:{'IBD?' if n is None else f'IBD{n}'}(pat={fmt(p)},mat={fmt(m)})")
    return ";".join(parts)


def is_ibd2(sharing, mild, severe, chrom, pos):
    return (sharing.shared(mild, severe, "paternal", chrom, pos) is True
            and sharing.shared(mild, severe, "maternal", chrom, pos) is True)


# ------------------------------------------------------------------ step 1

def compare_haplotypes(mild, severes, records, sharing, small_keys, sv_positions, sv_merge_dist):
    """{vid: hap-info dict} for every in-scope variant of one mild child."""
    out = {}
    for vid, rec in records.items():
        origin = rec["origin"]
        chrom, pos = rec["chrom"], rec["pos"]
        st_p, w_p = side_status(sharing, mild, severes, "paternal", chrom, pos)
        st_m, w_m = side_status(sharing, mild, severes, "maternal", chrom, pos)
        if origin == "paternal":
            status, with_ = st_p, w_p
        elif origin == "maternal":
            status, with_ = st_m, w_m
        elif origin == "both":
            # homozygous: one copy on each parental haplotype
            status = "mild_specific" if "mild_specific" in (st_p, st_m) else (
                "shared" if st_p == st_m == "shared" else "unknown")
            with_ = sorted(set(w_p) | set(w_m)) if status == "shared" else []
        else:
            # unknown origin: still resolvable when both sides agree (IBD0 / IBD2)
            status, with_ = (st_p, sorted(set(w_p) | set(w_m))) if st_p == st_m else ("unknown", [])

        present = in_severe(rec, small_keys, sv_positions, sv_merge_dist)
        if rec["de_novo"] and not present:
            category = "de_novo"
        elif status == "shared":
            category = "shared_haplotype"
        elif status == "mild_specific":
            category = "shared_allele" if present else "mild_specific"
        else:
            category = "unknown_haplotype"

        out[vid] = {
            "origin": origin, "source": rec["origin_source"],
            "status": status if not with_ else f"shared({','.join(with_)})",
            "category": category, "present_in_severe": present,
            "ibd": ibd_label(sharing, mild, severes, chrom, pos),
        }
    return out


# ------------------------------------------------------------------ step 2

def trans_by_readphase(ra, rb):
    return (ra["caller"] in LOCAL_PHASE_CALLERS and rb["caller"] in LOCAL_PHASE_CALLERS
            and ra["gt"] in ("0|1", "1|0") and rb["gt"] in ("0|1", "1|0")
            and ra["ps"] not in ("", ".") and ra["ps"] == rb["ps"] and ra["gt"] != rb["gt"])


def biallelic_genes(records):
    """{gene: (kind, [vids])} where this child has a qualifying hit on both
    haplotypes: homozygous, compound_het_parental (paternal + maternal by
    ORIGIN), or compound_het_readphase (trans by read phasing, when at
    least one partner's origin is unknown)."""
    by_gene = defaultdict(list)
    for vid, rec in records.items():
        for g in rec["genes"]:
            by_gene[g].append((vid, rec))
    out = {}
    for g, items in by_gene.items():
        homs = [v for v, r in items if r["origin"] == "both" or r["zygosity"] == "homo"]
        if homs:
            out[g] = ("homozygous", homs)
            continue
        pat = [v for v, r in items if r["origin"] == "paternal"]
        mat = [v for v, r in items if r["origin"] == "maternal"]
        if pat and mat:
            out[g] = ("compound_het_parental", pat + mat)
            continue
        pairs = set()
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                (va, ra), (vb, rb) = items[i], items[j]
                if "unknown" not in (ra["origin"], rb["origin"]):
                    continue        # both known and same side -> cis
                if trans_by_readphase(ra, rb):
                    pairs.update((va, vb))
        if pairs:
            out[g] = ("compound_het_readphase", sorted(pairs))
    return out


def side_summary(vids, hap):
    """'mild_specific' / 'shared(112_severe)' / ... for one haplotype's hits."""
    sts = sorted({hap[v]["status"] for v in vids})
    return "|".join(sts) if sts else "none"


# ------------------------------------------------------------------ run

def run(args):
    discordant, concordant = load_family_groups(args.family_json)
    manifest = load_vcf_manifest(args.vcf_manifest)

    children = {m for m, _ in discordant} | {s for _, ss in discordant for s in ss} | set(concordant)
    records = {}
    for c in sorted(children):
        if c not in manifest:
            logger.warning(f"{c}: no VCF in manifest, skipping")
            records[c] = {}
        else:
            records[c] = parse_consolidated_vcf(manifest[c])
            n_known = sum(r["origin"] != "unknown" for r in records[c].values())
            logger.info(f"{c}: parental haplotype known for {n_known:,}/{len(records[c]):,} in-scope variants")

    sharing = Sharing(args.sib_parent_segments)
    bi = {c: biallelic_genes(records[c]) for c in children}
    concordant_bi = set().union(*(bi[c] for c in concordant)) if concordant else set()
    concordant_records = [records[c] for c in concordant]

    hap_rows = []
    funnel = []
    dom_pass = {}   # (chrom,pos,ref,alt,gene) -> info
    rec_pass = {}
    dom_gene_samples = defaultdict(set)
    rec_gene_samples = defaultdict(lambda: defaultdict(set))

    for mild, severes in discordant:
        recs = records.get(mild, {})
        small_keys, sv_positions = build_severe_index([records[s] for s in severes] + concordant_records)
        hap = compare_haplotypes(mild, severes, recs, sharing, small_keys, sv_positions, args.sv_merge_dist)
        cats = Counter(h["category"] for h in hap.values())

        for vid, rec in recs.items():
            h = hap[vid]
            for g in sorted(rec["genes"]):
                hap_rows.append({
                    "mild": mild, "severe_sibs": ",".join(severes),
                    "#CHROM": rec["chrom"], "POS": rec["pos"], "ID": vid, "REF": rec["ref"], "ALT": rec["alt"],
                    "gene": g, "zygosity": rec["zygosity"], "GT": rec["gt"], "PS": rec["ps"],
                    "evidence_caller": rec["caller"],
                    "parent_of_origin": h["origin"], "origin_source": h["source"], "de_novo": rec["de_novo"],
                    "haplotype_status": h["status"], "present_in_severe": h["present_in_severe"],
                    "category": h["category"], "ibd_vs_severe": h["ibd"],
                })

        # dominant: variants on mild-specific haplotypes (or de novo), in
        # genes where no severe sib / concordant-severe sample carries ANY
        # in-scope variant. One hit is sufficient under this model, so a
        # severe sample with its own hit in the gene rules the gene out --
        # the haplotype filter alone can't do this, since a mild-specific
        # haplotype is by definition one severe didn't inherit, while
        # severe's other haplotype carries its own variants in the gene.
        severe_genes = set()
        for r in [records[s] for s in severes] + concordant_records:
            for srec in r.values():
                severe_genes.update(srec["genes"])
        n_dom = n_dom_gene_severe = 0
        for vid, rec in recs.items():
            h = hap[vid]
            if h["category"] not in ("mild_specific", "de_novo"):
                continue
            for g in rec["genes"]:
                if g in severe_genes:
                    n_dom_gene_severe += 1
                    continue
                key = (rec["chrom"], rec["pos"], rec["ref"], rec["alt"], g)
                info = dom_pass.setdefault(key, {"vid": vid, "gene": g, "rec": rec, "hap": h,
                                                 "reason": h["category"], "samples": set()})
                info["samples"].add(mild)
                dom_gene_samples[g].add(mild)
                n_dom += 1

        # recessive: biallelic in mild, not in severe, not IBD2 with a severe sib
        mild_bi = bi[mild]
        severe_bi = set(concordant_bi).union(*(bi[s] for s in severes))
        n_sev_bi = n_ibd2 = n_rec = 0
        for g, (kind, vids) in mild_bi.items():
            if g in severe_bi:
                n_sev_bi += 1
                continue
            if any(all(is_ibd2(sharing, mild, s, recs[v]["chrom"], recs[v]["pos"]) for v in vids) for s in severes):
                n_ibd2 += 1
                continue
            n_rec += 1
            if kind == "homozygous":
                reason = f"homozygous:{side_summary(vids, hap)}"
            elif kind == "compound_het_parental":
                pat = [v for v in vids if hap[v]["origin"] == "paternal"]
                mat = [v for v in vids if hap[v]["origin"] == "maternal"]
                reason = f"compound_het_parental:pat={side_summary(pat, hap)};mat={side_summary(mat, hap)}"
            else:
                reason = f"compound_het_readphase:{side_summary(vids, hap)}"
            for v in vids:
                rec = recs[v]
                partners = ",".join(x for x in vids if x != v)
                key = (rec["chrom"], rec["pos"], rec["ref"], rec["alt"], g)
                info = rec_pass.setdefault(key, {"vid": v, "gene": g, "rec": rec, "hap": hap[v],
                                                 "reason": reason + (f",partner={partners}" if partners else ""),
                                                 "samples": set()})
                info["samples"].add(mild)
            rec_gene_samples[g][kind].add(mild)

        logger.info(f"{mild} vs {', '.join(severes)}: {len(recs):,} in-scope variants; "
                    + ", ".join(f"{k}={n:,}" for k, n in cats.most_common()))
        logger.info(f"{mild}: dominant {n_dom:,} variant x gene "
                    f"(-{n_dom_gene_severe:,} in genes with a variant in severe); recessive {len(mild_bi):,} biallelic gene(s), "
                    f"-{n_sev_bi} biallelic in severe, -{n_ibd2} IBD2 with severe -> {n_rec}")
        funnel.append({"mild": mild, "severe_sibs": ",".join(severes), "in_scope_variants": len(recs),
                       **{f"cat_{k}": cats.get(k, 0) for k in
                          ("mild_specific", "de_novo", "shared_haplotype", "shared_allele", "unknown_haplotype")},
                       "dominant_dropped_gene_in_severe": n_dom_gene_severe,
                       "dominant_variant_gene_rows": n_dom, "biallelic_genes": len(mild_bi),
                       "dropped_biallelic_in_severe": n_sev_bi, "dropped_ibd2": n_ibd2,
                       "recessive_genes": n_rec})

    # recurrence across independent mild children
    k = args.min_occurrence
    dom_pass = {key: v for key, v in dom_pass.items() if len(dom_gene_samples[v["gene"]]) >= k}
    rec_pass = {key: v for key, v in rec_pass.items()
                if len(set().union(*rec_gene_samples[v["gene"]].values())) >= k}

    gnomad = GnomadNFEAnnotator(args.gnomad_chm13_vcf) if args.gnomad_chm13_vcf else None
    kg = None
    if args.kg_af_vcf and args.kg_af_vcf.rstrip("/").split("/")[-1] != "NO_FILE":
        kg = GnomadNFEAnnotator(args.kg_af_vcf, field="AF", label="1000 Genomes EUR")

    pre = args.output_prefix
    write_model(f"{pre}.dominant_modifier.tsv", dom_pass, gnomad, kg)
    write_model(f"{pre}.recessive_modifier.tsv", rec_pass, gnomad, kg)

    with open(f"{pre}.dominant_modifier.gene_summary.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "N_sibs_with_variant", "samples"])
        for g, ss in sorted(dom_gene_samples.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            if len(ss) >= k:
                w.writerow([g, len(ss), ",".join(sorted(ss))])
    with open(f"{pre}.recessive_modifier.gene_summary.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "N_sibs_with_variant", "N_homo_sibs", "N_compound_het_sibs", "samples"])
        rows = []
        for g, kinds in rec_gene_samples.items():
            allss = set().union(*kinds.values())
            if len(allss) < k:
                continue
            homo = kinds.get("homozygous", set())
            ch = kinds.get("compound_het_parental", set()) | kinds.get("compound_het_readphase", set())
            rows.append([g, len(allss), len(homo), len(ch), ",".join(sorted(allss))])
        for r in sorted(rows, key=lambda r: (-r[1], r[0])):
            w.writerow(r)

    hap_cols = ["mild", "severe_sibs", "#CHROM", "POS", "ID", "REF", "ALT", "gene", "zygosity", "GT", "PS",
                "evidence_caller", "parent_of_origin", "origin_source", "de_novo", "haplotype_status",
                "present_in_severe", "category", "ibd_vs_severe"]
    with open(f"{pre}.haplotype_variants.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hap_cols, delimiter="\t")
        w.writeheader()
        for r in sorted(hap_rows, key=lambda r: (r["mild"], r["#CHROM"], r["POS"], r["gene"])):
            w.writerow(r)
    with open(f"{pre}.haplotype_funnel.tsv", "w", newline="") as f:
        cols = list(funnel[0]) if funnel else ["mild"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(funnel)
    logger.info(f"final: {len(dom_pass):,} dominant and {len(rec_pass):,} recessive variant x gene rows "
                f"(>= {k} mild children)")


def write_model(path, passing, gnomad, kg):
    # AF comes from the VCF's INFO (bin/annotate_af.py, run once before this
    # step); a population file passed explicitly overrides it (legacy runs
    # whose VCFs weren't AF-annotated)
    af = annotate_final_pass(passing, gnomad) if passing and gnomad else \
        {key: info["rec"].get("gnomad_af") for key, info in passing.items()}
    kgaf = annotate_final_pass(passing, kg) if passing and kg else \
        {key: info["rec"].get("kg_af") for key, info in passing.items()}
    rows = []
    for key, info in passing.items():
        rec, h = info["rec"], info["hap"]
        rows.append({
            "#CHROM": rec["chrom"], "POS": rec["pos"], "ID": info["vid"], "REF": rec["ref"], "ALT": rec["alt"],
            "gene": info["gene"], "zygosity": rec["zygosity"], "reason": info["reason"],
            "evidence_caller": rec["caller"], "callers_detected": ",".join(rec["callers_all"]),
            "GT": rec["gt"], "PS": rec["ps"], "samples": ",".join(sorted(info["samples"])),
            "N_samples": len(info["samples"]),
            "gnomad_AF_NFE": "" if af.get(key) is None else af.get(key),
            "kg_AF_EUR": "" if kgaf.get(key) is None else kgaf.get(key),
            "parent_of_origin": h["origin"], "origin_source": h["source"], "haplotype_status": h["status"],
            "category": h["category"], "ibd_vs_severe": h["ibd"],
        })
    rows.sort(key=lambda r: (r["#CHROM"], r["POS"], r["gene"]))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=BASE_COLS + HAP_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    logger.info(f"wrote {len(rows):,} rows -> {path}")


def get_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-family_json", required=True)
    p.add_argument("-vcf_manifest", required=True,
                   help="childID<TAB>consolidated+snpEff VCF (with INFO/ORIGIN), one per child")
    p.add_argument("-sib_parent_segments", required=True, help="sib_haplotype_sharing.py .sib_parent_segments.tsv")
    p.add_argument("-sv_merge_dist", type=int, default=500)
    p.add_argument("-min_occurrence", type=int, default=1)
    p.add_argument("-gnomad_chm13_vcf", default=None,
                   help="Only for VCFs without INFO/gnomAD_AF_NFE (annotate_af.py); otherwise AF is read from the VCF")
    p.add_argument("-kg_af_vcf", default=None,
                   help="Only for VCFs without INFO/KG_AF_EUR (annotate_af.py); otherwise AF is read from the VCF")
    p.add_argument("-output_prefix", default="cohort")
    return p.parse_args()


if __name__ == "__main__":
    run(get_args())
