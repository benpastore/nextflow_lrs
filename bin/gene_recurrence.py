#!/usr/bin/env python3
"""
Count, per gene, how many mild children have it as a passing modifier gene
(modifier_variants.py <mild>.passing_genes.tsv), then, for the recurrent
genes (passing in >= -min_mild mild children), write EVERY variant in the
gene carried by those mild children or their severe sib(s)
(modifier_variants.py <mild>.gene_variants.tsv.gz) as two CHM13 VCFs for
AlphaGenome / AlphaMissense:

  candidate_variants.mild.chm13.vcf    variants carried by >= 1 such mild child
  candidate_variants.severe.chm13.vcf  variants carried by >= 1 such severe sib

Each variant appears once per VCF however many families carry it (no
rescoring across families). A variant carried on both sides is in both
VCFs. Only families where the gene passes contribute -- a severe sib's
variants in gene G are included only if G passes for its mild sib.

INFO (both VCFs; ID = chrom_pos_ref_alt as run_alphagenome.py /
run_alphamissense.py expect):
  GENE            recurrent gene(s) the variant lies in / is linked to
  SOURCE_MODEL    rule(s) the variant supports in a mild child
                  (homozygous / compound_het / regulatory), or 'context'
                  (in the gene, but not one of the supporting variants)
  REASON          family:gene:rule-or-context entries, pipe-separated
  MILD_SAMPLES    mild carriers, sample:zygosity, comma-separated ('.' = none)
  SEVERE_SAMPLES  severe carriers, sample:zygosity ('.' = none)
  SAMPLES         this VCF's side only (= MILD_SAMPLES or SEVERE_SAMPLES)
  CARRIERS        per family: family:mild_only / severe_only / shared /
                  shared_diff (see modifier_variants.py)

Symbolic-ALT SVs (<DEL>, <INS>, ...) are left out -- neither tool can score
them -- and counted in the log. With -max_af, a variant is also left out
when every annotated AF (gnomAD NFE, 1000 Genomes EUR) is >= max_af
(unannotated = kept); the gene counts never use AF.

Outputs:
  cohort.passing_genes.tsv              every family's passing-gene rows, concatenated
  cohort.gene_recurrence.tsv            one row per gene
  cohort.candidate_gene_variants.tsv.gz the gene_variants rows of the recurrent genes
  candidate_variants.{mild,severe}.chm13.vcf
"""
import argparse
import csv
import gzip
import sys
from collections import defaultdict

REASONS = ("homozygous", "compound_het", "regulatory")
SIDES = ("mild", "severe")


def log(msg):
    sys.stderr.write(f"[gene_recurrence] {msg}\n")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_af(v):
    try:
        af = float(v)
    except (TypeError, ValueError):
        return None
    return None if af != af else af


def chrom_key(c):
    c = c[3:] if c.startswith("chr") else c
    return (0, int(c), "") if c.isdigit() else (1, 0, c)


def fmt_carriers(d):
    return ",".join(f"{s}:{z}" for s, z in sorted(d.items())) or "."


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-passing_genes", nargs="+", required=True, help="<mild>.passing_genes.tsv files")
    p.add_argument("-gene_variants", nargs="+", required=True, help="<mild>.gene_variants.tsv.gz files")
    p.add_argument("-min_mild", type=int, default=2, help="a gene is recurrent when it passes in >= this many mild children")
    p.add_argument("-max_af", type=float, default=None,
                   help="leave a variant out of the AlphaGenome/AlphaMissense VCFs when every annotated AF is >= this")
    args = p.parse_args()

    # ---- gene counts
    header = None
    rows = []
    for path in sorted(args.passing_genes):
        with open(path) as f:
            reader = csv.reader(f, delimiter="\t")
            h = next(reader, None)
            if h is None:
                continue
            header = header or h
            rows.extend(reader)
    if header is None:
        sys.exit("[gene_recurrence] no passing-gene rows")
    col = {h: i for i, h in enumerate(header)}
    with open("cohort.passing_genes.tsv", "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    fams = defaultdict(set)
    milds = defaultdict(set)
    by_reason = defaultdict(lambda: defaultdict(set))
    for r in rows:
        g, fam, mild = r[col["gene"]], r[col["family"]], r[col["mild"]]
        fams[g].add(fam)
        milds[g].add(mild)
        for reason in r[col["reasons"]].split(","):
            if reason:
                by_reason[g][reason].add(mild)
    recurrent = {g for g, m in milds.items() if len(m) >= args.min_mild}
    with open("cohort.gene_recurrence.tsv", "w") as out:
        out.write("gene\tn_mild\tn_families\tmild_samples\tfamilies\t"
                  + "\t".join(f"n_mild_{x}" for x in REASONS) + f"\trecurrent_ge{args.min_mild}\n")
        for g in sorted(milds, key=lambda g: (-len(milds[g]), g)):
            out.write(f"{g}\t{len(milds[g])}\t{len(fams[g])}\t{','.join(sorted(milds[g]))}\t"
                      f"{','.join(sorted(fams[g]))}\t"
                      + "\t".join(str(len(by_reason[g][x])) for x in REASONS)
                      + f"\t{'yes' if g in recurrent else 'no'}\n")
    log(f"{len(milds):,} passing genes across {len({r[col['mild']] for r in rows})} mild children; "
        f"{len(recurrent):,} in >= {args.min_mild} mild children")

    # ---- every variant (mild or severe) in the recurrent genes
    variants = {}   # (chrom,pos,ref,alt) -> entry
    n_symbolic = set()
    n_common = set()
    out_header = None
    with gzip.open("cohort.candidate_gene_variants.tsv.gz", "wt", compresslevel=4) as tbl:
        for path in sorted(args.gene_variants):
            with opener(path) as f:
                reader = csv.DictReader(f, delimiter="\t")
                zcols = {side: [c for c in (reader.fieldnames or []) if c.startswith(f"{side}_zyg_")] for side in SIDES}
                for r in reader:
                    if r["gene"] not in recurrent:
                        continue
                    # per-family files have per-sample columns -> long table: one row per (gene, variant, family)
                    if out_header is None:
                        out_header = ["family", "mild", "gene", "gene_reasons", "chrom", "pos", "id", "ref", "alt",
                                      "svtype", "end", "carriers", "mild_carriers", "severe_carriers", "origin",
                                      "mild_supporting", "gnomAD_AF_NFE", "KG_AF_EUR"]
                        tbl.write("\t".join(out_header) + "\n")
                    carriers = {side: {c[len(side) + 5:]: r[c] for c in zcols[side] if r[c] != "absent"}
                                for side in SIDES}
                    tbl.write("\t".join([r[k] for k in out_header[:12]]
                                        + [fmt_carriers(carriers["mild"]), fmt_carriers(carriers["severe"])]
                                        + [r[k] for k in out_header[14:]]) + "\n")
                    key = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                    if r["alt"].startswith("<") or r["alt"] in (".", "*"):
                        n_symbolic.add(key)
                        continue
                    known = [a for a in (parse_af(r["gnomAD_AF_NFE"]), parse_af(r["KG_AF_EUR"])) if a is not None]
                    if args.max_af is not None and known and min(known) >= args.max_af:
                        n_common.add(key)
                        continue
                    e = variants.setdefault(key, {"models": set(), "genes": set(), "detail": set(),
                                                  "mild": {}, "severe": {}, "carriers": set()})
                    sup = [x for x in r["mild_supporting"].split(",") if x and x != "."]
                    e["models"].update(sup or ["context"])
                    e["genes"].add(r["gene"])
                    e["detail"].add(f"{r['family']}:{r['gene']}:{'+'.join(sup) or 'context'}")
                    e["carriers"].add(f"{r['family']}:{r['carriers']}")
                    for side in SIDES:
                        e[side].update(carriers[side])

    n_out = {}
    for side in SIDES:
        path = f"candidate_variants.{side}.chm13.vcf"
        keys = [k for k in variants if variants[k][side]]
        n_out[side] = len(keys)
        with open(path, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write(f"##candidate_side={side}\n")
            out.write('##INFO=<ID=SOURCE_MODEL,Number=.,Type=String,Description="Modifier rule(s) this variant supports in a mild child (homozygous / compound_het / regulatory), or context">\n')
            out.write('##INFO=<ID=GENE,Number=.,Type=String,Description="Recurrent modifier gene(s) this variant lies in / is linked to">\n')
            out.write('##INFO=<ID=REASON,Number=.,Type=String,Description="family:gene:rule-or-context entries, pipe-separated">\n')
            out.write('##INFO=<ID=MILD_SAMPLES,Number=.,Type=String,Description="Mild carriers, sample:zygosity (. = none)">\n')
            out.write('##INFO=<ID=SEVERE_SAMPLES,Number=.,Type=String,Description="Severe carriers, sample:zygosity (. = none)">\n')
            out.write(f'##INFO=<ID=SAMPLES,Number=.,Type=String,Description="{side.capitalize()} carriers (this VCF\'s side), sample:zygosity">\n')
            out.write('##INFO=<ID=CARRIERS,Number=.,Type=String,Description="family:mild_only/severe_only/shared/shared_diff">\n')
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for key in sorted(keys, key=lambda k: (chrom_key(k[0]), k[1], k[2], k[3])):
                chrom, pos, ref, alt = key
                e = variants[key]
                info = ";".join([f"SOURCE_MODEL={','.join(sorted(e['models']))}",
                                 f"GENE={','.join(sorted(e['genes']))}",
                                 f"REASON={'|'.join(sorted(e['detail']))}",
                                 f"MILD_SAMPLES={fmt_carriers(e['mild'])}",
                                 f"SEVERE_SAMPLES={fmt_carriers(e['severe'])}",
                                 f"SAMPLES={fmt_carriers(e[side])}",
                                 f"CARRIERS={','.join(sorted(e['carriers']))}"])
                out.write(f"{chrom}\t{pos}\t{chrom}_{pos}_{ref}_{alt}\t{ref}\t{alt}\t.\t.\t{info}\n")
    both = sum(1 for e in variants.values() if e["mild"] and e["severe"])
    log(f"{len(variants):,} unique variant(s) in recurrent genes: {n_out['mild']:,} -> mild VCF, "
        f"{n_out['severe']:,} -> severe VCF ({both:,} in both); "
        f"{len(n_symbolic):,} symbolic-ALT SV(s) left out"
        + (f"; {len(n_common):,} left out with every annotated AF >= {args.max_af}" if args.max_af is not None else ""))


if __name__ == "__main__":
    main()
