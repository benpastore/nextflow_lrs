#!/usr/bin/env python3
"""
Union the recessive- and dominant-modifier scripts' passing-variant tables
(bin/recessive_modifier_consolidated.py, bin/dominant_modifier_consolidated.py --
both emit the same #CHROM/POS/ID/REF/ALT/gene/zygosity/reason/samples schema)
into one small CHM13-coordinate VCF for AlphaGenome scoring, carrying full
traceability (which model(s), which gene(s)/reason(s)/sample(s)) through an
INFO field so it survives liftover and reattaches to the final output.

Symbolic-ALT structural variants (<DEL>, <INS>, ...) are excluded -- both
CrossMap liftover and AlphaGenome's variant scoring expect explicit ref/alt
sequence variants, not SV notation. Dropped count is logged, not silent.
"""
import argparse
import csv
import sys


def read_modifier_tsv(path, model_name):
    """Yields (chrom, pos, ref, alt, gene, reason, samples) per row, skipping
    symbolic-ALT SV rows."""
    if not path:
        return
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            alt = row["ALT"]
            if alt.startswith("<"):
                continue
            yield (
                row["#CHROM"], int(row["POS"]), row["REF"], alt,
                row["gene"], row["reason"], row["samples"],
            )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-recessive_tsv", type=str, default=None,
                    help="cohort.recessive_modifier.tsv (optional)")
    p.add_argument("-dominant_tsv", type=str, default=None,
                    help="cohort.dominant_modifier.tsv (optional)")
    p.add_argument("-output", type=str, required=True, help="Output VCF path")
    args = p.parse_args()

    if not args.recessive_tsv and not args.dominant_tsv:
        sys.exit("build_alphagenome_input.py: at least one of -recessive_tsv/-dominant_tsv is required")

    variants = {}  # (chrom,pos,ref,alt) -> {models:set, genes:set, reasons:set, samples:set}
    n_sv_skipped = 0

    for path, model_name in ((args.recessive_tsv, "recessive"), (args.dominant_tsv, "dominant")):
        if not path:
            continue
        for chrom, pos, ref, alt, gene, reason, samples in read_modifier_tsv(path, model_name):
            key = (chrom, pos, ref, alt)
            entry = variants.setdefault(key, {"models": set(), "genes": set(), "reasons": set(), "samples": set()})
            entry["models"].add(model_name)
            entry["genes"].add(gene)
            entry["reasons"].add(reason)
            entry["samples"].update(s for s in samples.split(",") if s)

    # count SV rows skipped separately, by re-scanning (cheap; these files are small)
    for path in (args.recessive_tsv, args.dominant_tsv):
        if not path:
            continue
        with open(path, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row["ALT"].startswith("<"):
                    n_sv_skipped += 1

    with open(args.output, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##INFO=<ID=SOURCE_MODEL,Number=.,Type=String,Description="Which modifier model(s) this variant passed under: recessive, dominant, or both">\n')
        out.write('##INFO=<ID=GENE,Number=.,Type=String,Description="Gene(s) this variant qualified through">\n')
        out.write('##INFO=<ID=REASON,Number=.,Type=String,Description="Modifier-script reason string(s), pipe-separated">\n')
        out.write('##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Mild sample ID(s) this variant passed in">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for (chrom, pos, ref, alt) in sorted(variants, key=lambda k: (k[0], k[1])):
            entry = variants[(chrom, pos, ref, alt)]
            vid = f"{chrom}_{pos}_{ref}_{alt}"
            info = ";".join([
                f"SOURCE_MODEL={','.join(sorted(entry['models']))}",
                f"GENE={','.join(sorted(entry['genes']))}",
                f"REASON={'|'.join(sorted(entry['reasons']))}",
                f"SAMPLES={','.join(sorted(entry['samples']))}",
            ])
            out.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{info}\n")

    sys.stderr.write(
        f"[build_alphagenome_input] {len(variants)} unique candidate variant(s) written to {args.output} "
        f"({n_sv_skipped} symbolic-ALT SV row(s) skipped -- not scored by this pass)\n"
    )


if __name__ == "__main__":
    main()
