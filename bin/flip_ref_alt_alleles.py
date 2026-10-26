#!/usr/bin/env python3
"""
Recover CrossMap's REF==ALT liftover records (run with --no-comp-alleles).

At these sites the hg38 reference base IS the candidate's ALT allele --
CHM13 and hg38 just carry different alleles there. CrossMap rewrites REF
from hg38, so REF ends up equal to ALT and the record describes no change.
Rewrite it as hg38 REF (= the candidate allele) > original CHM13 REF, and
tag INFO ALLELE_FLIPPED so downstream (run_alphagenome.py,
build_variant_report.py) knows the scored change runs opposite to the
candidate allele.

Strand: CrossMap reverse-complements ALT on minus-strand chain hits, so for
an SNV the lifted ALT differs from the original ALT exactly when the hit is
on the minus strand (a base is never its own complement); the original REF
is complemented to match in that case. Non-SNV REF==ALT records can't be
oriented this way -- they're dropped to the .unmap file as
Fail(REF==ALT_unflippable) (not expected in practice: REF==ALT needs
equal-length alleles).

Stdlib only -- runs in the crossmap container.
"""
import argparse
import sys

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

FLIPPED_HEADER = ('##INFO=<ID=ALLELE_FLIPPED,Number=0,Type=Flag,Description="hg38 reference base is the '
                  'candidate (CHM13 ALT) allele; record rewritten as hg38 REF > CHM13 REF, i.e. the scored '
                  'change runs opposite to the candidate allele">\n')


def read_original(path):
    """ID -> (REF, ALT) from the pre-liftover candidate VCF."""
    alleles = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            alleles[fields[2]] = (fields[3], fields[4])
    return alleles


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-original_vcf", required=True, help="Pre-liftover candidate VCF (ID = original chrom_pos_ref_alt)")
    p.add_argument("-lifted_vcf", required=True, help="CrossMap --no-comp-alleles output")
    p.add_argument("-unmap", required=True, help="CrossMap .unmap file to append unflippable records to")
    p.add_argument("-output", required=True)
    args = p.parse_args()

    original = read_original(args.original_vcf)
    n_kept = n_flipped = n_dropped = 0

    with open(args.lifted_vcf) as f, open(args.output, "w") as out, open(args.unmap, "a") as unmap:
        for line in f:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                out.write(FLIPPED_HEADER)
                out.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            ref, alt = fields[3], fields[4]
            if ref.upper() != alt.upper():
                out.write(line)
                n_kept += 1
                continue

            orig_ref, orig_alt = original.get(fields[2], ("", ""))
            if len(ref) != 1 or len(orig_ref) != 1 or len(orig_alt) != 1:
                unmap.write(line.rstrip("\n") + "\tFail(REF==ALT_unflippable)\n")
                n_dropped += 1
                continue

            minus = alt.upper() != orig_alt.upper()
            fields[4] = orig_ref.translate(COMPLEMENT) if minus else orig_ref
            fields[7] = "ALLELE_FLIPPED" if fields[7] in ("", ".") else fields[7] + ";ALLELE_FLIPPED"
            out.write("\t".join(fields) + "\n")
            n_flipped += 1

    sys.stderr.write(
        f"[flip_ref_alt_alleles] {n_kept} lifted as-is, {n_flipped} REF==ALT record(s) flipped "
        f"(hg38 ref = candidate allele), {n_dropped} unflippable non-SNV REF==ALT record(s) sent to {args.unmap}\n"
    )


if __name__ == "__main__":
    main()
