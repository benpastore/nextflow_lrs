#!/usr/bin/env python3
"""
Recover CrossMap's Fail(REF==ALT) liftover records.

At these sites the hg38 reference base IS the candidate's ALT allele --
CHM13 and hg38 just carry different alleles there. CrossMap (0.7.4)
rewrites REF from hg38, finds ALT == REF, and sends the record to .unmap --
unconditionally: it drops ALT alleles equal to the new REF *before* its
--no-comp-alleles option is consulted, so that flag can't keep them. The
.unmap line is the original CHM13 record, without the lifted coordinates.

So the coordinates are lifted separately with `CrossMap bed` (pure
coordinate liftover, no allele logic -- see LIFTOVER_ALPHAGENOME_VARIANTS),
and each recovered record is rewritten as
    REF = hg38 base (= the candidate allele), ALT = original CHM13 REF
both complemented when the chain hit is on the minus strand, and tagged
INFO ALLELE_FLIPPED so downstream (run_alphagenome.py,
build_variant_report.py) knows the scored change runs opposite to the
candidate allele.

Only SNVs with a single, 1bp bed hit are recovered; when -reference is
given, the hg38 base at the new position must also equal the expected REF.
Anything else stays in the .unmap (multi-hit -> Fail(REF==ALT_multi_hit),
base mismatch -> Fail(REF==ALT_ref_mismatch), non-SNV unchanged).

Also handles REF==ALT records that do reach the lifted VCF (other CrossMap
versions honour --no-comp-alleles), the same way. Needs only the stdlib,
plus pysam (a CrossMap dependency) for the optional reference check.
"""
import argparse
import os
import sys
from collections import defaultdict

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

FLIPPED_HEADER = ('##INFO=<ID=ALLELE_FLIPPED,Number=0,Type=Flag,Description="hg38 reference base is the '
                  'candidate (CHM13 ALT) allele; record rewritten as hg38 REF > CHM13 REF, i.e. the scored '
                  'change runs opposite to the candidate allele">\n')


def comp(base, minus):
    return base.translate(COMPLEMENT) if minus else base


def add_flag(info):
    return "ALLELE_FLIPPED" if info in ("", ".") else info + ";ALLELE_FLIPPED"


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


def read_lifted_bed(path):
    """ID -> [(chrom, start, end, strand)] from `CrossMap bed` output of
    chrom/start/end/ID/score/strand records (several entries = split hit)."""
    hits = defaultdict(list)
    if path and os.path.exists(path):
        with open(path) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 6:
                    hits[fields[3]].append((fields[0], int(fields[1]), int(fields[2]), fields[5]))
    return hits


def open_reference(path):
    if not path:
        return None
    try:
        import pysam
    except ImportError:
        sys.stderr.write("[flip_ref_alt_alleles] pysam unavailable -- skipping the hg38 base check\n")
        return None
    return pysam.FastaFile(path)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-original_vcf", required=True, help="Pre-liftover candidate VCF (ID = original chrom_pos_ref_alt)")
    p.add_argument("-lifted_vcf", required=True, help="CrossMap vcf output")
    p.add_argument("-unmap_in", required=True, help="CrossMap vcf .unmap (original line + tab + Fail(reason))")
    p.add_argument("-lifted_bed", default="", help="CrossMap bed output for the Fail(REF==ALT) SNV positions")
    p.add_argument("-reference", default="", help="hg38 FASTA (+.fai) to confirm the recovered REF base")
    p.add_argument("-output", required=True, help="Lifted VCF incl. recovered records")
    p.add_argument("-unmap_out", required=True, help="Remaining unmapped records")
    args = p.parse_args()

    original = read_original(args.original_vcf)
    bed_hits = read_lifted_bed(args.lifted_bed)
    ref_fa = open_reference(args.reference)
    n_kept = n_flipped_vcf = n_recovered = 0
    unmap_reasons = defaultdict(int)

    with open(args.output, "w") as out, open(args.unmap_out, "w") as unmap:
        # pass 1: CrossMap's lifted records
        with open(args.lifted_vcf) as f:
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
                    unmap_reasons["Fail(REF==ALT_unflippable)"] += 1
                    continue
                # CrossMap complements ALT on minus-strand hits, so the
                # lifted ALT differs from the original exactly then
                minus = alt.upper() != orig_alt.upper()
                fields[4] = comp(orig_ref, minus)
                fields[7] = add_flag(fields[7])
                out.write("\t".join(fields) + "\n")
                n_flipped_vcf += 1

        # pass 2: Fail(REF==ALT) records from the .unmap, via bed coordinates
        with open(args.unmap_in) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                reason = fields[-1]
                rec = fields[:-1]
                if reason != "Fail(REF==ALT)" or len(rec[3]) != 1 or len(rec[4]) != 1:
                    unmap.write(line if line.endswith("\n") else line + "\n")
                    unmap_reasons[reason] += 1
                    continue

                hits = bed_hits.get(rec[2], [])
                if len(hits) != 1 or hits[0][2] - hits[0][1] != 1:
                    unmap.write("\t".join(rec) + "\tFail(REF==ALT_multi_hit)\n")
                    unmap_reasons["Fail(REF==ALT_multi_hit)"] += 1
                    continue
                chrom, start, _end, strand = hits[0]
                minus = strand == "-"
                new_ref = comp(rec[4].upper(), minus)   # hg38 base = candidate allele
                new_alt = comp(rec[3].upper(), minus)   # original CHM13 REF
                if ref_fa is not None:
                    try:
                        hg38_base = ref_fa.fetch(chrom, start, start + 1).upper()
                    except (KeyError, ValueError):
                        hg38_base = ""
                    if hg38_base != new_ref:
                        unmap.write("\t".join(rec) + "\tFail(REF==ALT_ref_mismatch)\n")
                        unmap_reasons["Fail(REF==ALT_ref_mismatch)"] += 1
                        continue

                info = rec[7] if len(rec) > 7 else "."
                out.write("\t".join([chrom, str(start + 1), rec[2], new_ref, new_alt,
                                     rec[5] if len(rec) > 5 else ".", rec[6] if len(rec) > 6 else ".",
                                     add_flag(info)]) + "\n")
                n_recovered += 1

    remaining = ", ".join(f"{n} {r}" for r, n in sorted(unmap_reasons.items())) or "none"
    sys.stderr.write(
        f"[flip_ref_alt_alleles] {n_kept} lifted as-is, {n_flipped_vcf + n_recovered} REF==ALT record(s) flipped "
        f"(hg38 ref = candidate allele; {n_recovered} recovered from .unmap via bed liftover); "
        f"still unmapped: {remaining}\n"
    )


if __name__ == "__main__":
    main()
