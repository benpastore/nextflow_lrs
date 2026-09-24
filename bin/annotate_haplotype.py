#!/usr/bin/env python3
"""
Post-hoc haplotype annotation for callers that have no native phasing mode
(Straglr repeat expansions, Spectre CNVs) -- unlike Clair3 (phased via
longphase/longphase_sv) and Sniffles (phased natively via --phase), neither
tool reads the HP:i:1/HP:i:2 tags already present on reads in the haplotagged
BAM they were called against.

Rather than re-calling on haplotype-split BAMs (which halves read depth and
hurts sensitivity, especially for large Spectre CNVs), this tallies HP tags
on reads near each call's boundaries in the existing haplotagged BAM and
stamps a haplotype assignment onto the call without re-genotyping it. Only
the flanking windows around start/end are sampled, not the whole call span --
for a multi-kb/Mb CNV the breakpoint-flanking reads are what's actually
informative about haplotype; reads deep inside the call carry no more signal
and would be expensive to pull for a large region.

Two input modes, since the two callers have different table layouts:
  straglr: TSV, chrom/start/end in the first three columns (straglr's own
           locus-report convention), optional "#"-prefixed header line.
  spectre: BED, chrom/start/end in the first three columns, no header.

Either way, output is the same input table with one "haplotype" column
appended (HP1 / HP2 / AMBIGUOUS).
"""

import argparse
import subprocess


def count_hp_reads(bam, chrom, pos, flank, hp_value):
    start = max(1, pos - flank)
    end = pos + flank
    region = f"{chrom}:{start}-{end}"
    # capture_output=/text= need Python >=3.7 -- the rnaseq container's
    # `rnaseq` conda env runs 3.6, so use the equivalent stdout/stderr pipes.
    result = subprocess.run(
        ["samtools", "view", "-c", "-d", f"HP:{hp_value}", bam, region],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True, check=True,
    )
    return int(result.stdout.strip())


def assign_haplotype(bam, chrom, start, end, flank, majority_ratio):
    hp1 = count_hp_reads(bam, chrom, start, flank, 1) + count_hp_reads(bam, chrom, end, flank, 1)
    hp2 = count_hp_reads(bam, chrom, start, flank, 2) + count_hp_reads(bam, chrom, end, flank, 2)
    total = hp1 + hp2

    if total == 0:
        return "AMBIGUOUS"

    if hp1 / total >= majority_ratio:
        return "HP1"
    if hp2 / total >= majority_ratio:
        return "HP2"
    return "AMBIGUOUS"


def iter_calls(path, mode):
    """Yield (raw_line, chrom, start, end) for each call, skipping/passing through header lines."""
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                yield (line, None, None, None)
                continue
            fields = line.split("\t")
            if mode == "straglr" and len(fields) < 3:
                yield (line, None, None, None)
                continue
            chrom, start, end = fields[0], fields[1], fields[2]
            try:
                yield (line, chrom, int(start), int(end))
            except ValueError:
                # non-numeric start/end -- treat as a header we didn't
                # recognize by "#" and pass it through unannotated.
                yield (line, None, None, None)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True, help="Haplotagged BAM (HP:i:1/HP:i:2 tags), indexed")
    p.add_argument("--calls", required=True, help="Straglr TSV or Spectre BED (chrom/start/end in first 3 columns)")
    p.add_argument("--mode", required=True, choices=["straglr", "spectre"])
    p.add_argument("--output", required=True)
    p.add_argument("--flank", type=int, default=1000, help="bp window around each call's start/end to tally HP tags in")
    p.add_argument("--majority-ratio", type=float, default=0.8,
                    help="Fraction of HP-tagged reads one haplotype must reach to call HP1/HP2 vs AMBIGUOUS")
    args = p.parse_args()

    with open(args.output, "w") as out:
        for line, chrom, start, end in iter_calls(args.calls, args.mode):
            if chrom is None:
                if args.mode == "straglr" and line.startswith("#"):
                    out.write(f"{line}\thaplotype\n")
                else:
                    out.write(f"{line}\n")
                continue

            haplotype = assign_haplotype(args.bam, chrom, start, end, args.flank, args.majority_ratio)
            out.write(f"{line}\t{haplotype}\n")


if __name__ == "__main__":
    main()
