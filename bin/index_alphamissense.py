#!/usr/bin/env python3
"""
One-time prep: turn DeepMind's AlphaMissense_hg38.tsv.gz (plain gzip, ~71M
rows) into a bgzipped + tabix-indexed copy, so run_alphamissense.py can
fetch each candidate position directly instead of streaming the whole file.

Rows are split per chromosome into temp files while streaming (tabix needs
each chromosome contiguous and positions ascending); a chromosome is only
sorted in memory if its rows turn out to be out of order. The comment /
#CHROM header lines are kept, so the column names travel with the index.
Works for any of the AlphaMissense position-level tables (hg38, isoforms
hg38) since nothing here depends on columns past CHROM/POS.

    index_alphamissense.py -input AlphaMissense_hg38.tsv.gz -output AlphaMissense_hg38.tabix.tsv
      -> AlphaMissense_hg38.tabix.tsv.gz + .tbi

Needs pysam (a CrossMap dependency -- runs in the crossmap container).
"""
import argparse
import gzip
import os
import shutil
import sys
import tempfile

import pysam


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-input", required=True, help="AlphaMissense_*.tsv.gz as downloaded (gzip or plain)")
    p.add_argument("-output", required=True, help="Uncompressed output path; becomes <output>.gz + <output>.gz.tbi")
    args = p.parse_args()

    tmpdir = tempfile.mkdtemp(prefix="am_split_", dir=os.path.dirname(os.path.abspath(args.output)))
    header = []
    chrom_order = []
    handles = {}
    last_pos = {}
    unsorted = set()
    n = 0

    opener = gzip.open if args.input.endswith(".gz") else open
    with opener(args.input, "rt") as f:
        for line in f:
            if line.startswith("#"):
                header.append(line)
                continue
            chrom, pos, _ = line.split("\t", 2)
            pos = int(pos)
            if chrom not in handles:
                chrom_order.append(chrom)
                handles[chrom] = open(os.path.join(tmpdir, f"{len(chrom_order)}.tsv"), "w")
                last_pos[chrom] = pos
            elif pos < last_pos[chrom]:
                unsorted.add(chrom)
            last_pos[chrom] = max(last_pos[chrom], pos)
            handles[chrom].write(line)
            n += 1
    for h in handles.values():
        h.close()

    if not any(l.startswith("#CHROM") for l in header):
        sys.exit(f"index_alphamissense.py: no #CHROM header line in {args.input} -- not an AlphaMissense position table?")

    with open(args.output, "w") as out:
        out.writelines(header)
        for i, chrom in enumerate(chrom_order, start=1):
            part = os.path.join(tmpdir, f"{i}.tsv")
            if chrom in unsorted:
                with open(part) as f:
                    rows = f.readlines()
                rows.sort(key=lambda l: int(l.split("\t", 2)[1]))
                out.writelines(rows)
            else:
                with open(part) as f:
                    shutil.copyfileobj(f, out)
            os.unlink(part)
    os.rmdir(tmpdir)

    # compresses <output> to <output>.gz (bgzf) and removes the plain copy
    pysam.tabix_index(args.output, seq_col=0, start_col=1, end_col=1, meta_char="#", force=True)
    sys.stderr.write(
        f"[index_alphamissense] {n:,} rows over {len(chrom_order)} chromosome(s)"
        f"{f' ({len(unsorted)} re-sorted)' if unsorted else ''} -> {args.output}.gz (+ .tbi)\n"
    )


if __name__ == "__main__":
    main()
