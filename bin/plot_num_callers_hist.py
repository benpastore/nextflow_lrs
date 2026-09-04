#!/usr/bin/env python3
"""
Summarize and plot the NUM_CALLERS distribution from one or more
consolidate_variants.py master VCFs.

Usage:
    plot_num_callers_hist.py S1.master.vcf [S2.master.vcf ...] -o hist.png

    # or a directory of *.master.vcf(.gz) files:
    plot_num_callers_hist.py --input-dir out/ -o hist.png

Prints a count table (NUM_CALLERS -> number of sites) to stdout, and
writes a histogram PNG.
"""

import argparse
import gzip
import os
import re
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def opener(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


NUM_CALLERS_RE = re.compile(r"(?:^|;)NUM_CALLERS=([0-9]+)")


def count_num_callers(path):
    counts = Counter()
    with opener(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = f[7]
            m = NUM_CALLERS_RE.search(info)
            if not m:
                continue
            counts[int(m.group(1))] += 1
    return counts


def gather_inputs(args):
    if args.input_dir:
        paths = []
        for root, _, files in os.walk(args.input_dir):
            for fname in files:
                if fname.endswith(".master.vcf") or fname.endswith(".master.vcf.gz"):
                    paths.append(os.path.join(root, fname))
        if not paths:
            sys.exit(f"No *.master.vcf(.gz) files found under {args.input_dir}")
        return sorted(paths)
    return args.vcf


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf", nargs="*", help="Master VCF path(s) (.vcf or .vcf.gz)")
    p.add_argument("--input-dir", help="Directory to scan for *.master.vcf(.gz) files")
    p.add_argument("-o", "--output", default="num_callers_hist.png", help="Output PNG path")
    p.add_argument("--per-sample", action="store_true", help="Plot one bar group per input file instead of pooling all sites together")
    args = p.parse_args()

    paths = gather_inputs(args)
    if not paths:
        p.error("provide VCF path(s) or --input-dir")

    per_file = {}
    for path in paths:
        sample = os.path.basename(path).split(".master.vcf")[0]
        per_file[sample] = count_num_callers(path)

    total = Counter()
    for c in per_file.values():
        total.update(c)

    max_n = max(total) if total else 0
    print(f"{'NUM_CALLERS':>12}  {'SITES':>10}  {'PCT':>6}")
    grand_total = sum(total.values())
    for n in range(1, max_n + 1):
        n_sites = total.get(n, 0)
        pct = 100 * n_sites / grand_total if grand_total else 0
        print(f"{n:>12}  {n_sites:>10}  {pct:5.1f}%")
    print(f"{'TOTAL':>12}  {grand_total:>10}")
    multi = sum(v for k, v in total.items() if k > 1)
    print(f"\n{multi} / {grand_total} sites ({100*multi/grand_total:.1f}%) found by more than one caller"
          if grand_total else "\nNo sites found.")

    fig, ax = plt.subplots(figsize=(6, 4))
    if args.per_sample and len(per_file) > 1:
        samples = list(per_file.keys())
        bins = list(range(1, max_n + 1))
        width = 0.8 / len(samples)
        for i, sample in enumerate(samples):
            heights = [per_file[sample].get(n, 0) for n in bins]
            xs = [b + i * width for b in bins]
            ax.bar(xs, heights, width=width, label=sample)
        ax.set_xticks([b + width * (len(samples) - 1) / 2 for b in bins])
        ax.set_xticklabels(bins)
        ax.legend(fontsize=8)
    else:
        bins = list(range(1, max_n + 1))
        heights = [total.get(n, 0) for n in bins]
        ax.bar(bins, heights, color="#4C72B0")
        ax.set_xticks(bins)

    ax.set_xlabel("NUM_CALLERS")
    ax.set_ylabel("Number of variant sites")
    ax.set_title("Variant sites by number of supporting callers")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
