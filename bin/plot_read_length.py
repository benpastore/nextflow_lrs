#!/usr/bin/env python3
"""
Publication-quality read-length histogram (number of reads vs read length)
with the N50 marked, for one sample's FASTQ -- a more readable alternative
to NanoPlot's length plots.

Log-scaled x axis with log-spaced bins (ONT read lengths span 3-4 orders of
magnitude, so linear bins squash everything into the first few bars), one
series, N50 as a labelled reference line, and a small summary block (reads,
bases, median, mean, N50). Writes PNG (300 dpi) + PDF (vector, editable text)
plus the stats and the exact length -> count table, so the figure can be
restyled or combined across samples later without re-reading the FASTQ.

    plot_read_length.py -fastq S1.fastq.gz -sample S1 -prefix S1.read_length
"""
import argparse
import gzip
import shutil
import subprocess
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullFormatter

# Static print figure: light surface only. Values from the dataviz reference
# palette (series slot 1 + neutral chrome/ink tokens).
SERIES = "#2a78d6"
SURFACE = "#ffffff"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRIDLINE = "#e1e0d9"
AXIS = "#c3c2b7"


def read_lengths(path):
    """Counter of read length -> number of reads. Sequence is every 4th line
    of the FASTQ; decompresses with pigz when available (much faster on
    tens of GB of reads), else gzip."""
    counts = Counter()
    if path.endswith(".gz") and shutil.which("pigz"):
        proc = subprocess.Popen(["pigz", "-dc", path], stdout=subprocess.PIPE)
        fh = proc.stdout
    else:
        proc = None
        fh = gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")
    for i, line in enumerate(fh):
        if i % 4 == 1:
            counts[len(line) - 1] += 1  # minus the trailing newline
    fh.close()
    if proc is not None and proc.wait() != 0:
        sys.exit(f"plot_read_length.py: pigz failed on {path}")
    return counts


def summarize(counts):
    lengths = np.array(sorted(counts), dtype=np.int64)
    n = np.array([counts[l] for l in lengths], dtype=np.int64)
    n_reads = int(n.sum())
    bases = lengths * n
    total = int(bases.sum())

    # N50: length L such that reads >= L hold at least half of all bases
    desc_cum = np.cumsum(bases[::-1])
    n50 = int(lengths[::-1][np.searchsorted(desc_cum, total / 2)])
    median = int(lengths[np.searchsorted(np.cumsum(n), (n_reads + 1) / 2)])
    return {
        "reads": n_reads, "bases": total, "mean": total / n_reads, "median": median,
        "n50": n50, "min": int(lengths[0]), "max": int(lengths[-1]),
    }, lengths, n


def fmt_bp(x, _pos=None):
    if x >= 1e6:
        return f"{x / 1e6:g} Mb"
    if x >= 1e3:
        return f"{x / 1e3:g} kb"
    return f"{x:g} bp"


def fmt_count(y, _pos=None):
    if y >= 1e6:
        return f"{y / 1e6:g}M"
    if y >= 1e3:
        return f"{y / 1e3:g}K"
    return f"{y:g}"


def fmt_bases(b):
    for unit, scale in (("Gb", 1e9), ("Mb", 1e6), ("kb", 1e3)):
        if b >= scale:
            return f"{b / scale:.1f} {unit}"
    return f"{b} bp"


def fmt_len(l):
    return f"{l / 1e3:.1f} kb" if l >= 1e3 else f"{l:.0f} bp"


def plot(lengths, n, stats, sample, prefix, n_bins, min_length, width, height, font_size):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": font_size,
        "pdf.fonttype": 42,  # embed TrueType so text stays editable in Illustrator/Inkscape
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    lo = max(min_length, stats["min"], 1)
    edges = np.logspace(np.log10(lo), np.log10(stats["max"] + 1), n_bins + 1)
    hist, _ = np.histogram(lengths, bins=edges, weights=n)

    fig, ax = plt.subplots(figsize=(width, height))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # thin surface-colored edge = the gap between touching bars
    ax.bar(edges[:-1], hist, width=np.diff(edges), align="edge",
           color=SERIES, edgecolor=SURFACE, linewidth=0.4, zorder=2)

    ax.set_xscale("log")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, hist.max() * 1.12)
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_bp))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))  # clean 1/2/5 ticks, no 1.25K
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_count))

    ax.set_xlabel("Read length", color=INK_PRIMARY)
    ax.set_ylabel("Number of reads", color=INK_PRIMARY)
    if sample:
        ax.set_title(sample, color=INK_PRIMARY, loc="left", fontweight="bold")

    # recessive chrome: hairline y grid only, no top/right spines
    ax.grid(axis="y", color=GRIDLINE, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK_SECONDARY, which="both", length=3, width=0.8)
    ax.tick_params(which="minor", length=1.5)

    # N50 reference line, labelled directly
    ax.axvline(stats["n50"], color=INK_PRIMARY, linewidth=1.0, zorder=3)
    ax.annotate(f"N50 = {fmt_len(stats['n50'])}", xy=(stats["n50"], hist.max() * 1.08),
                xytext=(4, 0), textcoords="offset points", ha="left", va="center",
                color=INK_PRIMARY, fontsize=font_size)

    summary = (f"Reads  {stats['reads']:,}\n"
               f"Bases  {fmt_bases(stats['bases'])}\n"
               f"Median  {fmt_len(stats['median'])}\n"
               f"Mean  {fmt_len(stats['mean'])}")
    # put the summary on whichever side the N50 label isn't
    right_side = stats["n50"] < np.sqrt(edges[0] * edges[-1])
    ax.text(0.98 if right_side else 0.02, 0.80, summary, transform=ax.transAxes,
            ha="right" if right_side else "left", va="top",
            color=INK_SECONDARY, fontsize=font_size * 0.9, linespacing=1.5)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{prefix}.{ext}", dpi=300, facecolor=SURFACE)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-fastq", required=True, help="Reads (.fastq or .fastq.gz)")
    p.add_argument("-sample", default="", help="Sample ID, used as the plot title")
    p.add_argument("-prefix", required=True, help="Output prefix (.png/.pdf/.stats.tsv/.counts.tsv.gz)")
    p.add_argument("-bins", type=int, default=80, help="Number of log-spaced bins (default 80)")
    p.add_argument("-min_length", type=int, default=100,
                   help="Left edge of the x axis, bp (shorter reads still count in the stats; default 100)")
    p.add_argument("-width", type=float, default=5.0, help="Figure width, inches (default 5)")
    p.add_argument("-height", type=float, default=3.4, help="Figure height, inches (default 3.4)")
    p.add_argument("-font_size", type=float, default=10, help="Base font size, pt (default 10)")
    args = p.parse_args()

    counts = read_lengths(args.fastq)
    if not counts:
        sys.exit(f"plot_read_length.py: no reads in {args.fastq}")
    stats, lengths, n = summarize(counts)

    with open(f"{args.prefix}.stats.tsv", "w") as out:
        out.write("sample\treads\tbases\tmean\tmedian\tN50\tmin\tmax\n")
        out.write(f"{args.sample}\t{stats['reads']}\t{stats['bases']}\t{stats['mean']:.1f}\t"
                  f"{stats['median']}\t{stats['n50']}\t{stats['min']}\t{stats['max']}\n")
    with gzip.open(f"{args.prefix}.counts.tsv.gz", "wt") as out:
        out.write("length\tn_reads\n")
        for l, c in zip(lengths, n):
            out.write(f"{l}\t{c}\n")

    plot(lengths, n, stats, args.sample, args.prefix, args.bins, args.min_length,
         args.width, args.height, args.font_size)
    sys.stderr.write(f"[plot_read_length] {args.sample}: {stats['reads']:,} reads, "
                     f"{fmt_bases(stats['bases'])}, N50 {fmt_len(stats['n50'])} -> {args.prefix}.png/.pdf\n")


if __name__ == "__main__":
    main()
