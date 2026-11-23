#!/usr/bin/env python3
"""
Diagnostic: smooth WhatsHap's --recombination-list into windows and compare
siblings, as a cross-check on hap-ibd.

WhatsHap reports every switch in the transmitted parental haplotype, which on
ONT data includes many short flip-flops (genotype / phase noise). This takes,
per child and parent, the transmitted haplotype (0/1) along each chromosome,
calls each -window_mb window by majority (bp-weighted) and then:

  * counts state changes between consecutive called windows (~ crossovers;
    expect roughly 25-30 paternal and 40-45 maternal per child genome-wide)
  * for every pair of children, the fraction of windows where they carry the
    same paternal / maternal haplotype and IBD0/1/2 (expect ~25/50/25)

Caveat: the 0/1 labels are relative to the parent's phase, so a parental
phase-block switch looks like a crossover in EVERY child at the same spot --
those show up as changes shared by all children (column n_shared_changes).

    recombination_windows.py -recomb FAM.recombination.tsv [-window_mb 1] [-out prefix]
"""
import argparse
import sys
from collections import defaultdict


def log(msg):
    sys.stderr.write(f"[recombination_windows] {msg}\n")


def read_events(path):
    """{(child, chrom): [(pos1, pos2, (f_before, f_after), (m_before, m_after))]}"""
    ev = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            x = line.split()
            if len(x) < 8:
                continue
            ev[(x[0], x[1])].append((int(x[2]), int(x[3]), (x[4], x[5]), (x[6], x[7])))
    for v in ev.values():
        v.sort()
    return ev


def intervals(events, which):
    """[(start, end, state)] for one parent (0 = father, 1 = mother); the
    uncertain stretch between pos1 and pos2 of an event is left out."""
    out = []
    first = events[0]
    cur_state = first[2 + which][0]
    cur_start = 0
    for p1, p2, fs, ms in events:
        before, after = (fs, ms)[which]
        if before == after:
            continue                     # the other parent switched here
        out.append((cur_start, p1, before))
        cur_start, cur_state = p2, after
    out.append((cur_start, None, cur_state))
    return out


def window_states(iv, win, chrom_end):
    """{window index: majority state} (bp-weighted) over the known stretches;
    the open last stretch runs to chrom_end (last event on this chromosome
    in any child)."""
    bp = defaultdict(lambda: defaultdict(int))
    for s, e, st in iv:
        e = e if e is not None else chrom_end
        w = s // win
        while w * win < e:
            lo, hi = max(s, w * win), min(e, (w + 1) * win)
            if hi > lo:
                bp[w][st] += hi - lo
            w += 1
    return {w: max(d, key=d.get) for w, d in bp.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-recomb", required=True)
    p.add_argument("-window_mb", type=float, default=1.0)
    p.add_argument("-out", default="recombination_windows")
    args = p.parse_args()
    win = int(args.window_mb * 1e6)

    ev = read_events(args.recomb)
    children = sorted({c for c, _ in ev})
    chroms = sorted({ch for _, ch in ev})
    log(f"{sum(len(v) for v in ev.values()):,} raw events, children {children}, {len(chroms)} chromosomes")

    # states[child][parent][chrom] = {window: state}
    chrom_end = defaultdict(int)
    for (c, ch), events in ev.items():
        chrom_end[ch] = max(chrom_end[ch], max(e[1] for e in events))
    states = {c: {0: {}, 1: {}} for c in children}
    for (c, ch), events in ev.items():
        for which in (0, 1):
            states[c][which][ch] = window_states(intervals(events, which), win, chrom_end[ch])

    names = ("paternal", "maternal")
    with open(f"{args.out}.windows.tsv", "w") as out:
        out.write("child\tparent\tchrom\twindow_start\tstate\n")
        for c in children:
            for which in (0, 1):
                for ch in chroms:
                    for w, st in sorted(states[c][which].get(ch, {}).items()):
                        out.write(f"{c}\t{names[which]}\t{ch}\t{w * win}\t{st}\n")

    # crossovers per child / parent, and how many sit at the same window in every child
    change_sets = defaultdict(set)
    with open(f"{args.out}.crossovers.tsv", "w") as out:
        out.write("child\tparent\tn_raw_switches\tn_window_changes\n")
        for c in children:
            for which in (0, 1):
                raw = sum(1 for (cc, ch), e in ev.items() if cc == c
                          for r in e if r[2 + which][0] != r[2 + which][1])
                n = 0
                for ch in chroms:
                    ws = sorted(states[c][which].get(ch, {}).items())
                    for (w0, s0), (w1, s1) in zip(ws, ws[1:]):
                        if s0 != s1:
                            n += 1
                            change_sets[(which, ch, w1)].add(c)
                out.write(f"{c}\t{names[which]}\t{raw}\t{n}\n")
                log(f"{c} {names[which]}: {raw:,} raw switches -> {n} changes between {args.window_mb:g} Mb windows")
    for which in (0, 1):
        shared = sum(1 for (wh, _, _), cs in change_sets.items() if wh == which and len(cs) == len(children) > 1)
        log(f"{names[which]}: {shared} window changes seen in ALL children (parental phase switches, not crossovers)")

    # sibling sharing
    with open(f"{args.out}.sib_sharing.tsv", "w") as out:
        out.write("child1\tchild2\tn_windows\tfrac_pat_shared\tfrac_mat_shared\tfrac_IBD0\tfrac_IBD1\tfrac_IBD2\n")
        for i, a in enumerate(children):
            for b in children[i + 1:]:
                n = pat = mat = 0
                ibd = [0, 0, 0]
                for ch in chroms:
                    pa, pb = states[a][0].get(ch, {}), states[b][0].get(ch, {})
                    ma, mb = states[a][1].get(ch, {}), states[b][1].get(ch, {})
                    for w in set(pa) & set(pb) & set(ma) & set(mb):
                        n += 1
                        sp, sm = pa[w] == pb[w], ma[w] == mb[w]
                        pat += sp
                        mat += sm
                        ibd[sp + sm] += 1
                if n:
                    out.write(f"{a}\t{b}\t{n}\t{pat / n:.3f}\t{mat / n:.3f}\t"
                              f"{ibd[0] / n:.3f}\t{ibd[1] / n:.3f}\t{ibd[2] / n:.3f}\n")
                    log(f"{a} vs {b}: {n} windows; paternal shared {pat / n:.1%}, maternal shared {mat / n:.1%}; "
                        f"IBD0/1/2 {ibd[0] / n:.1%} / {ibd[1] / n:.1%} / {ibd[2] / n:.1%}")


if __name__ == "__main__":
    main()
