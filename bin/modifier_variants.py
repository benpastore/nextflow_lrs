#!/usr/bin/env python3
"""
Modifier-variant filtering for the whole cohort in one run: every mild
child against its severe sibling(s), from their annotated master VCFs (consolidate_variants.py
-> snpEff -> annotate_af.py: GT:PS per caller, INFO/ORIGIN, ANN,
gnomAD_AF_NFE, KG_AF_EUR). Plain genotype rules only -- no statistical model,
no allele-frequency filter (AF is carried through as annotation).

A variant is GENIC when snpEff puts it in a protein-coding gene with one of
GENIC_TERMS (CDS-affecting, splice, exon, intron, UTR); otherwise it is
REGULATORY and is linked to the nearest protein-coding gene on each side
(from -gtf, strand ignored, with distances in bp).

Families come from -family_json (children grouped by father + mother):
  discordant   >= 1 mild and >= 1 severe child -> the rules are applied to
               each mild child against its severe sib(s)
  concordant   severe children only -> their severe children form the
               CONCORDANT SEVERE group, indexed once and used only as extra
               exclusion controls (below); nothing about them is written out
  (families with mild children only are skipped; a severe child without
  parents in family.json counts as its own concordant family). With no
  concordant family the concordant step is skipped entirely and the rules
  are the plain sib-vs-sib ones.

Matching a mild variant in a severe sib: small variants by exact
chrom/pos/ref/alt; SVs by chrom + SVTYPE with POS within -sv_match_dist bp.
"In severe" = present (het or hom) in ANY severe sib of the family OR in
ANY concordant severe child: every exclusion below applies to both, so a
variant a concordant severe child carries the same way can't make a gene
pass. Outputs only describe the family's own severe sib(s).

Gene-level rules (per gene, for the mild child):
  homozygous    the mild child is homozygous for a genic variant and no
                severe sib / concordant severe child is homozygous for that
                same variant (severe het or absent is fine)
  compound_het  the mild child has >= 1 het on the paternal copy and >= 1 het
                on the maternal copy of the gene (INFO/ORIGIN), and at least
                one of those hets is absent from every severe sib and every
                concordant severe child. Hets with no ORIGIN fall back to
                phase: two such hets in the same phase block (same caller,
                same PS) on opposite haplotypes, at least one absent likewise.
                SVs count like any variant.
  regulatory    a non-genic mild variant that no severe sib / concordant
                severe child carries with the same zygosity; both flanking
                protein-coding genes pass.

Every output variant gets, per mild x severe pair, the hap-ibd sharing state
at its position (cohort.sib_ibd_segments.tsv): IBD (IBD0/1/2), IBD_PAT and
IBD_MAT (shared / not_shared); '.' where the pair wasn't scanned (e.g. chrX).
With one severe sib the INFO keys are IBD / IBD_PAT / IBD_MAT; with several,
each is suffixed with the severe sib's ID (IBD_<severe>, ...).

The rules only SELECT genes. A gene passing doesn't mean the severe sib's
copy is intact -- it may carry its own variants there, and what differs
between the sibs is what those variants do. So for every passing gene the
full variant content of BOTH sides is recorded (<mild>.gene_variants.tsv.gz):
every variant the mild child or any severe sib carries that is genic in the
gene, plus the mild child's passing regulatory sites linked to it. Each is
labelled by who carries it:
  mild_only     mild child only
  severe_only   severe sib(s) only
  shared        mild child and every severe sib, same zygosity
  shared_diff   mild child and >= 1 severe sib, zygosity differs somewhere
                (incl. a severe sib without it when there are several)
gene_recurrence.py turns this table into the AlphaGenome / AlphaMissense
input, so both sibs' variants in the gene get scored and can be compared.
Severe-only regulatory variants near a gene are not included (the
regulatory rule is defined from the mild child's sites).

Outputs (in <family>/ under the working directory, family = <father>_<mother>;
<mild> = the mild child's ID):
  <mild>.passing_genes.tsv          one row per passing gene
  <mild>.modifier_variants.tsv.gz   one row per (gene, supporting variant)
  <mild>.gene_variants.tsv.gz       one row per (passing gene, variant in mild
                                    OR severe), with carriers + zygosities
  <mild>.mild.vcf                   the mild child's supporting variants
  <severe>.severe.vcf               per severe sib: its variants in every
                                    passing gene, plus its calls at the
                                    passing regulatory sites
  <mild>.modifier_summary.tsv       counts at each step (incl. how many
                                    variants the concordant group excluded)

    modifier_variants.py -family_json family.json \\
        -sample_vcf S1=S1.snpeff.af.vcf.gz -sample_vcf S2=S2.snpeff.af.vcf.gz ... \\
        -ibd_segments cohort.sib_ibd_segments.tsv -gtf genes.gtf
"""
import argparse
import bisect
import gzip
import json
import os
import re
import sys
from collections import defaultdict

from _modifier_common import CALLER_PRIORITY, classify_gt, get_ann, get_info_field, opener
from _cnv_repeat_common import load_gtf, present

# snpEff terms that place a variant inside a gene (protein-coding transcripts).
GENIC_TERMS = {
    # CDS-affecting
    "missense_variant", "synonymous_variant", "stop_gained", "stop_lost", "start_lost",
    "start_retained_variant", "stop_retained_variant", "initiator_codon_variant",
    "frameshift_variant", "inframe_insertion", "inframe_deletion",
    "disruptive_inframe_insertion", "disruptive_inframe_deletion",
    "conservative_inframe_insertion", "conservative_inframe_deletion",
    "coding_sequence_variant", "protein_altering_variant", "rare_amino_acid_variant",
    # splice
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "splice_donor_5th_base_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    # exon / intron
    "exon_variant", "exon_loss_variant", "intron_variant",
    # UTR
    "5_prime_UTR_variant", "3_prime_UTR_variant",
    "5_prime_UTR_premature_start_codon_gain_variant",
    "5_prime_UTR_truncation", "3_prime_UTR_truncation",
    # structural variants overlapping a gene
    "transcript_ablation", "feature_ablation", "gene_fusion", "bidirectional_gene_fusion",
    "transcript_amplification", "feature_elongation", "feature_truncation",
}

ZYG_RANK = {"absent": 0, "het": 1, "homo": 2}
REASONS = ("homozygous", "compound_het", "regulatory")


def log(msg):
    sys.stderr.write(f"[modifier_variants] {msg}\n")


def safe_key(s):
    return re.sub(r"[^0-9A-Za-z_.]", "_", s)


# ------------------------------------------------------------------ VCF parsing

class Rec(object):
    __slots__ = ("idx", "chrom", "pos", "ref", "alt", "vid", "svtype", "end", "zyg", "gt", "ps",
                 "caller", "origin", "genic", "gnomad", "kg")


def read_header(path):
    header, callers = [], []
    with opener(path) as f:
        for line in f:
            if line.startswith("##"):
                header.append(line)
            elif line.startswith("#CHROM"):
                callers = line.rstrip("\n").split("\t")[9:]
                break
    return header, callers


def genic_genes(info, coding_genes):
    """Protein-coding genes this variant is genic in, from snpEff ANN."""
    ann = get_ann(info)
    genes = set()
    if not ann:
        return genes
    for a in ann:
        if len(a) < 8 or not a[3]:
            continue
        if not any(t in GENIC_TERMS for t in a[1].split("&")):
            continue
        # snpEff built from a UCSC GTF can leave the biotype empty -- fall
        # back to "the GTF has a CDS for this gene"
        if a[7] == "protein_coding" or (not a[7] and a[3] in coding_genes):
            genes.update(g for g in a[3].split("&") if g)
    return genes


def iter_records(path, coding_genes, with_genic=True):
    """Yields (Rec, raw_line) for every called (het/hom) record. Rec.idx is the
    record's ordinal among called records -- stable across passes over a file."""
    _, callers = read_header(path)
    order = sorted(range(len(callers)),
                   key=lambda i: CALLER_PRIORITY.index(callers[i]) if callers[i] in CALLER_PRIORITY else 99)
    n = 0
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 10:
                continue
            fmt = x[8].split(":")
            zyg = gt = ps = caller = None
            for i in order:
                d = dict(zip(fmt, x[9 + i].split(":")))
                z = classify_gt(d.get("GT", ""))
                if z:
                    zyg, gt, ps, caller = z, d.get("GT", ""), d.get("PS", "."), callers[i]
                    break
            if zyg is None:
                continue
            info = x[7]
            r = Rec()
            r.idx = n
            n += 1
            r.chrom, r.pos, r.vid, r.ref, r.alt = x[0], int(x[1]), x[2], x[3], x[4]
            r.svtype = get_info_field(info, "SVTYPE")
            span_end = r.pos + max(len(r.ref), 1) - 1
            try:
                r.end = int(get_info_field(info, "END") or span_end)
            except ValueError:
                r.end = span_end
            r.zyg, r.gt, r.ps, r.caller = zyg, gt, ps, caller
            r.origin = get_info_field(info, "ORIGIN") or "unknown"
            r.gnomad = get_info_field(info, "gnomAD_AF_NFE") or "."
            r.kg = get_info_field(info, "KG_AF_EUR") or "."
            r.genic = genic_genes(info, coding_genes) if with_genic else set()
            yield r, line


def small_key(r):
    return (r.chrom, r.pos, r.ref, r.alt)


class SiteIndex(object):
    """Value lookup for a set of variants: exact key for small variants,
    chrom + SVTYPE + POS within sv_dist for SVs. merge(old, new) resolves
    two entries landing on the same small key / an SV window."""

    def __init__(self, items, sv_dist, merge=None):
        self.sv_dist = sv_dist
        self.merge = merge or (lambda old, new: old)
        self.small = {}
        sv = defaultdict(list)
        for r, val in items:
            if r.svtype:
                sv[(r.chrom, r.svtype)].append((r.pos, val))
            else:
                k = small_key(r)
                self.small[k] = self.merge(self.small[k], val) if k in self.small else val
        self.sv = {}
        for k, lst in sv.items():
            lst.sort(key=lambda t: t[0])
            self.sv[k] = ([p for p, _ in lst], [v for _, v in lst])

    def find(self, r, default=None):
        if not r.svtype:
            return self.small.get(small_key(r), default)
        entry = self.sv.get((r.chrom, r.svtype))
        if not entry:
            return default
        pos, vals = entry
        lo = bisect.bisect_left(pos, r.pos - self.sv_dist)
        hi = bisect.bisect_right(pos, r.pos + self.sv_dist)
        if lo >= hi:
            return default
        best = vals[lo]
        for v in vals[lo + 1:hi]:
            best = self.merge(best, v)
        return best


class SeenSet(object):
    """Growing set of variants with the same matching as SiteIndex (exact
    small-variant key; SV = chrom + SVTYPE + POS within sv_dist)."""

    def __init__(self, sv_dist):
        self.sv_dist = sv_dist
        self.small = set()
        self.sv = defaultdict(list)

    def add(self, r):
        if r.svtype:
            self.sv[(r.chrom, r.svtype)].append(r.pos)
        else:
            self.small.add(small_key(r))

    def __contains__(self, r):
        if not r.svtype:
            return small_key(r) in self.small
        return any(abs(p - r.pos) <= self.sv_dist for p in self.sv.get((r.chrom, r.svtype), ()))


def carrier_class(zyg, mild, severe):
    """mild_only / severe_only / shared / shared_diff (see module docstring)."""
    mz = zyg[mild]
    sz = [zyg[s] for s in severe]
    if mz == "absent":
        return "severe_only"
    if all(z == "absent" for z in sz):
        return "mild_only"
    return "shared" if all(z == mz for z in sz) else "shared_diff"


CARRIER_CLASSES = ("mild_only", "severe_only", "shared", "shared_diff")


def stronger(a, b):
    return a if ZYG_RANK[a] >= ZYG_RANK[b] else b


def index_severe(sample, path, coding_genes, sv_dist):
    """{variant -> 'het'/'homo'} for one severe sib (hom wins on collisions)."""
    items = [(r, r.zyg) for r, _ in iter_records(path, coding_genes, with_genic=False)]
    log(f"{sample}: {len(items):,} called records indexed")
    return SiteIndex(items, sv_dist, merge=stronger)


class ConcordantIndex(object):
    """Pooled concordant-severe calls: variant -> set of zygosities seen in
    ANY concordant severe child (empty set = nobody carries it). Built one
    sample at a time so only the pooled index is held in memory."""

    def __init__(self, sv_dist):
        self.sv_dist = sv_dist
        self.small = defaultdict(set)
        self.sv = defaultdict(list)     # (chrom, svtype) -> [(pos, zyg)], sorted by find()
        self._sorted = True

    def add_sample(self, sample, path, coding_genes):
        n = 0
        for r, _ in iter_records(path, coding_genes, with_genic=False):
            n += 1
            if r.svtype:
                self.sv[(r.chrom, r.svtype)].append((r.pos, r.zyg))
                self._sorted = False
            else:
                self.small[small_key(r)].add(r.zyg)
        log(f"concordant severe {sample}: {n:,} called records pooled")

    def find(self, r):
        if not r.svtype:
            return self.small.get(small_key(r), frozenset())
        if not self._sorted:
            for k in self.sv:
                self.sv[k].sort()
            self.pos = {k: [p for p, _ in v] for k, v in self.sv.items()}
            self._sorted = True
        lst = self.sv.get((r.chrom, r.svtype))
        if not lst:
            return frozenset()
        pos = self.pos[(r.chrom, r.svtype)]
        lo = bisect.bisect_left(pos, r.pos - self.sv_dist)
        hi = bisect.bisect_right(pos, r.pos + self.sv_dist)
        return frozenset(z for _, z in lst[lo:hi])


def group_families(fam, samples):
    """{family: {'mild': [...], 'severe': [...]}} for the children given,
    family = <father>_<mother> (a child without both parents is its own family)."""
    groups = defaultdict(lambda: {"mild": [], "severe": []})
    for s in samples:
        info = fam.get(s, {})
        ph = info.get("phenotype")
        if ph not in ("mild", "severe"):
            continue
        fid = f"{info['father']}_{info['mother']}" if info.get("father") and info.get("mother") else s
        groups[fid][ph].append(s)
    for g in groups.values():
        g["mild"].sort()
        g["severe"].sort()
    return groups


# ------------------------------------------------------------------ genes / IBD

class FlankIndex(object):
    """Nearest protein-coding gene on each side of an interval (0-based half-open)."""

    def __init__(self, loci):
        by_start, by_end = defaultdict(list), defaultdict(list)
        for g in loci:
            by_start[g.chrom].append((g.start, g.gene))
            by_end[g.chrom].append((g.end, g.gene))
        self.starts, self.ends = {}, {}
        for c in by_start:
            by_start[c].sort()
            by_end[c].sort()
            self.starts[c] = ([s for s, _ in by_start[c]], [g for _, g in by_start[c]])
            self.ends[c] = ([e for e, _ in by_end[c]], [g for _, g in by_end[c]])

    def flanks(self, chrom, start, end):
        """((left_gene, dist_bp) or None, (right_gene, dist_bp) or None)."""
        left = right = None
        if chrom in self.ends:
            ends, genes = self.ends[chrom]
            i = bisect.bisect_right(ends, start) - 1
            if i >= 0:
                left = (genes[i], start - ends[i])
        if chrom in self.starts:
            starts, genes = self.starts[chrom]
            j = bisect.bisect_left(starts, end)
            if j < len(starts):
                right = (genes[j], starts[j] - end)
        return left, right


def load_ibd(path, mild, severe):
    """{severe: {chrom: (starts, ends, [(ibd, pat, mat)])}} for this mild child."""
    if not present(path):
        log("no IBD segments file -- IBD fields will be '.'")
        return {s: {} for s in severe}
    raw = {s: defaultdict(list) for s in severe}
    with opener(path) as f:
        col = {h: i for i, h in enumerate(f.readline().rstrip("\n").split("\t"))}
        for line in f:
            x = line.rstrip("\n").split("\t")
            if x[col["mild"]] != mild or x[col["severe"]] not in raw:
                continue
            raw[x[col["severe"]]][x[col["chrom"]]].append(
                (int(x[col["start"]]), int(x[col["end"]]),
                 (x[col["ibd"]], x[col["paternal"]], x[col["maternal"]])))
    idx = {}
    for s, by_chrom in raw.items():
        idx[s] = {}
        for c, segs in by_chrom.items():
            segs.sort(key=lambda t: t[0])
            idx[s][c] = ([a for a, _, _ in segs], [b for _, b, _ in segs], [v for _, _, v in segs])
        log(f"IBD {mild} vs {s}: {sum(len(v[0]) for v in idx[s].values()):,} segments")
    return idx


def ibd_at(idx, chrom, pos):
    entry = idx.get(chrom)
    if entry:
        starts, ends, vals = entry
        i = bisect.bisect_right(starts, pos) - 1
        if i >= 0 and pos <= ends[i]:
            return vals[i]
    return (".", ".", ".")


# ------------------------------------------------------------------ rules

def phase_side(r):
    """0/1 = which haplotype of the phase block carries ALT, or None if unphased."""
    if "|" not in r.gt or r.ps in ("", "."):
        return None
    return 0 if r.gt.split("|")[0] != "0" else 1


def compound_het(hets, novel):
    """hets: [Rec] of one gene. Returns {rec.idx: haplotype label} of the
    supporting hets (empty = the gene doesn't pass)."""
    support = {}
    pat = [r for r in hets if r.origin == "paternal"]
    mat = [r for r in hets if r.origin == "maternal"]
    if pat and mat and any(novel[r.idx] for r in pat + mat):
        support.update((r.idx, "paternal") for r in pat)
        support.update((r.idx, "maternal") for r in mat)
    # fallback for hets without ORIGIN: same phase block, opposite haplotypes
    blocks = defaultdict(lambda: ([], []))
    for r in hets:
        if r.origin not in ("paternal", "maternal"):
            side = phase_side(r)
            if side is not None:
                blocks[(r.caller, r.ps)][side].append(r)
    for (caller, ps), (s0, s1) in blocks.items():
        if s0 and s1 and any(novel[r.idx] for r in s0 + s1):
            support.update((r.idx, f"{caller}_PS{ps}_hap1") for r in s0)
            support.update((r.idx, f"{caller}_PS{ps}_hap2") for r in s1)
    return support


def fmt_flank(side):
    return f"{side[0]}:{side[1]}" if side else "."


# ------------------------------------------------------------------ output

def vcf_header(header, extra_info, last_line):
    lines = [h for h in header if not h.startswith("##fileformat")]
    out = ["##fileformat=VCFv4.2\n"] + lines
    out += [f'##INFO=<ID={k},Number={n},Type=String,Description="{d}">\n' for k, n, d in extra_info]
    out.append(last_line)
    return "".join(out)


def ibd_info_defs(severe, ibd_keys):
    defs = []
    for s in severe:
        k_ibd, k_pat, k_mat = ibd_keys(s)
        defs += [(k_ibd, "1", f"hap-ibd state vs severe sib {s} at this position: IBD0/IBD1/IBD2 (. = not scanned)"),
                 (k_pat, "1", f"Paternal copy shared with severe sib {s}: shared/not_shared"),
                 (k_mat, "1", f"Maternal copy shared with severe sib {s}: shared/not_shared")]
    return defs


def ibd_info(r, severe, ibd_idx, ibd_keys):
    fields = []
    for s in severe:
        ibd, pat, mat = ibd_at(ibd_idx[s], r.chrom, r.pos)
        k_ibd, k_pat, k_mat = ibd_keys(s)
        fields += [f"{k_ibd}={ibd}", f"{k_pat}={pat}", f"{k_mat}={mat}"]
    return fields


def add_info(line, fields):
    x = line.rstrip("\n").split("\t")
    x[7] = ";".join(([x[7]] if x[7] not in ("", ".") else []) + fields)
    return "\t".join(x) + "\n"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-family_json", required=True)
    p.add_argument("-sample_vcf", action="append", required=True,
                   help="SAMPLE=annotated master VCF, once per mild / severe child of the cohort")
    p.add_argument("-ibd_segments", default="", help="hapibd_segments.py cohort.sib_ibd_segments.tsv (or NO_FILE)")
    p.add_argument("-gtf", required=True, help="gene model (protein-coding = transcripts with a CDS)")
    p.add_argument("-sv_match_dist", type=int, default=500)
    args = p.parse_args()

    with open(args.family_json) as f:
        fam = json.load(f)
    # absolute paths: each family's outputs are written from inside <family>/
    vcfs = {s: os.path.abspath(v) for s, v in (x.split("=", 1) for x in args.sample_vcf)}
    if present(args.ibd_segments):
        args.ibd_segments = os.path.abspath(args.ibd_segments)

    groups = group_families(fam, vcfs)
    discordant = {f: g for f, g in groups.items() if g["mild"] and g["severe"]}
    concordant = sorted(s for g in groups.values() if not g["mild"] for s in g["severe"])
    mild_only = sorted(f for f, g in groups.items() if g["mild"] and not g["severe"])
    log(f"{len(discordant)} discordant famil{'y' if len(discordant) == 1 else 'ies'}: "
        + "; ".join(f"{f} (mild {g['mild']}, severe {g['severe']})" for f, g in sorted(discordant.items())))
    if concordant:
        log(f"{len(concordant)} concordant severe child(ren) as extra exclusion controls: {concordant}")
    else:
        log("no concordant (severe-only) families in family.json -- concordant exclusion skipped")
    if mild_only:
        log(f"skipped, no severe sib: {mild_only}")
    if not discordant:
        sys.exit("[modifier_variants] no family with >= 1 mild and >= 1 severe child")

    _, loci = load_gtf(args.gtf)
    coding = [g for g in loci if g.coding]
    coding_genes = {g.gene for g in coding}
    flanks = FlankIndex(coding)
    conc_idx = None                     # None = no concordant group: sib-only rules
    if concordant:
        conc_idx = ConcordantIndex(args.sv_match_dist)
        for s in concordant:
            conc_idx.add_sample(s, vcfs[s], coding_genes)

    top = os.getcwd()
    for family in sorted(discordant):
        severe = discordant[family]["severe"]
        sev_idx = {s: index_severe(s, vcfs[s], coding_genes, args.sv_match_dist) for s in severe}
        multi = len(severe) > 1

        def ibd_keys(s, multi=multi):
            sfx = f"_{safe_key(s)}" if multi else ""
            return f"IBD{sfx}", f"IBD_PAT{sfx}", f"IBD_MAT{sfx}"

        os.makedirs(family, exist_ok=True)
        os.chdir(family)
        try:
            # one mild child per family is expected; more are each run on their own
            for m in discordant[family]["mild"]:
                run_mild(args, family, m, vcfs, severe, sev_idx, conc_idx, coding_genes, flanks, ibd_keys)
        finally:
            os.chdir(top)


def run_mild(args, family, m, vcfs, severe, sev_idx, conc_idx, coding_genes, flanks, ibd_keys):
    counts = defaultdict(int)
    ibd_idx = load_ibd(args.ibd_segments, m, severe)

    # ---- pass 1 over the mild VCF: apply the rules
    kept = {}                         # idx -> Rec (anything that may be written)
    sev_status = {}                   # idx -> {severe: zygosity}
    in_conc = {}                      # idx -> carried by a concordant severe child (exclusion only)
    hets_by_gene = defaultdict(list)
    var_genes = defaultdict(lambda: defaultdict(set))   # idx -> gene -> {reason}
    var_flank = {}                    # idx -> (left, right)
    hom_genes = set()
    for r, _ in iter_records(vcfs[m], coding_genes):
        counts["mild_called_records"] += 1
        st = {s: sev_idx[s].find(r, "absent") for s in severe}
        conc = conc_idx.find(r) if conc_idx is not None else frozenset()
        if r.genic:
            counts["genic_records"] += 1
            if r.zyg == "homo":
                counts["genic_hom"] += 1
                if any(z == "homo" for z in st.values()):
                    continue
                if "homo" in conc:
                    counts["genic_hom_excluded_by_concordant"] += 1
                    continue
                counts["genic_hom_passing"] += 1
                for g in r.genic:
                    var_genes[r.idx][g].add("homozygous")
                    hom_genes.add(g)
            else:
                counts["genic_het"] += 1
                for g in r.genic:
                    hets_by_gene[g].append(r)
        else:
            counts["regulatory_records"] += 1
            if any(z == r.zyg for z in st.values()):
                continue
            if r.zyg in conc:
                counts["regulatory_excluded_by_concordant"] += 1
                continue
            left, right = flanks.flanks(r.chrom, r.pos - 1, r.end)
            if left is None and right is None:
                counts["regulatory_no_flanking_gene"] += 1
                continue
            counts["regulatory_passing"] += 1
            var_flank[r.idx] = (left, right)
            for side in (left, right):
                if side:
                    var_genes[r.idx][side[0]].add("regulatory")
        kept[r.idx] = r
        sev_status[r.idx] = st
        in_conc[r.idx] = bool(conc)
    log(f"{m}: {counts['mild_called_records']:,} called records "
        f"({counts['genic_records']:,} genic, {counts['regulatory_records']:,} regulatory)")

    # absent from every severe sib AND every concordant severe child
    novel = {idx: all(z == "absent" for z in st.values()) and not in_conc[idx] for idx, st in sev_status.items()}
    counts["hets_absent_in_sibs_but_in_concordant"] = sum(
        1 for idx, st in sev_status.items() if in_conc[idx] and all(z == "absent" for z in st.values()))
    var_hap = {}
    chet_genes = set()
    for g, hets in hets_by_gene.items():
        support = compound_het(hets, novel)
        if support:
            chet_genes.add(g)
            for idx, hap in support.items():
                var_genes[idx][g].add("compound_het")
                var_hap[idx] = hap
    # het records not supporting any compound het are no longer needed
    for idx in [i for i in kept if i not in var_genes]:
        del kept[idx]

    # ---- passing genes
    gene_rows = defaultdict(lambda: defaultdict(int))
    for idx, genes in var_genes.items():
        for g, reasons in genes.items():
            for reason in reasons:
                gene_rows[g][reason] += 1
            if novel[idx]:
                gene_rows[g]["absent_in_severe"] += 1
    sev_str = ",".join(severe)
    with open(f"{m}.passing_genes.tsv", "w") as out:
        out.write("family\tmild\tsevere\tgene\treasons\tn_homozygous_variants\tn_compound_het_variants\t"
                  "n_regulatory_variants\tn_variants_absent_in_severe\n")
        for g in sorted(gene_rows):
            row = gene_rows[g]
            reasons = ",".join(k for k in REASONS if row[k])
            out.write(f"{family}\t{m}\t{sev_str}\t{g}\t{reasons}\t{row['homozygous']}\t"
                      f"{row['compound_het']}\t{row['regulatory']}\t{row['absent_in_severe']}\n")

    # every variant either side carries in a passing gene (gene_variants.tsv.gz):
    # (Rec, {sample: zygosity}, {gene}, {gene: {supporting reason}})
    passing = set(gene_rows)
    gene_vars = []
    seen = SeenSet(args.sv_match_dist)

    # ---- pass 2 over the mild VCF: variant table + mild VCF + mild side of gene_vars
    header, callers = read_header(vcfs[m])
    extra = [("MOD_GENE", ".", "Modifier gene(s) this variant supports"),
             ("MOD_REASON", ".", "gene:reason pairs (homozygous / compound_het / regulatory)"),
             ("MOD_HAP", "1", "Compound-het copy: paternal / maternal (INFO/ORIGIN) or <caller>_PS<block>_hap1/2"),
             ("MOD_FLANK_LEFT", "1", "Regulatory: nearest protein-coding gene on the left, gene:distance_bp"),
             ("MOD_FLANK_RIGHT", "1", "Regulatory: nearest protein-coding gene on the right, gene:distance_bp"),
             ("MOD_SEVERE_GT", ".", "Severe sib genotype at this variant, sample:homo/het/absent")]
    extra += ibd_info_defs(severe, ibd_keys)
    chrom_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(callers) + "\n"
    sev_cols = "\t".join(f"severe_{s}" for s in severe)
    ibd_cols = "\t".join(f"IBD_{s}\tIBD_PAT_{s}\tIBD_MAT_{s}" for s in severe)
    with gzip.open(f"{m}.modifier_variants.tsv.gz", "wt", compresslevel=4) as tsv, \
            open(f"{m}.mild.vcf", "w") as vcf:
        tsv.write("family\tmild\tgene\treason\tchrom\tpos\tid\tref\talt\tsvtype\tend\tzygosity\torigin\t"
                  f"haplotype\tgnomAD_AF_NFE\tKG_AF_EUR\tflank_left\tflank_right\t{sev_cols}\t{ibd_cols}\n")
        vcf.write(vcf_header(header, extra, chrom_line))
        for r, line in iter_records(vcfs[m], coding_genes):
            in_genes = r.genic & passing
            if r.idx in var_flank:
                in_genes |= {x[0] for x in var_flank[r.idx] if x}
            if in_genes:
                zyg = {m: r.zyg}
                zyg.update({s: sev_idx[s].find(r, "absent") for s in severe})
                gene_vars.append((r, zyg, in_genes, var_genes.get(r.idx, {})))
                seen.add(r)
            if r.idx not in var_genes:
                continue
            genes = var_genes[r.idx]
            st = sev_status[r.idx]
            left, right = var_flank.get(r.idx, (None, None))
            hap = var_hap.get(r.idx, ".")
            ibd = [ibd_at(ibd_idx[s], r.chrom, r.pos) for s in severe]
            fields = [f"MOD_GENE={','.join(sorted(genes))}",
                      "MOD_REASON=" + ",".join(f"{g}:{x}" for g in sorted(genes) for x in sorted(genes[g])),
                      "MOD_SEVERE_GT=" + ",".join(f"{s}:{st[s]}" for s in severe)]
            if hap != ".":
                fields.append(f"MOD_HAP={hap}")
            if r.idx in var_flank:
                fields += [f"MOD_FLANK_LEFT={fmt_flank(left)}", f"MOD_FLANK_RIGHT={fmt_flank(right)}"]
            fields += ibd_info(r, severe, ibd_idx, ibd_keys)
            vcf.write(add_info(line, fields))
            counts["mild_vcf_records"] += 1
            common = [r.chrom, str(r.pos), r.vid, r.ref, r.alt, r.svtype or ".", str(r.end), r.zyg, r.origin, hap,
                      r.gnomad, r.kg, fmt_flank(left), fmt_flank(right)] + [st[s] for s in severe] + \
                     [v for t in ibd for v in t]
            for g in sorted(genes):
                for reason in sorted(genes[g]):
                    tsv.write("\t".join([family, m, g, reason] + common) + "\n")

    # ---- severe VCFs (+ severe-only side of gene_vars)
    mild_genes = passing
    reg_sites = SiteIndex([(kept[i], i) for i in var_flank], args.sv_match_dist)
    mild_zyg = index_severe(m, vcfs[m], coding_genes, args.sv_match_dist)   # all of the mild child's calls
    for s in severe:
        sev_extra = [("MOD_GENE", ".", f"Modifier gene(s) of mild sib {m} this variant falls in / flanks"),
                     ("MOD_REGULATORY_SITE", "0", f"At a passing regulatory site of mild sib {m}"),
                     ("MOD_MILD_GT", "1", f"Mild sib {m} genotype at this variant: homo/het/absent")]
        sev_extra += ibd_info_defs([s], ibd_keys)
        sev_header, sev_callers = read_header(vcfs[s])
        sev_chrom = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sev_callers) + "\n"
        with open(f"{s}.severe.vcf", "w") as vcf:
            vcf.write(vcf_header(sev_header, sev_extra, sev_chrom))
            for r, line in iter_records(vcfs[s], coding_genes):
                genes = sorted(r.genic & mild_genes)
                if genes and mild_zyg.find(r, "absent") == "absent" and r not in seen:
                    zyg = {m: "absent"}
                    zyg.update({s2: (r.zyg if s2 == s else sev_idx[s2].find(r, "absent")) for s2 in severe})
                    gene_vars.append((r, zyg, set(genes), {}))
                    seen.add(r)
                reg = reg_sites.find(r)
                if not genes and reg is None:
                    continue
                fields = []
                if reg is not None:
                    left, right = var_flank[reg]
                    genes = sorted(set(genes) | {x[0] for x in (left, right) if x})
                    fields.append("MOD_REGULATORY_SITE")
                fields = [f"MOD_GENE={','.join(genes)}"] + fields + \
                         [f"MOD_MILD_GT={mild_zyg.find(r, 'absent')}"] + ibd_info(r, [s], ibd_idx, ibd_keys)
                vcf.write(add_info(line, fields))
                counts[f"severe_vcf_records_{s}"] += 1

    # ---- gene_variants.tsv.gz: one row per (passing gene, variant in mild or severe)
    gene_reasons = {g: ",".join(k for k in REASONS if gene_rows[g][k]) for g in gene_rows}
    gv_rows = []
    for r, zyg, in_genes, support in gene_vars:
        cls = carrier_class(zyg, m, severe)
        counts[f"gene_variants_{cls}"] += 1
        ibd = [v for s in severe for v in ibd_at(ibd_idx[s], r.chrom, r.pos)]
        for g in sorted(in_genes):
            gv_rows.append((r.chrom, r.pos, g, [
                family, m, g, gene_reasons.get(g, ""), r.chrom, str(r.pos), r.vid, r.ref, r.alt,
                r.svtype or ".", str(r.end), cls, zyg[m]] + [zyg[s] for s in severe] + [
                r.origin, ",".join(sorted(support.get(g, ()))) or ".", r.gnomad, r.kg] + ibd))
    gv_rows.sort(key=lambda t: (t[0], t[1], t[2]))
    with gzip.open(f"{m}.gene_variants.tsv.gz", "wt", compresslevel=4) as out:
        out.write("family\tmild\tgene\tgene_reasons\tchrom\tpos\tid\tref\talt\tsvtype\tend\tcarriers\t"
                  f"mild_zyg_{m}\t" + "\t".join(f"severe_zyg_{s}" for s in severe) + "\torigin\tmild_supporting\t"
                  f"gnomAD_AF_NFE\tKG_AF_EUR\t{ibd_cols}\n")
        for *_, row in gv_rows:
            out.write("\t".join(row) + "\n")
    counts["gene_variant_rows"] = len(gv_rows)

    # ---- summary
    counts["genes_homozygous"] = len(hom_genes)
    counts["genes_compound_het"] = len(chet_genes)
    counts["genes_regulatory"] = sum(1 for g in gene_rows if gene_rows[g]["regulatory"])
    counts["passing_genes"] = len(gene_rows)
    with open(f"{m}.modifier_summary.tsv", "w") as out:
        out.write("family\tmild\tsevere\tstep\tcount\n")
        for k in ("mild_called_records", "genic_records", "genic_hom", "genic_hom_excluded_by_concordant",
                  "genic_hom_passing", "genic_het", "hets_absent_in_sibs_but_in_concordant",
                  "regulatory_records", "regulatory_excluded_by_concordant", "regulatory_passing",
                  "regulatory_no_flanking_gene",
                  "genes_homozygous", "genes_compound_het", "genes_regulatory", "passing_genes",
                  "mild_vcf_records", "gene_variant_rows") + tuple(f"gene_variants_{c}" for c in CARRIER_CLASSES):
            out.write(f"{family}\t{m}\t{sev_str}\t{k}\t{counts[k]}\n")
        for s in severe:
            out.write(f"{family}\t{m}\t{sev_str}\tsevere_vcf_records_{s}\t{counts['severe_vcf_records_' + s]}\n")
    log(f"{m}: {len(hom_genes):,} homozygous genes, {len(chet_genes):,} compound-het genes, "
        f"{counts['genes_regulatory']:,} regulatory genes -> {len(gene_rows):,} passing genes; "
        f"{counts['mild_vcf_records']:,} mild VCF records")


if __name__ == "__main__":
    main()
