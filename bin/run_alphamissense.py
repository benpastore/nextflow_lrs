#!/usr/bin/env python3
"""
AlphaMissense pathogenicity for the recessive/dominant-modifier candidates
that land in protein-coding sequence.

AlphaMissense (Cheng et al. 2023, Science) precomputes a pathogenicity score
for every possible missense SNV in hg38 canonical transcripts
(AlphaMissense_hg38.tsv.gz) or all isoforms
(AlphaMissense_isoforms_hg38.tsv.gz). The table only contains missense
changes, so a lifted SNV with no row there is, by construction, not a
missense change in an AlphaMissense transcript (synonymous, stop
gain/loss, splice, UTR, intronic or intergenic) -- the lookup doubles as the
coding-missense filter. Per candidate (one row each, same VID as
build_alphagenome_input.py / the final report):

  status
    scored                 >=1 AlphaMissense row; columns describe the
                           highest-scoring transcript, am_all_hits lists all
    not_missense           lifted SNV, no AlphaMissense row (see above)
    not_snv                indel / MNV / SV -- AlphaMissense scores SNVs only
    hg38_ref_is_candidate  CrossMap Fail(REF==ALT): hg38's reference base IS
                           the candidate allele, so in hg38 the candidate
                           residue is the reference residue -- nothing to score
    liftover_failed        no hg38 coordinates (note = CrossMap's reason)

  am_class follows the AlphaMissense paper's calibrated cut-offs:
  likely_benign < 0.34 <= ambiguous <= 0.564 < likely_pathogenic.

The lookup uses tabix (pysam) when <alphamissense>.tbi exists (see
bin/index_alphamissense.py), else streams the whole table once.
"""
import argparse
import csv
import gzip
import os
import sys
from collections import Counter, defaultdict

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BENIGN_MAX, PATHOGENIC_MIN = 0.34, 0.564

COLUMNS = [
    "VID", "CHROM", "POS", "REF", "ALT", "GENE", "SOURCE_MODEL",
    "HG38_CHROM", "HG38_POS", "HG38_REF", "HG38_ALT",
    "status", "note",
    "am_pathogenicity", "am_class", "am_protein_variant", "am_transcript", "am_uniprot",
    "am_gene_mean_pathogenicity", "am_n_hits", "am_all_hits",
]


def log(msg):
    sys.stderr.write(f"[run_alphamissense] {msg}\n")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def present(path):
    return bool(path) and os.path.basename(path) != "NO_FILE" and os.path.exists(path)


def info_value(info, key):
    for field in info.split(";"):
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return ""


def read_vcf(path):
    """VID -> (chrom, pos, ref, alt, info) for every record."""
    records = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8:
                records[fields[2]] = (fields[0], fields[1], fields[3], fields[4], fields[7])
    return records


def read_unmap(path):
    """VID -> CrossMap failure reason (last column of the .unmap)."""
    reasons = {}
    if present(path):
        with opener(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 3:
                    reasons[fields[2]] = fields[-1]
    return reasons


def classify(score):
    if score < BENIGN_MAX:
        return "likely_benign"
    return "likely_pathogenic" if score > PATHOGENIC_MIN else "ambiguous"


def header_columns(lines):
    for line in lines:
        if line.startswith("#CHROM"):
            return line.lstrip("#").rstrip("\n").split("\t")
    sys.exit("run_alphamissense.py: no #CHROM header in the AlphaMissense table")


def parse_hit(fields, col):
    score = float(fields[col["am_pathogenicity"]])
    return {
        "score": score,
        "class": fields[col["am_class"]] if "am_class" in col else classify(score),
        "protein_variant": fields[col["protein_variant"]],
        "transcript": fields[col["transcript_id"]],
        "uniprot": fields[col["uniprot_id"]] if "uniprot_id" in col else "",
    }


def lookup(path, wanted):
    """{(chrom, pos, ref, alt): [hit, ...]} for the wanted hg38 SNVs."""
    hits = defaultdict(list)
    if os.path.exists(path + ".tbi"):
        try:
            import pysam
        except ImportError:
            pysam = None
            log("pysam unavailable -- streaming the table instead of using its tabix index")
        if pysam is not None:
            tbx = pysam.TabixFile(path)
            col = {c: i for i, c in enumerate(header_columns(tbx.header))}
            contigs = set(tbx.contigs)
            for chrom, pos, ref, alt in sorted(wanted):
                if chrom not in contigs:
                    continue
                for line in tbx.fetch(chrom, int(pos) - 1, int(pos)):
                    fields = line.split("\t")
                    if fields[1] == pos and fields[2] == ref and fields[3] == alt:
                        hits[(chrom, pos, ref, alt)].append(parse_hit(fields, col))
            log(f"tabix lookup: {len(hits)}/{len(wanted)} lifted SNV(s) have an AlphaMissense row")
            return hits

    n = 0
    header = []
    col = None
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                header.append(line)
                continue
            if col is None:
                col = {c: i for i, c in enumerate(header_columns(header))}
            n += 1
            key = tuple(line.split("\t", 4)[:4])
            if key in wanted:
                hits[key].append(parse_hit(line.rstrip("\n").split("\t"), col))
    log(f"streamed {n:,} rows: {len(hits)}/{len(wanted)} lifted SNV(s) have an AlphaMissense row")
    return hits


def read_gene_means(path):
    """transcript_id -> mean_am_pathogenicity (AlphaMissense_gene_hg38.tsv.gz),
    keyed both with and without the version suffix."""
    means = {}
    if not present(path):
        return means
    with opener(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("transcript_id"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                means[fields[0]] = fields[1]
                means.setdefault(fields[0].split(".")[0], fields[1])
    return means


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-candidate_vcf", required=True, help="CHM13 (or native hg38) candidates, ID = VID (build_alphagenome_input.py)")
    p.add_argument("-hg38_vcf", required=True, help="Candidates lifted to hg38, ID preserved (same file as -candidate_vcf if already hg38)")
    p.add_argument("-unmap", default="", help="CrossMap .unmap for the liftover")
    p.add_argument("-alphamissense", required=True, help="AlphaMissense_hg38.tsv.gz (tabix-indexed copy preferred)")
    p.add_argument("-gene_tsv", default="", help="AlphaMissense_gene_hg38.tsv.gz (optional, per-transcript mean pathogenicity)")
    p.add_argument("-output", required=True)
    args = p.parse_args()

    candidates = read_vcf(args.candidate_vcf)
    lifted = read_vcf(args.hg38_vcf)
    unmap = read_unmap(args.unmap)
    gene_means = read_gene_means(args.gene_tsv)

    wanted = {rec[:4] for rec in lifted.values()
              if len(rec[2]) == 1 and len(rec[3]) == 1 and "ALLELE_FLIPPED" not in rec[4].split(";")}
    hits = lookup(args.alphamissense, wanted)

    status_counts = Counter()
    class_counts = Counter()
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for vid, (chrom, pos, ref, alt, info) in candidates.items():
            row = {"VID": vid, "CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt,
                   "GENE": info_value(info, "GENE"), "SOURCE_MODEL": info_value(info, "SOURCE_MODEL")}

            if len(ref) != 1 or len(alt) != 1 or alt.startswith("<"):
                row["status"] = "not_snv"
            elif vid in lifted:
                h_chrom, h_pos, h_ref, h_alt, h_info = lifted[vid]
                row.update({"HG38_CHROM": h_chrom, "HG38_POS": h_pos, "HG38_REF": h_ref, "HG38_ALT": h_alt})
                # CrossMap complements both alleles on minus-strand chain hits
                minus = h_alt.upper() != alt.upper()
                if "ALLELE_FLIPPED" in h_info.split(";"):
                    row["status"] = "hg38_ref_is_candidate"
                elif (h_chrom, h_pos, h_ref, h_alt) in hits:
                    vhits = sorted(hits[(h_chrom, h_pos, h_ref, h_alt)], key=lambda h: -h["score"])
                    best = vhits[0]
                    row.update({
                        "status": "scored",
                        "am_pathogenicity": f"{best['score']:.4f}",
                        "am_class": best["class"],
                        "am_protein_variant": best["protein_variant"],
                        "am_transcript": best["transcript"],
                        "am_uniprot": best["uniprot"],
                        "am_gene_mean_pathogenicity": gene_means.get(best["transcript"],
                                                                     gene_means.get(best["transcript"].split(".")[0], "")),
                        "am_n_hits": len(vhits),
                        "am_all_hits": "|".join(f"{h['transcript']}:{h['protein_variant']}={h['score']:.4f}"
                                                for h in vhits),
                    })
                    class_counts[best["class"]] += 1
                else:
                    row["status"] = "not_missense"
                if h_ref.upper() != (ref.translate(COMPLEMENT) if minus else ref).upper():
                    row["note"] = "hg38_ref_differs"  # scored as hg38 REF > candidate ALT
            elif unmap.get(vid) == "Fail(REF==ALT)":
                row["status"] = "hg38_ref_is_candidate"
            else:
                row["status"] = "liftover_failed"
                row["note"] = unmap.get(vid, "not_in_liftover_output")

            status_counts[row["status"]] += 1
            writer.writerow(row)

    log(f"{len(candidates)} candidate(s): " + ", ".join(f"{n} {s}" for s, n in status_counts.most_common()))
    if class_counts:
        log("scored: " + ", ".join(f"{n} {c}" for c, n in class_counts.most_common()))


if __name__ == "__main__":
    main()
