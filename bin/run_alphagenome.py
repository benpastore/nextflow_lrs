#!/usr/bin/env python3
"""
Score AlphaGenome variant-effect predictions (RNA-seq expression + splicing)
for the recessive/dominant-modifier candidate variants, restricted to
motor-neuron-relevant tissue (this is an SMA study).

Local inference only -- against a pre-downloaded Kaggle AlphaGenome weights
directory (google-deepmind/alphagenome_research), not the hosted DeepMind
API the public `alphagenome` pip package defaults to by default.

Input: the CHM13->hg38-lifted candidate VCF (build_alphagenome_input.py +
LIFTOVER_ALPHAGENOME_VARIANTS, see modules/alphagenome/main.nf) -- its
SOURCE_MODEL/GENE/REASON/SAMPLES INFO fields are carried straight through
to the output table.

*** First-real-run risk, flagged deliberately (see the plan this was built
from): the exact shape of what model.score_variant(...) returns (which
AnnData attribute holds quantile_score -- .layers vs a per-scorer .X vs
something else) is inferred from documentation, not from having actually
run this. summarize_scores() below is the one function to adjust if the
real return shape differs -- everything else (ontology resolution, VCF
I/O, CLI) doesn't depend on that detail and shouldn't need changes.
"""
import argparse
import gzip
import os
import sys


def opener(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def get_info_field(info_str, key):
    prefix = f"{key}="
    for field in info_str.split(";"):
        if field.startswith(prefix):
            return field[len(prefix):]
    return ""


def parse_vcf(path):
    """Yields (vid, chrom, pos, ref, alt, source_model, gene, reason, samples) per
    record. vid is the VCF ID column -- build_alphagenome_input.py sets it to
    the original CHM13 chrom_pos_ref_alt, and CrossMap preserves the ID
    column through liftover unchanged, so it's what ties a lifted hg38
    record back to its original CHM13 candidate (CHROM/POS/REF/ALT
    themselves change across liftover, so they can't be used for that)."""
    with opener(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, vid, ref, alt, _qual, _filt, info = fields[:8]
            yield (
                vid, chrom, int(pos), ref, alt,
                get_info_field(info, "SOURCE_MODEL"),
                get_info_field(info, "GENE"),
                get_info_field(info, "REASON"),
                get_info_field(info, "SAMPLES"),
            )


# Ontology terms scored, in priority order -- motor neuron is the primary,
# non-negotiable target for this SMA study; spinal cord / general neuron /
# skeletal muscle are secondary context columns (muscle included since SMA's
# clinical phenotype is driven as much by denervation atrophy of skeletal
# muscle as by the motor neuron loss itself). Resolved against whatever
# CURIEs are actually present in this AlphaGenome release via
# resolve_ontology_terms() below -- these are *labels* to search
# output_metadata for, not hardcoded CURIEs, since the exact CURIE set can
# vary by release. UBERON:0001134 (skeletal muscle tissue) is confirmed as
# a documented example UBERON term in AlphaGenome's own docs.
ONTOLOGY_TARGETS = {
    "motor_neuron": {"cl_id": "CL:0000100", "name_contains": ["motor neuron"]},
    "spinal_cord": {"cl_id": None, "name_contains": ["spinal cord"]},
    "neuron": {"cl_id": "CL:0000540", "name_contains": ["neuron"]},
    "muscle": {"cl_id": "UBERON:0001134", "name_contains": ["skeletal muscle", "muscle"]},
}

# output_metadata's per-modality attributes, per AlphaGenome's documented
# OutputType list -- accessed as output_metadata.<name> (confirmed pattern:
# "output_metadata.rna_seq" in the quick-start docs).
OUTPUT_METADATA_ATTRS = [
    "atac", "cage", "dnase", "rna_seq", "chip_histone", "chip_tf",
    "splice_sites", "splice_site_usage", "splice_junctions", "contact_maps", "procap",
]


def resolve_ontology_terms(model, dna_client):
    """
    Query this AlphaGenome release's actual supported ontology terms and
    resolve ONTOLOGY_TARGETS to real CURIEs present in it. Hard-fails only
    if motor neuron itself can't be found (the one non-negotiable term for
    this SMA study); logs a warning and drops any other target that's
    missing rather than failing the whole run over it.
    """
    metadata = model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)

    dfs = []
    for attr in OUTPUT_METADATA_ATTRS:
        df = getattr(metadata, attr, None)
        if df is not None and hasattr(df, "columns") and "ontology_curie" in df.columns:
            dfs.append(df)

    resolved = {}
    for key, spec in ONTOLOGY_TARGETS.items():
        found_curie = None
        for df in dfs:
            if found_curie:
                break
            if spec["cl_id"] is not None and (df["ontology_curie"] == spec["cl_id"]).any():
                found_curie = spec["cl_id"]
                break
            if "biosample_name" in df.columns:
                for needle in spec["name_contains"]:
                    hit = df[df["biosample_name"].str.contains(needle, case=False, na=False)]
                    if len(hit):
                        found_curie = hit.iloc[0]["ontology_curie"]
                        break
        resolved[key] = found_curie
        if found_curie:
            sys.stderr.write(f"[run_alphagenome] ontology term resolved: {key} -> {found_curie}\n")
        else:
            sys.stderr.write(f"[run_alphagenome] WARNING: no ontology term found for '{key}' in this AlphaGenome release\n")

    if not resolved["motor_neuron"]:
        sys.exit(
            "[run_alphagenome] FATAL: motor neuron (CL:0000100) not found in this AlphaGenome "
            "release's output_metadata -- cannot score against the target tissue for this study."
        )

    return resolved


def load_model(weights_dir):
    """Load AlphaGenome for local inference against a pre-downloaded Kaggle weights cache."""
    os.environ.setdefault("KAGGLEHUB_CACHE", weights_dir)
    from alphagenome_research.model import dna_model
    from alphagenome.models import dna_client
    model = dna_model.create_from_kaggle("all_folds")
    return model, dna_client


# AlphaGenome requires one of a small set of fixed interval sizes. Default
# matches the ~100kb window documented in the quick-start guide; if the
# installed package exposes its own supported-lengths constant, that's
# preferred (see resize_interval()) so this default is only a fallback, not
# something silently trusted over the real model constraint.
DEFAULT_WINDOW = 131_072


def resize_interval(genome_mod, dna_client, chrom, pos, window):
    supported = getattr(dna_client, "SUPPORTED_SEQUENCE_LENGTHS", None)
    if supported:
        window = min(supported, key=lambda length: abs(length - window))
    half = window // 2
    return genome_mod.Interval(chromosome=chrom, start=max(0, pos - half), end=pos + half)


def score_variant(model, genome_mod, dna_client, ontology_curies, chrom, pos, ref, alt, window):
    interval = resize_interval(genome_mod, dna_client, chrom, pos, window)
    variant = genome_mod.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)
    outputs = [
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.SPLICE_SITE_USAGE,
        dna_client.OutputType.SPLICE_JUNCTIONS,
    ]
    return model.score_variant(
        interval=interval,
        variant=variant,
        ontology_terms=ontology_curies,
        requested_outputs=outputs,
    )


def summarize_scores(scores, ontology_resolved):
    """
    scores: AnnData-like (genes x tracks) with track metadata in .var
    (including ontology_curie) and a quantile_score layer -- see the
    first-real-run-risk note at the top of this file. Returns, per
    resolved ontology term, the largest-magnitude quantile_score across
    that term's tracks.
    """
    result = {
        "motor_neuron_quantile_score": None,
        "spinal_cord_quantile_score": None,
        "neuron_quantile_score": None,
        "muscle_quantile_score": None,
        "best_output_type": None,
        "best_track": None,
    }

    var = getattr(scores, "var", None)
    layers = getattr(scores, "layers", None)
    if var is None or layers is None or "quantile_score" not in layers or "ontology_curie" not in var.columns:
        return result

    quantile = layers["quantile_score"]
    best_abs = -1.0
    for key, curie in ontology_resolved.items():
        if not curie:
            continue
        track_mask = (var["ontology_curie"] == curie).values
        if not track_mask.any():
            continue
        sub = quantile[:, track_mask]
        if sub.size == 0:
            continue
        flat = sub.flatten()
        max_flat_idx = abs(flat).argmax()
        max_val = float(flat[max_flat_idx])
        result[f"{key}_quantile_score"] = max_val
        if key == "motor_neuron" and abs(max_val) > best_abs:
            best_abs = abs(max_val)
            gene_idx, track_col = divmod(max_flat_idx, sub.shape[1])
            track_idx = track_mask.nonzero()[0][track_col]
            result["best_output_type"] = str(var.iloc[track_idx].get("output_type", ""))
            result["best_track"] = str(var.iloc[track_idx].get("name", var.index[track_idx]))

    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-candidate_vcf", required=True,
                   help="Original candidate VCF (build_alphagenome_input.py's output, before any "
                        "liftover) -- every variant in this file gets exactly one output row, "
                        "whether or not it lifted/scored successfully, so the full "
                        "recessive/dominant-modifier PASS set is always fully represented here. "
                        "In T2T mode this is CHM13 coordinates; in HG38 mode it's identical to "
                        "-vcf (see params.genome_version / the alphagenome subworkflow).")
    p.add_argument("-vcf", required=True,
                   help="Candidate VCF in hg38 coordinates -- lifted (T2T mode) or the same file "
                        "as -candidate_vcf (HG38 mode, liftover skipped).")
    p.add_argument("-reference", required=True, help="hg38 reference FASTA")
    p.add_argument("-weights", required=True, help="Pre-downloaded local Kaggle AlphaGenome weights directory")
    p.add_argument("-window", type=int, default=DEFAULT_WINDOW, help="Sequence context window (bp), resized to the nearest AlphaGenome-supported length")
    p.add_argument("-quantile_threshold", type=float, default=0.9, help="|quantile_score| at/above which a variant is flagged large_impact_motor_neuron")
    p.add_argument("-output", required=True, help="Output TSV path")
    args = p.parse_args()

    model, dna_client = load_model(args.weights)
    from alphagenome.models import genome as genome_mod
    ontology_resolved = resolve_ontology_terms(model, dna_client)
    ontology_curies = [c for c in ontology_resolved.values() if c]

    # Reconcile against the full pre-liftover candidate list, keyed by VCF
    # ID (CrossMap preserves it unchanged through liftover; CHROM/POS/REF/ALT
    # themselves change, so they can't be used to match records across the
    # two files). This is what guarantees every recessive/dominant-modifier
    # PASS variant gets exactly one row below, even ones that failed to
    # lift over or failed AlphaGenome scoring -- neither case silently
    # drops a variant from the final table.
    candidates_by_vid = {}
    for vid, chrom, pos, ref, alt, source_model, gene, reason, samples in parse_vcf(args.candidate_vcf):
        candidates_by_vid[vid] = {
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "source_model": source_model, "gene": gene, "reason": reason, "samples": samples,
        }

    lifted_by_vid = {}
    for vid, chrom, pos, ref, alt, *_ in parse_vcf(args.vcf):
        lifted_by_vid[vid] = (chrom, pos, ref, alt)

    columns = [
        "VID", "ORIG_CHROM", "ORIG_POS", "ORIG_REF", "ORIG_ALT",
        "HG38_CHROM", "HG38_POS", "HG38_REF", "HG38_ALT",
        "SOURCE_MODEL", "GENE", "REASON", "SAMPLES", "status", "note",
        "motor_neuron_quantile_score", "spinal_cord_quantile_score", "neuron_quantile_score", "muscle_quantile_score",
        "best_output_type", "best_track", "large_impact_motor_neuron",
    ]

    n_scored = 0
    n_liftover_failed = 0
    n_scoring_failed = 0
    with open(args.output, "w") as out:
        out.write("\t".join(columns) + "\n")
        for vid in sorted(candidates_by_vid):
            meta = candidates_by_vid[vid]
            base_row = [
                vid, meta["chrom"], str(meta["pos"]), meta["ref"], meta["alt"],
            ]

            if vid not in lifted_by_vid:
                n_liftover_failed += 1
                out.write("\t".join(base_row + ["", "", "", ""] + [
                    meta["source_model"], meta["gene"], meta["reason"], meta["samples"],
                    "liftover_failed", "", "", "", "", "", "", "", "",
                ]) + "\n")
                continue

            chrom, pos, ref, alt = lifted_by_vid[vid]
            base_row += [chrom, str(pos), ref, alt, meta["source_model"], meta["gene"], meta["reason"], meta["samples"]]
            try:
                scores = score_variant(model, genome_mod, dna_client, ontology_curies, chrom, pos, ref, alt, args.window)
                summary = summarize_scores(scores, ontology_resolved)
            except Exception as e:
                sys.stderr.write(f"[run_alphagenome] WARNING: scoring failed for {vid} ({chrom}:{pos}:{ref}>{alt}): {e}\n")
                n_scoring_failed += 1
                note = " ".join(str(e).split())  # tabs/newlines would otherwise break the TSV row
                out.write("\t".join(base_row + ["scoring_failed", note, "", "", "", "", "", "", ""]) + "\n")
                continue

            mn_score = summary["motor_neuron_quantile_score"]
            large_impact = mn_score is not None and abs(mn_score) >= args.quantile_threshold

            out.write("\t".join(base_row + [
                "scored", "",
                "" if summary["motor_neuron_quantile_score"] is None else f"{summary['motor_neuron_quantile_score']:.4f}",
                "" if summary["spinal_cord_quantile_score"] is None else f"{summary['spinal_cord_quantile_score']:.4f}",
                "" if summary["neuron_quantile_score"] is None else f"{summary['neuron_quantile_score']:.4f}",
                "" if summary["muscle_quantile_score"] is None else f"{summary['muscle_quantile_score']:.4f}",
                summary["best_output_type"] or "",
                summary["best_track"] or "",
                str(large_impact),
            ]) + "\n")
            n_scored += 1

    sys.stderr.write(
        f"[run_alphagenome] {len(candidates_by_vid)} total candidate(s): {n_scored} scored, "
        f"{n_liftover_failed} failed liftover, {n_scoring_failed} failed scoring -> {args.output}\n"
    )


if __name__ == "__main__":
    main()
