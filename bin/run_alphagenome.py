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

model.score_variant(...) returns a list of AnnData (one per variant
scorer); summarize_scores() flattens it with the package's own
variant_scorers.tidy_scores() rather than poking at AnnData internals
(quantiles live in .layers['quantiles'], track metadata in .var).
"""
import argparse
import gzip
import os
import sys
import warnings
from collections import Counter


def opener(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def get_info_field(info_str, key):
    prefix = f"{key}="
    for field in info_str.split(";"):
        if field.startswith(prefix):
            return field[len(prefix):]
    return ""


def parse_vcf(path):
    """Yields (vid, chrom, pos, ref, alt, source_model, gene, reason, samples,
    flipped) per record -- flipped is the ALLELE_FLIPPED INFO flag set by
    flip_ref_alt_alleles.py on lifted records where hg38's reference base is
    the candidate allele. vid is the VCF ID column -- build_alphagenome_input.py sets it to
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
                "ALLELE_FLIPPED" in info.split(";"),
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

# Scorers run per variant (keys of variant_scorers.RECOMMENDED_VARIANT_SCORERS,
# which for these are also the OutputType names). Each gets its own output
# column per tissue. Expression/splicing: RNA_SEQ (gene-level expression
# LFC), SPLICE_SITE_USAGE, SPLICE_JUNCTIONS. Chromatin: ATAC and DNASE
# (accessibility) and CHIP_HISTONE (histone marks) -- center-mask scorers,
# i.e. the local effect in a ~501bp (2001bp for histone) window around the
# variant. Order is also ontology-resolution priority (see below).
VARIANT_SCORER_KEYS = ["RNA_SEQ", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS", "DNASE", "ATAC", "CHIP_HISTONE"]

# output_metadata's per-modality attributes (output_metadata.rna_seq etc.) --
# only the modalities actually scored, so a tissue term can only resolve if
# some scored track carries it (otherwise its columns would be silently empty).
OUTPUT_METADATA_ATTRS = [k.lower() for k in VARIANT_SCORER_KEYS]

# Scorers whose score is signed (ALT vs REF direction). For ALLELE_FLIPPED
# variants the scored change is hg38 REF (= candidate allele) -> CHM13 REF,
# i.e. the reverse of carrying the candidate allele, so these are negated.
# The splice scorers report an unsigned magnitude of change, symmetric
# under swapping REF/ALT, so they're left as-is.
SIGNED_SCORERS = {"RNA_SEQ", "DNASE", "ATAC", "CHIP_HISTONE"}

# Scorers whose rows are per gene (tidy_scores gene_name/gene_id), feeding
# the impacted_genes list. The chromatin scorers are local to the variant,
# not attributed to a gene.
GENE_SCORERS = {"RNA_SEQ", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"}


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
        # A name match must not reuse a term another target already took --
        # otherwise e.g. "neuron" matches "motor neuron" and the neuron
        # column just duplicates motor_neuron.
        taken = {c for c in resolved.values() if c}
        found_curie = None
        for df in dfs:
            if found_curie:
                break
            if spec["cl_id"] is not None and (df["ontology_curie"] == spec["cl_id"]).any():
                found_curie = spec["cl_id"]
                break
            if "biosample_name" in df.columns:
                for needle in spec["name_contains"]:
                    hit = df[df["biosample_name"].str.contains(needle, case=False, na=False)
                             & ~df["ontology_curie"].isin(taken)]
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


def find_checkpoint(weights_dir, model_version="all_folds"):
    """Resolve the local checkpoint dir inside a pre-downloaded weights dir.

    Accepts either a kagglehub cache root (<dir>/models/google/alphagenome/
    jax/all_folds/<N>/, newest <N> wins) or the checkpoint dir itself.
    """
    model_root = os.path.join(weights_dir, "models", "google", "alphagenome", "jax", model_version)
    if os.path.isdir(model_root):
        versions = [d for d in os.listdir(model_root)
                    if d.isdigit() and os.path.isdir(os.path.join(model_root, d))]
        if not versions:
            sys.exit(f"[run_alphagenome] FATAL: no numbered version dirs under {model_root}")
        return os.path.join(model_root, max(versions, key=int))
    return weights_dir


# OrganismSettings field -> filename, matching the basenames of
# alphagenome_research's default Google Cloud Storage URLs (see
# bin/download_alphagenome_resources.sh, which fetches exactly these).
RESOURCE_FILES = {
    "gtf_feather_path": "gencode.v46.annotation.gtf.gz.feather",
    "pas_feather_path": "polyadb_human_v3_exon3_contiguous_gtfv46.feather",
    "splice_site_starts_feather_path": "gencode.v46.splice_sites_starts.feather",
    "splice_site_ends_feather_path": "gencode.v46.splice_sites_ends.feather",
    "calibration_path": "calibration_scores.pb",
}


def build_organism_settings(dna_model, dna_client, reference, resources_dir):
    """Start from alphagenome_research's default hg38 settings (all remote
    GCS URLs) and swap in local files: the FASTA is always the pipeline's
    own hg38 reference; the annotation/calibration files come from
    resources_dir when given. Any file missing from resources_dir keeps its
    remote default, so a partial download degrades to network access rather
    than silently disabling the scorers that depend on it."""
    import dataclasses
    human = dna_client.Organism.HOMO_SAPIENS
    settings = dict(dna_model.default_organism_settings())
    overrides = {"fasta_path": reference}
    if resources_dir:
        for field, fname in RESOURCE_FILES.items():
            path = os.path.join(resources_dir, fname)
            if os.path.isfile(path):
                overrides[field] = path
            else:
                sys.stderr.write(f"[run_alphagenome] WARNING: {path} not found -- {field} falls back to remote default\n")
    settings[human] = dataclasses.replace(settings[human], **overrides)
    return settings


def load_model(weights_dir, reference, resources_dir):
    """Load AlphaGenome for local inference against a pre-downloaded Kaggle weights cache.

    Calls dna_model.create() on the local checkpoint directly rather than
    create_from_kaggle(), which calls kagglehub.login() (an interactive
    prompt) whenever no Kaggle credentials are configured -- even when the
    weights are already cached -- and so crashes on a batch compute node.
    """
    from alphagenome_research.model import dna_model
    from alphagenome.models import dna_client
    checkpoint = find_checkpoint(weights_dir)
    sys.stderr.write(f"[run_alphagenome] loading checkpoint: {checkpoint}\n")
    organism_settings = build_organism_settings(dna_model, dna_client, reference, resources_dir)
    model = dna_model.create(checkpoint, organism_settings=organism_settings)
    return model, dna_client


def load_api_client():
    """Hosted DeepMind API backend (backup for when local GPU inference
    isn't available). Same score_variant()/output_metadata() interface as
    the local model, so everything downstream is shared. Sequence and
    annotations come from DeepMind's side -- -reference is only used for
    chromosome lengths (window clamping). Needs internet from the compute
    node and ALPHAGENOME_API_KEY in the environment (a Nextflow secret,
    see ALPHAGENOME_API in modules/alphagenome/main.nf)."""
    from alphagenome.models import dna_client
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        sys.exit("[run_alphagenome] FATAL: -backend api needs ALPHAGENOME_API_KEY set in the environment")
    sys.stderr.write("[run_alphagenome] using hosted AlphaGenome API\n")
    return dna_client.create(api_key), dna_client


# AlphaGenome requires one of a small set of fixed interval sizes. Default
# matches the ~100kb window documented in the quick-start guide; if the
# installed package exposes its own supported-lengths constant, that's
# preferred (see resize_interval()) so this default is only a fallback, not
# something silently trusted over the real model constraint.
DEFAULT_WINDOW = 131_072


def load_chrom_sizes(fai_path):
    sizes = {}
    with open(fai_path) as f:
        for line in f:
            fields = line.split("\t")
            sizes[fields[0]] = int(fields[1])
    return sizes


def resize_interval(genome_mod, dna_client, chrom_sizes, chrom, pos, window):
    # dict of name -> length, e.g. {'SEQUENCE_LENGTH_100KB': 131072, ...}
    supported = getattr(dna_client, "SUPPORTED_SEQUENCE_LENGTHS", None)
    if supported:
        window = min(supported.values(), key=lambda length: abs(length - window))
    # The model only accepts exactly-supported lengths, so near a chromosome
    # end shift the window inward rather than truncating it (a truncated
    # window breaks the model with incompatible-shape errors).
    start = pos - window // 2
    chrom_len = chrom_sizes.get(chrom)
    if chrom_len is not None:
        start = min(start, chrom_len - window)
    start = max(0, start)
    return genome_mod.Interval(chromosome=chrom, start=start, end=start + window)


def score_variant(model, genome_mod, dna_client, variant_scorers, chrom_sizes, chrom, pos, ref, alt, window):
    """Returns (scores, failed_scorer_notes). Tries all scorers in one pass;
    if that raises (e.g. SpliceJunctionScorer errors when no junctions fall
    in the window), retries each scorer alone so one scorer's failure
    doesn't throw away the others' results. Raises only if every scorer fails."""
    interval = resize_interval(genome_mod, dna_client, chrom_sizes, chrom, pos, window)
    variant = genome_mod.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)

    def run(scorers):
        return model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=scorers,
            organism=dna_client.Organism.HOMO_SAPIENS,
        )

    try:
        return run(variant_scorers), []
    except Exception:
        if len(variant_scorers) == 1:
            raise

    scores, failed, last_err = [], [], None
    for key, scorer in zip(VARIANT_SCORER_KEYS, variant_scorers):
        try:
            scores.extend(run([scorer]))
        except Exception as e:
            failed.append(f"{key}: {' '.join(str(e).split())}")
            last_err = e
    if not scores:
        raise last_err
    return scores, failed


def score_columns():
    """One column per tissue x scorer, e.g. motor_neuron_RNA_SEQ."""
    return [f"{tissue}_{scorer}" for tissue in ONTOLOGY_TARGETS for scorer in VARIANT_SCORER_KEYS]


def gene_label(row):
    # tidy_scores leaves missing gene fields as NaN, not None
    return next((g for g in (row.get("gene_name"), row.get("gene_id")) if isinstance(g, str) and g), "")


def impacted_genes(df, ontology_resolved, candidate_genes, threshold):
    """'GENE*(tissue.SCORER=q|...);GENE(...)' for every gene with a
    gene-level score at/above threshold in a target tissue -- * marks the
    variant's candidate gene(s) from the modifier search. Within a gene,
    the strongest score per tissue x scorer, strongest first; candidate
    genes first, then by strongest score. Returns (string, n_genes)."""
    hits = {}  # gene -> {(tissue, scorer): q}
    for tissue, curie in ontology_resolved.items():
        if not curie:
            continue
        sub = df[(df["ontology_curie"] == curie) & df["scorer"].isin(GENE_SCORERS)
                 & (df["quantile_score"].abs() >= threshold)]
        for _, r in sub.iterrows():
            gene = gene_label(r)
            if not gene:
                continue
            key = (tissue, r["scorer"])
            q = float(r["quantile_score"])
            cur = hits.setdefault(gene, {})
            if key not in cur or abs(q) > abs(cur[key]):
                cur[key] = q

    def gene_sort(item):
        gene, scores = item
        return (gene not in candidate_genes, -max(abs(q) for q in scores.values()))

    parts = []
    for gene, scores in sorted(hits.items(), key=gene_sort):
        detail = "|".join(f"{t}.{sc}={q:.3f}" for (t, sc), q in sorted(scores.items(), key=lambda kv: -abs(kv[1])))
        parts.append(f"{gene}{'*' if gene in candidate_genes else ''}({detail})")
    return ";".join(parts), len(hits)


def summarize_scores(scores, ontology_resolved, variant_scorers_lib, flipped=False, candidate_genes=(), threshold=0.99):
    """
    scores: list[AnnData] from model.score_variant(), one per scorer.
    Returns {tissue_SCORER: signed quantile_score} -- the largest-magnitude
    quantile across that tissue's tracks for that scorer (all genes in the
    window for gene-level scorers, all marks for CHIP_HISTONE) -- plus the
    single strongest motor-neuron hit across all scorers, and the
    impacted_genes list (see impacted_genes()). flipped: negate the
    SIGNED_SCORERS first, so every output describes the candidate allele.
    """
    result = {col: None for col in score_columns()}
    result.update({"motor_neuron_max_abs": None, "best_output_type": None, "best_track": None,
                   "impacted_genes": "", "n_impacted_genes": 0})

    df = variant_scorers_lib.tidy_scores(scores)
    if df is None or "quantile_score" not in df.columns or "ontology_curie" not in df.columns:
        return result
    df = df.dropna(subset=["quantile_score"])
    # output_type may come back as an OutputType enum or "OutputType.RNA_SEQ"
    # -- normalize to the bare name so it matches VARIANT_SCORER_KEYS
    df = df.assign(scorer=df["output_type"].map(lambda o: str(getattr(o, "name", o)).split(".")[-1].upper()))
    if flipped:
        signed = df["scorer"].isin(SIGNED_SCORERS)
        df.loc[signed, "quantile_score"] = -df.loc[signed, "quantile_score"]

    for tissue, curie in ontology_resolved.items():
        if not curie:
            continue
        sub = df[df["ontology_curie"] == curie]
        for scorer, grp in sub.groupby("scorer"):
            col = f"{tissue}_{scorer}"
            if col in result:
                result[col] = float(grp.loc[grp["quantile_score"].abs().idxmax(), "quantile_score"])
        if tissue == "motor_neuron" and not sub.empty:
            best = sub.loc[sub["quantile_score"].abs().idxmax()]
            result["motor_neuron_max_abs"] = abs(float(best["quantile_score"]))
            result["best_output_type"] = best["scorer"]
            gene = gene_label(best)
            track = str(best.get("track_name", ""))
            result["best_track"] = f"{track}|{gene}" if gene else track

    result["impacted_genes"], result["n_impacted_genes"] = impacted_genes(df, ontology_resolved, candidate_genes, threshold)
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
    p.add_argument("-reference", required=True, help="hg38 reference FASTA (with .fai alongside) -- used as the model's sequence source")
    p.add_argument("-backend", choices=["local", "api"], default="local",
                   help="local = GPU inference against -weights; api = hosted DeepMind API "
                        "(needs ALPHAGENOME_API_KEY and internet; -weights/-resources ignored)")
    p.add_argument("-weights", help="Pre-downloaded local Kaggle AlphaGenome weights directory (required for -backend local)")
    p.add_argument("-resources", default="",
                   help="Directory of pre-downloaded AlphaGenome annotation/calibration files "
                        "(bin/download_alphagenome_resources.sh). Unset = fetch them from GCS at startup.")
    p.add_argument("-window", type=int, default=DEFAULT_WINDOW, help="Sequence context window (bp), resized to the nearest AlphaGenome-supported length")
    p.add_argument("-quantile_threshold", type=float, default=0.99,
                   help="|quantile_score| at/above which a gene counts as impacted (impacted_genes) and "
                        "large_impact_motor_neuron is flagged")
    p.add_argument("-output", required=True, help="Output TSV path")
    args = p.parse_args()

    # alphagenome's own tidy_scores() pd.concat trips this on every call
    # with newer pandas -- harmless, and it floods the log
    warnings.filterwarnings("ignore", category=FutureWarning,
                            message="The behavior of DataFrame concatenation with empty or all-NA entries")

    if args.backend == "api":
        model, dna_client = load_api_client()
    else:
        if not args.weights:
            p.error("-weights is required for -backend local")
        model, dna_client = load_model(args.weights, args.reference, args.resources)
    from alphagenome.data import genome as genome_mod
    from alphagenome.models import variant_scorers as variant_scorers_lib
    variant_scorers = [variant_scorers_lib.RECOMMENDED_VARIANT_SCORERS[k] for k in VARIANT_SCORER_KEYS]
    ontology_resolved = resolve_ontology_terms(model, dna_client)
    chrom_sizes = load_chrom_sizes(args.reference + ".fai")

    # Reconcile against the full pre-liftover candidate list, keyed by VCF
    # ID (CrossMap preserves it unchanged through liftover; CHROM/POS/REF/ALT
    # themselves change, so they can't be used to match records across the
    # two files). This is what guarantees every recessive/dominant-modifier
    # PASS variant gets exactly one row below, even ones that failed to
    # lift over or failed AlphaGenome scoring -- neither case silently
    # drops a variant from the final table.
    candidates_by_vid = {}
    for vid, chrom, pos, ref, alt, source_model, gene, reason, samples, _flipped in parse_vcf(args.candidate_vcf):
        candidates_by_vid[vid] = {
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "source_model": source_model, "gene": gene, "reason": reason, "samples": samples,
        }

    lifted_by_vid = {}
    for vid, chrom, pos, ref, alt, *_, flipped in parse_vcf(args.vcf):
        lifted_by_vid[vid] = (chrom, pos, ref, alt, flipped)

    columns = [
        "VID", "ORIG_CHROM", "ORIG_POS", "ORIG_REF", "ORIG_ALT",
        "HG38_CHROM", "HG38_POS", "HG38_REF", "HG38_ALT",
        "SOURCE_MODEL", "GENE", "REASON", "SAMPLES", "status", "allele_flipped", "note",
        *score_columns(),
        "n_impacted_genes", "impacted_genes",
        "motor_neuron_max_abs", "best_output_type", "best_track", "large_impact_motor_neuron",
    ]

    def fmt(v):
        if v is None:
            return ""
        return f"{v:.4f}" if isinstance(v, float) else str(v)

    def write_row(out, row):
        out.write("\t".join(fmt(row.get(c)) for c in columns) + "\n")

    n_scored = 0
    n_flipped = 0
    partial_counts = Counter()
    n_liftover_failed = 0
    n_scoring_failed = 0
    with open(args.output, "w") as out:
        out.write("\t".join(columns) + "\n")
        for vid in sorted(candidates_by_vid):
            meta = candidates_by_vid[vid]
            row = {
                "VID": vid, "ORIG_CHROM": meta["chrom"], "ORIG_POS": meta["pos"],
                "ORIG_REF": meta["ref"], "ORIG_ALT": meta["alt"],
                "SOURCE_MODEL": meta["source_model"], "GENE": meta["gene"],
                "REASON": meta["reason"], "SAMPLES": meta["samples"],
            }

            if vid not in lifted_by_vid:
                n_liftover_failed += 1
                write_row(out, {**row, "status": "liftover_failed"})
                continue

            chrom, pos, ref, alt, flipped = lifted_by_vid[vid]
            row.update({"HG38_CHROM": chrom, "HG38_POS": pos, "HG38_REF": ref, "HG38_ALT": alt,
                        "allele_flipped": flipped})
            try:
                scores, failed_scorers = score_variant(model, genome_mod, dna_client, variant_scorers, chrom_sizes,
                                                       chrom, pos, ref, alt, args.window)
                summary = summarize_scores(scores, ontology_resolved, variant_scorers_lib, flipped=flipped,
                                           candidate_genes=set(filter(None, meta["gene"].split(","))),
                                           threshold=args.quantile_threshold)
            except Exception as e:
                sys.stderr.write(f"[run_alphagenome] WARNING: scoring failed for {vid} ({chrom}:{pos}:{ref}>{alt}): {e}\n")
                n_scoring_failed += 1
                note = " ".join(str(e).split())  # tabs/newlines would otherwise break the TSV row
                write_row(out, {**row, "status": "scoring_failed", "note": note})
                continue

            mn_max = summary["motor_neuron_max_abs"]
            notes = []
            if flipped:
                # signed scores already inverted inside summarize_scores()
                n_flipped += 1
                notes.append("hg38 ref = candidate allele; scored hg38ref>CHM13ref, signed scores inverted")
            if failed_scorers:
                # one summary line at the end instead of a log line per variant
                partial_counts.update(f.split(":")[0] for f in failed_scorers)
                notes.append("partial; " + "; ".join(failed_scorers))
            note = " | ".join(notes)

            write_row(out, {
                **row, **summary, "status": "scored", "note": note,
                "large_impact_motor_neuron": mn_max is not None and mn_max >= args.quantile_threshold,
            })
            n_scored += 1
            if n_scored % 500 == 0:
                sys.stderr.write(f"[run_alphagenome] progress: {n_scored} scored / {len(candidates_by_vid)} candidates\n")

    for scorer, n in sorted(partial_counts.items()):
        sys.stderr.write(f"[run_alphagenome] {scorer} failed for {n} variant(s) scored with the remaining scorers "
                         f"(reason in each row's note column)\n")
    sys.stderr.write(
        f"[run_alphagenome] {len(candidates_by_vid)} total candidate(s): {n_scored} scored "
        f"({n_flipped} allele-flipped), "
        f"{n_liftover_failed} failed liftover, {n_scoring_failed} failed scoring -> {args.output}\n"
    )


if __name__ == "__main__":
    main()
