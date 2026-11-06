#!/usr/bin/env python3
"""
Final per-variant report for the recessive/dominant-modifier PASS
candidates: one row per unique variant (chrom/pos/ref/alt, union of both
models), joining everything the pipeline knows about it:

  - modifier-search info (bin/recessive_modifier_consolidated.py,
    bin/dominant_modifier_consolidated.py): model(s), gene(s), zygosity,
    reason, calling evidence, the mild sibs it passed in, gnomAD AF_nfe
  - cohort carriers: every sample in the VCF manifest whose consolidated+
    snpEff VCF has the variant called (non-ref GT from its preferred caller,
    same CALLER_PRIORITY the modifier search trusts), split mild/severe by
    family.json phenotype
  - snpEff: the most severe ANN entry for the variant's candidate gene(s)
  - AlphaMissense (optional): bin/run_alphamissense.py's per-candidate
    result -- status (scored / not_missense / not_snv / ...), the
    highest-scoring transcript's pathogenicity + class, and whether its
    protein change matches one of snpEff's (CHM13 RefSeq) HGVS.p calls for
    the variant (alphamissense_snpeff_agree)
  - AlphaGenome (optional): every tissue x scorer column from
    bin/run_alphagenome.py, from the local run, falling back per variant to
    the API run when the local one didn't score it (alphagenome_source says
    which)

Optional inputs that are absent (or the pipeline's NO_FILE placeholder)
just leave their columns empty.

Also writes a variant funnel (-funnel_output): how many variants enter
and leave each step (modifier search -> AlphaGenome input -> liftover ->
AlphaGenome / AlphaMissense -> report), counted from each step's actual output files rather
than by re-applying its filters, so it reflects what really happened.
"""
import argparse
import csv
import json
import os
import sys
from collections import defaultdict

from _modifier_common import (
    CALLER_PRIORITY, classify_gt, get_ann, get_info_field, load_vcf_manifest, logger, opener,
)

IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}

# run_alphagenome.py columns that are identity/provenance, not results --
# everything else in its header is carried into the report as alphagenome_*
ALPHAGENOME_ID_COLUMNS = {
    "VID", "ORIG_CHROM", "ORIG_POS", "ORIG_REF", "ORIG_ALT",
    "HG38_CHROM", "HG38_POS", "HG38_REF", "HG38_ALT",
    "SOURCE_MODEL", "GENE", "REASON", "SAMPLES",
}

# modifier TSV column -> report column; values unioned across both models
MODIFIER_FIELDS = [
    ("gene", "gene"), ("zygosity", "zygosity"), ("reason", "reason"),
    ("evidence_caller", "evidence_caller"), ("callers_detected", "callers_detected"),
    ("GT", "GT"), ("PS", "PS"),
]


def present(path):
    return bool(path) and os.path.basename(path) != "NO_FILE" and os.path.exists(path)


def read_modifier_tsvs(paths):
    """({(chrom,pos,ref,alt): aggregated info} across both models' TSVs,
    {model: number of rows}) -- rows are per variant x gene, so the row
    count can exceed that model's unique-variant count."""
    variants = {}
    n_rows = {}
    for model, path in paths:
        n_rows[model] = 0
        with open(path, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                n_rows[model] += 1
                key = (row["#CHROM"], row["POS"], row["REF"], row["ALT"])
                v = variants.setdefault(key, {
                    "models": set(), "samples": set(), "gnomad_AF_NFE": "", "kg_AF_EUR": "",
                    **{out: set() for _, out in MODIFIER_FIELDS},
                })
                v["models"].add(model)
                for src, out in MODIFIER_FIELDS:
                    for val in (row.get(src) or "").split(","):
                        if val:
                            v[out].add(val)
                v["samples"].update(s for s in (row.get("samples") or "").split(",") if s)
                for af_col in ("gnomad_AF_NFE", "kg_AF_EUR"):
                    if not v[af_col] and row.get(af_col):
                        v[af_col] = row[af_col]
    return variants, n_rows


def load_phenotypes(family_json):
    with open(family_json) as f:
        return {child: info.get("phenotype", "") for child, info in json.load(f).items()}


def scan_sample_vcfs(manifest, wanted):
    """One pass per sample VCF. Returns (carriers, ann):
    carriers: key -> [(sample, zygosity)], ann: key -> parsed ANN list
    (taken from the first VCF the variant appears in -- snpEff annotation
    depends only on the site, not the sample)."""
    carriers = defaultdict(list)
    ann = {}
    for sample, vcf in sorted(manifest.items()):
        caller_col_idx = {}
        n_hit = 0
        with opener(vcf) as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    cols = line.rstrip("\n").split("\t")
                    caller_col_idx = {name: i for i, name in enumerate(cols) if i >= 9}
                    continue
                # cheap key check before splitting the whole line
                head = line.split("\t", 5)
                if len(head) < 6:
                    continue
                key = (head[0], head[1], head[3], head[4])
                if key not in wanted:
                    continue

                fields = line.rstrip("\n").split("\t")
                info, fmt = fields[7], fields[8]
                if key not in ann:
                    ann[key] = get_ann(info) or []

                callers_field = get_info_field(info, "CALLERS")
                callers_present = callers_field.split(",") if callers_field else []
                chosen = next((c for c in CALLER_PRIORITY if c in callers_present and c in caller_col_idx), None)
                if chosen is None:
                    continue
                d = dict(zip(fmt.split(":"), fields[caller_col_idx[chosen]].split(":")))
                zygosity = classify_gt(d.get("GT", ""))
                if zygosity:
                    carriers[key].append((sample, zygosity))
                    n_hit += 1
        logger.info(f"{sample}: {n_hit} candidate variant(s) called in {os.path.basename(vcf)}")
    return carriers, ann


def best_snpeff(ann_entries, genes):
    """Most severe ANN entry, preferring the candidate gene(s) the modifier
    search qualified the variant through."""
    entries = [a for a in ann_entries if len(a) > 3]
    in_gene = [a for a in entries if a[3] in genes]
    pool = in_gene or entries
    if not pool:
        return {}
    a = max(pool, key=lambda e: IMPACT_RANK.get(e[2] if len(e) > 2 else "", -1))
    get = lambda i: a[i] if len(a) > i else ""
    return {
        "snpeff_effect": get(1), "snpeff_impact": get(2), "snpeff_gene": get(3),
        "snpeff_transcript": get(6), "snpeff_hgvs_c": get(9), "snpeff_hgvs_p": get(10),
    }


def read_hg38_vcf(path):
    """VID -> (chrom, pos, ref, alt, flipped) from the lifted (or
    already-hg38) candidate VCF; ID is the original chrom_pos_ref_alt (see
    build_alphagenome_input.py), preserved by CrossMap. flipped = the
    ALLELE_FLIPPED INFO flag from flip_ref_alt_alleles.py."""
    lifted = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8:
                lifted[fields[2]] = (fields[0], fields[1], fields[3], fields[4],
                                     "ALLELE_FLIPPED" in fields[7].split(";"))
    return lifted


# run_alphamissense.py column -> report column
ALPHAMISSENSE_FIELDS = [
    ("status", "alphamissense_status"), ("am_pathogenicity", "alphamissense_pathogenicity"),
    ("am_class", "alphamissense_class"), ("am_protein_variant", "alphamissense_protein_variant"),
    ("am_transcript", "alphamissense_transcript"), ("am_uniprot", "alphamissense_uniprot"),
    ("am_gene_mean_pathogenicity", "alphamissense_gene_mean_pathogenicity"),
    ("am_n_hits", "alphamissense_n_hits"), ("am_all_hits", "alphamissense_all_hits"),
]

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Sec": "U", "Pyl": "O", "Ter": "*",
}


def read_alphamissense(path):
    """VID -> row from a run_alphamissense.py TSV."""
    with open(path, newline="") as f:
        return {row["VID"]: row for row in csv.DictReader(f, delimiter="\t")}


def hgvs_p_short(hgvs_p):
    """snpEff HGVS.p missense (p.Val600Glu) -> AlphaMissense style (V600E);
    "" for anything else (synonymous, frameshift, ...)."""
    p = hgvs_p[2:] if hgvs_p.startswith("p.") else hgvs_p
    if len(p) < 7 or not p[3:-3].isdigit():
        return ""
    ref, alt = AA3TO1.get(p[:3]), AA3TO1.get(p[-3:])
    return f"{ref}{p[3:-3]}{alt}" if ref and alt else ""


def snpeff_agrees(ann_entries, am_all_hits):
    """yes = one of AlphaMissense's protein changes equals a snpEff missense
    HGVS.p for the variant (same residue, position and substitution -- a
    consistency check on the liftover and gene); no = snpEff has missense
    calls but none match (usually a different isoform's numbering); "" =
    snpEff calls no missense change here."""
    snpeff = {hgvs_p_short(a[10]) for a in ann_entries if len(a) > 10 and "missense" in a[1]} - {""}
    if not snpeff:
        return ""
    am = {h.split(":", 1)[1].split("=", 1)[0] for h in am_all_hits.split("|") if ":" in h}
    return "yes" if snpeff & am else "no"


def read_vcf_ids(path):
    ids = set()
    with opener(path) as f:
        for line in f:
            if not line.startswith("#"):
                ids.add(line.split("\t", 3)[2])
    return ids


def count_unmap_reasons(path):
    """CrossMap .unmap (plus flip_ref_alt_alleles.py's appended records):
    the failure reason is the last column, e.g. Fail(Unmap)."""
    reasons = defaultdict(int)
    with opener(path) as f:
        for line in f:
            if not line.startswith("#") and line.strip():
                reasons[line.rstrip("\n").split("\t")[-1]] += 1
    return reasons


def build_funnel(variants, n_rows, candidate_ids, lifted, unmap_reasons, ag_tables, am_rows, n_report):
    """[(step, metric, count, description)] -- see module docstring."""
    rows = []
    add = lambda step, metric, n, desc: rows.append((step, metric, n, desc))

    by_model = defaultdict(int)
    for v in variants.values():
        for m in v["models"]:
            by_model[m] += 1
    for model in n_rows:
        add("1_modifier_search", f"{model}_rows", n_rows[model], f"rows in cohort.{model}_modifier.tsv (one per variant x gene)")
        add("1_modifier_search", f"{model}_unique_variants", by_model[model], f"unique variants passing the {model} model")
    add("1_modifier_search", "unique_variants", len(variants), "START: unique variants across both models (one report row each)")

    n_sv = sum(1 for k in variants if k[3].startswith("<"))
    if candidate_ids is not None:
        n_af = sum(1 for k in variants if not k[3].startswith("<") and "_".join(k) not in candidate_ids)
        add("2_alphagenome_input", "removed_symbolic_sv", n_sv, "symbolic-ALT SVs, not liftable/scoreable")
        add("2_alphagenome_input", "removed_common_af", n_af,
            "gnomAD AF_nfe and 1000 Genomes EUR AF both >= alphagenome_max_af where annotated (unannotated kept)")
        add("2_alphagenome_input", "submitted", len(candidate_ids), "sent on to liftover / AlphaGenome")

    if lifted:
        n_flipped = sum(1 for hg in lifted.values() if hg[4])
        add("3_liftover", "lifted_as_is", len(lifted) - n_flipped, "lifted to hg38 unchanged (or already hg38)")
        add("3_liftover", "lifted_allele_flipped", n_flipped, "REF==ALT sites recovered by flipping (hg38 ref = candidate allele)")
        for reason, n in sorted(unmap_reasons.items()):
            add("3_liftover", f"removed_{reason}", n, "failed liftover (CrossMap .unmap)")
        add("3_liftover", "lifted_total", len(lifted), "candidates with hg38 coordinates")
        if candidate_ids is not None and len(lifted) + sum(unmap_reasons.values()) != len(candidate_ids):
            logger.warning(f"Funnel: liftover accounts for {len(lifted) + sum(unmap_reasons.values())} variant(s) "
                           f"but {len(candidate_ids)} were submitted -- check the liftover step")

    for source, ag_rows in ag_tables:
        status = defaultdict(int)
        n_partial = 0
        for r in ag_rows.values():
            status[r.get("status", "")] += 1
            if r.get("status") == "scored" and "partial;" in (r.get("note") or ""):
                n_partial += 1
        step = f"4_alphagenome_{source}"
        add(step, "input_rows", len(ag_rows), "candidates in this AlphaGenome run (one row per submitted variant)")
        add(step, "scored", status["scored"], "scored by AlphaGenome")
        add(step, "scored_partial", n_partial, "  of which >=1 scorer failed (see note column)")
        add(step, "scoring_failed", status["scoring_failed"], "every scorer failed")
        add(step, "liftover_failed", status["liftover_failed"], "no hg38 coordinates, not scored")

    if am_rows:
        status = defaultdict(int)
        am_class = defaultdict(int)
        for r in am_rows.values():
            status[r.get("status", "")] += 1
            if r.get("status") == "scored":
                am_class[r.get("am_class", "")] += 1
        step = "4_alphamissense"
        add(step, "input_rows", len(am_rows), "candidates looked up (same SV/AF pre-filter as AlphaGenome)")
        for s, desc in (("not_snv", "indel/MNV/SV, AlphaMissense scores SNVs only"),
                        ("liftover_failed", "no hg38 coordinates"),
                        ("hg38_ref_is_candidate", "hg38 reference base is the candidate allele, nothing to score"),
                        ("not_missense", "lifted SNV, not a missense change in an AlphaMissense transcript"),
                        ("scored", "missense SNVs with an AlphaMissense score")):
            add(step, s, status[s], desc)
        for c in ("likely_pathogenic", "ambiguous", "likely_benign"):
            add(step, f"scored_{c}", am_class[c], f"  of which am_class = {c}")
    add("5_report", "report_rows", n_report, "END: rows in cohort.variant_report.tsv (should equal unique_variants)")
    return rows


def read_alphagenome(path):
    """(VID -> row, result columns) from a run_alphagenome.py TSV."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        cols = [c for c in reader.fieldnames if c not in ALPHAGENOME_ID_COLUMNS]
        return {row["VID"]: row for row in reader}, cols


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-recessive_tsv", required=True)
    p.add_argument("-dominant_tsv", required=True)
    p.add_argument("-family_json", required=True)
    p.add_argument("-vcf_manifest", required=True, help="sampleID<TAB>snpeff_vcf, VCFs staged alongside")
    p.add_argument("-hg38_vcf", default="", help="Lifted (or already-hg38) candidate VCF, ID = original VID")
    p.add_argument("-alphamissense_tsv", default="", help="run_alphamissense.py output")
    p.add_argument("-alphagenome_tsv", default="", help="run_alphagenome.py output, local backend")
    p.add_argument("-alphagenome_api_tsv", default="", help="run_alphagenome.py output, API backend")
    p.add_argument("-candidate_vcf", default="", help="AlphaGenome input VCF (build_alphagenome_input.py), for the funnel")
    p.add_argument("-unmap", default="", help="Liftover .unmap file, for the funnel")
    p.add_argument("-output", required=True)
    p.add_argument("-funnel_output", required=True, help="Per-step variant counts TSV")
    args = p.parse_args()

    variants, n_rows = read_modifier_tsvs([("recessive", args.recessive_tsv), ("dominant", args.dominant_tsv)])
    logger.info(f"{len(variants)} unique candidate variant(s) across recessive + dominant")

    phenotypes = load_phenotypes(args.family_json)
    carriers, ann = scan_sample_vcfs(load_vcf_manifest(args.vcf_manifest), set(variants))

    lifted = read_hg38_vcf(args.hg38_vcf) if present(args.hg38_vcf) else {}
    vid_of = {key: "_".join(key) for key in variants}

    am_rows = read_alphamissense(args.alphamissense_tsv) if present(args.alphamissense_tsv) else {}

    ag_tables = []
    ag_cols = []
    for source, path in (("local", args.alphagenome_tsv), ("api", args.alphagenome_api_tsv)):
        if present(path):
            rows, cols = read_alphagenome(path)
            ag_tables.append((source, rows))
            ag_cols = ag_cols or cols

    columns = [
        "VID", "CHROM", "POS", "REF", "ALT", "HG38_CHROM", "HG38_POS", "HG38_REF", "HG38_ALT", "allele_flipped",
        "source_model", *[out for _, out in MODIFIER_FIELDS],
        "modifier_samples", "N_modifier_samples",
        "N_mild_carriers", "mild_carriers", "N_severe_carriers", "severe_carriers", "other_carriers",
        "gnomad_AF_NFE", "kg_AF_EUR",
        "snpeff_effect", "snpeff_impact", "snpeff_gene", "snpeff_transcript", "snpeff_hgvs_c", "snpeff_hgvs_p",
        "alphamissense_status", "alphamissense_pathogenicity", "alphamissense_class", "alphamissense_protein_variant",
        "alphamissense_snpeff_agree", "alphamissense_transcript", "alphamissense_uniprot",
        "alphamissense_gene_mean_pathogenicity", "alphamissense_n_hits", "alphamissense_all_hits",
        "alphagenome_source", *[f"alphagenome_{c}" for c in ag_cols],
    ]

    def sort_key(key):
        return (key[0], int(key[1]) if key[1].isdigit() else 0, key[2], key[3])

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for key in sorted(variants, key=sort_key):
            v = variants[key]
            vid = vid_of[key]
            row = {
                "VID": vid, "CHROM": key[0], "POS": key[1], "REF": key[2], "ALT": key[3],
                "source_model": ",".join(sorted(v["models"])),
                **{out: ",".join(sorted(v[out])) for _, out in MODIFIER_FIELDS},
                "modifier_samples": ",".join(sorted(v["samples"])),
                "N_modifier_samples": len(v["samples"]),
                "gnomad_AF_NFE": v["gnomad_AF_NFE"], "kg_AF_EUR": v["kg_AF_EUR"],
            }

            by_pheno = defaultdict(list)
            for sample, zyg in carriers.get(key, []):
                by_pheno[phenotypes.get(sample, "")].append(f"{sample}({zyg})")
            row.update({
                "N_mild_carriers": len(by_pheno["mild"]), "mild_carriers": ",".join(by_pheno["mild"]),
                "N_severe_carriers": len(by_pheno["severe"]), "severe_carriers": ",".join(by_pheno["severe"]),
                "other_carriers": ",".join(c for ph, cs in by_pheno.items() if ph not in ("mild", "severe") for c in cs),
            })

            row.update(best_snpeff(ann.get(key, []), v["gene"]))

            if vid in lifted:
                hg = lifted[vid]
                row.update({"HG38_CHROM": hg[0], "HG38_POS": hg[1], "HG38_REF": hg[2], "HG38_ALT": hg[3],
                            "allele_flipped": hg[4]})

            am = am_rows.get(vid)
            if am:
                row.update({out: am.get(src, "") for src, out in ALPHAMISSENSE_FIELDS})
                if am.get("status") == "scored":
                    row["alphamissense_snpeff_agree"] = snpeff_agrees(ann.get(key, []), am.get("am_all_hits", ""))
                # AlphaMissense lifts every candidate (AlphaGenome skips common ones)
                if not row.get("HG38_CHROM") and am.get("HG38_CHROM"):
                    row.update({c: am[c] for c in ("HG38_CHROM", "HG38_POS", "HG38_REF", "HG38_ALT")})
            elif am_rows:
                row["alphamissense_status"] = "not_submitted"

            # local result wins when it scored; otherwise the first table
            # that did, otherwise whichever has a row (to surface its status)
            ag_rows = [(src, rows[vid]) for src, rows in ag_tables if vid in rows]
            pick = next(((s, r) for s, r in ag_rows if r.get("status") == "scored"), ag_rows[0] if ag_rows else None)
            if pick:
                row["alphagenome_source"] = pick[0]
                row.update({f"alphagenome_{c}": pick[1].get(c, "") for c in ag_cols})
            elif ag_tables:
                # never sent to AlphaGenome: SV, or dropped by the gnomAD AF
                # pre-filter (build_alphagenome_input.py -max_af)
                row["alphagenome_status"] = "not_submitted"

            writer.writerow(row)

    logger.info(f"Wrote {len(variants)} variant(s) -> {args.output}")

    funnel = build_funnel(
        variants, n_rows,
        read_vcf_ids(args.candidate_vcf) if present(args.candidate_vcf) else None,
        lifted,
        count_unmap_reasons(args.unmap) if present(args.unmap) else {},
        ag_tables, am_rows, len(variants),
    )
    with open(args.funnel_output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["step", "metric", "count", "description"])
        writer.writerows(funnel)
    logger.info("Variant funnel:\n" + "\n".join(f"  {step:<22} {metric:<32} {n:>8}" for step, metric, n, _ in funnel))


if __name__ == "__main__":
    main()
