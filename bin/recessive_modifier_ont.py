#!/usr/bin/env python3

import pandas as pd
import numpy as np
import time
import os
import logging
import colorlog
import itertools
import multiprocessing as mp
import psutil
import argparse
import gc
from collections import defaultdict
import datetime

timestamp = datetime.datetime.now().strftime("%y%m%d_%H%M%S")

# Configure logging
level = "INFO"
logger = logging.getLogger("color_logger")

if not logger.handlers:
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        "%(asctime)s - %(log_color)s%(levelname)s:%(reset)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        log_colors={
            'DEBUG': 'blue',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'bold_red',
        }
    ))
    logger.addHandler(handler)
    logger.setLevel(getattr(logging, level, logging.INFO))


def timing(f):
    """Decorator for timing functions."""
    def wrap(*args, **kwargs):
        time1 = time.time()
        ret = f(*args, **kwargs)
        time2 = time.time()
        logger.info(f"{f.__name__} function took {time2 - time1:.3f} s")
        return ret
    return wrap


def get_available_cpus():
    """Get a list of available CPU cores."""
    try:
        return list(os.sched_getaffinity(0))  # Linux
    except AttributeError:
        return list(range(os.cpu_count()))


def get_memory_usage():
    """Get current memory usage in GB."""
    process = psutil.Process(os.getpid())
    mem_gb = process.memory_info().rss / 1024 / 1024 / 1024
    return mem_gb


def get_ann(info_str):
    """Extract and parse the ANN field from VCF INFO string."""
    for field in info_str.split(';'):
        if field.startswith('ANN='):
            ann = field.split('=', 1)[1]
            return [x.split("|") for x in ann.split(",")]
    return False


@timing
def parse_vcf(vcf, model):
    """
    Parse a VCF file into variant structures.

    Args:
        vcf: Path to VCF file
        model: 'dominant' or 'recessive'

    Returns:
        dict: Variant dictionary with structure:
              {var_id: (patient_list, zygosity, ann_list, chrom, pos, ref, alt)}
    """
    results = {}
    het_variants = {}
    het_gene_variant_map = defaultdict(set)
    patient = os.path.basename(vcf).split(".")[0]

    # Define allowed annotation types
    allowed_annotations = {
        'intron_variant',
        'intron',
        '5_prime_UTR_variant',
        'five_prime_UTR_variant',
        '3_prime_UTR_variant',
        'three_prime_UTR_variant',
        'exon_variant',
        'exon',
        'synonymous_variant',
        'missense_variant',
        'stop_gained',
        'stop_lost',
        'start_lost',
        'splice_donor_variant',
        'splice_acceptor_variant',
        'splice_region_variant',
        'frameshift_variant',
        'inframe_insertion',
        'inframe_deletion',
        'protein_altering_variant'
    }

    variant_count = 0
    filtered_count = 0

    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue

            info = line.rstrip("\n").split("\t")
            if "ANN" not in info[7]:
                continue

            var_id = info[2]
            zygosity = var_id.split("|")[3]
            chrom = info[0]
            pos = info[1]
            ref = info[3]
            alt = info[4]

            ann = get_ann(info[7])
            if not ann:
                continue

            # Filter annotations to only keep allowed types
            # ANN format: Allele|Annotation|Impact|Gene|... (index 1 is annotation type)
            filtered_ann = []
            for a in ann:
                if len(a) > 1:
                    ann_type = a[1]
                    if any(allowed in ann_type for allowed in allowed_annotations):
                        filtered_ann.append(a)

            # Skip variant if no annotations match our criteria
            if not filtered_ann:
                filtered_count += 1
                continue

            variant_count += 1

            # Store variant info
            if model == "dominant":
                results[var_id] = ([patient], zygosity, filtered_ann, chrom, pos, ref, alt)

            elif model == "recessive":
                if zygosity == "homo":
                    results[var_id] = ([patient], zygosity, filtered_ann, chrom, pos, ref, alt)
                elif zygosity == "het":
                    het_variants[var_id] = ([patient], zygosity, filtered_ann, chrom, pos, ref, alt)
                    genes = {a[3] for a in filtered_ann if len(a) > 3 and a[3]}
                    for gene in genes:
                        het_gene_variant_map[gene].add(var_id)

    # For recessive model, add het variants that could be compound hets:
    # keep a variant's FULL annotation list whenever any of its genes has
    # >=2 distinct het variants (not just this one annotation entry).
    if model == "recessive":
        for gene, var_id_set in het_gene_variant_map.items():
            if len(var_id_set) >= 2:
                for var_id in var_id_set:
                    results[var_id] = het_variants[var_id]

    logger.debug(f"Parsed {os.path.basename(vcf)}: kept {variant_count} variants, "
                f"filtered {filtered_count} (wrong annotation type)")
    return results


def process_vcf(vcf, model):
    """Process a VCF file or list of VCF files."""
    if isinstance(vcf, list):
        merged = {}
        for v in vcf:
            r = parse_vcf(v, model)
            for var_id, (patients, zyg, ann, chrom, pos, ref, alt) in r.items():
                if var_id in merged:
                    merged[var_id][0].extend(patients)
                else:
                    merged[var_id] = [patients, zyg, ann, chrom, pos, ref, alt]
        return merged
    else:
        return parse_vcf(vcf, model)


def process_wrapper(args):
    """Wrapper that sets CPU affinity before running the target function."""
    func, vcf_args, assigned_cpus = args
    pid = os.getpid()
    try:
        psutil.Process(pid).cpu_affinity(assigned_cpus)
    except AttributeError:
        pass
    logger.debug(f"Processing on PID {pid} using CPUs {assigned_cpus}")
    return func(*vcf_args)


@timing
def parallel_read_vcf_concordant(vcfs, model):
    """Read concordant severe VCF files in parallel."""
    if not vcfs:
        return {}

    available_cpus = get_available_cpus()
    num_parallel_tasks = min(len(vcfs), 6)
    cpus_per_task = max(1, len(available_cpus) // num_parallel_tasks)

    logger.info(f"Reading {len(vcfs)} concordant VCFs: {num_parallel_tasks} parallel tasks")
    logger.info(f"Memory before parallel processing: {get_memory_usage():.2f} GB")

    assigned_cpu_groups = itertools.cycle(
        [available_cpus[i:i + cpus_per_task] for i in range(0, len(available_cpus), cpus_per_task)]
    )

    task_args = [(process_vcf, [v, model], next(assigned_cpu_groups)) for v in vcfs]

    with mp.Pool(processes=num_parallel_tasks) as pool:
        results = pool.map(process_wrapper, task_args)

    # Merge all results
    merged = {}
    for d in results:
        for var_id, (patients, zyg, ann, chrom, pos, ref, alt) in d.items():
            if var_id in merged:
                merged[var_id][0].extend(patients)
            else:
                merged[var_id] = [patients, zyg, ann, chrom, pos, ref, alt]

    logger.info(f"Loaded {len(merged)} unique variants from concordant VCFs")
    logger.info(f"Memory after parallel processing: {get_memory_usage():.2f} GB")

    gc.collect()
    return merged


@timing
def read_vcfs_discordant(mild, severe, model):
    """Read discordant sibling pair VCFs."""
    logger.info(f"Reading discordant pair: {os.path.basename(mild)} (mild) vs "
                f"{[os.path.basename(s) for s in severe]} (severe)")

    mild_result = process_vcf(mild, model)
    severe_result = process_vcf(severe, model)

    return mild_result, severe_result


def dominant_model_from_dict(mild, sev):
    """
    Apply dominant model filtering.

    Filters out:
    1. Variants in mild that also appear in severe (same variant)
    2. ANY variant in mild for genes with variants in severe
    """
    result = {}

    # Build set of exact variants in severe
    sev_keys = set(sev.keys())

    # Build set of genes with ANY variants in severe
    sev_genes = set()
    for _, (__, ___, annlist, ____, _____, ______, _______) in sev.items():
        for ann in annlist:
            sev_genes.add(ann[3])  # ann[3] is gene name

    logger.info(f"Severe has {len(sev_keys)} variants across {len(sev_genes)} genes")

    # Filter mild variants
    filtered_same_variant = 0
    filtered_gene_in_severe = 0

    for var, (patients, zygosity, annlist, chrom, pos, ref, alt) in mild.items():
        # Filter if exact variant is in severe
        if var in sev_keys:
            filtered_same_variant += 1
            continue

        # Get genes for this variant
        genes = {ann[3] for ann in annlist}

        # Filter if this variant is in a gene with any variant in severe
        if genes & sev_genes:
            filtered_gene_in_severe += 1
            continue

        # Keep this variant
        result[var] = (patients, zygosity, annlist, chrom, pos, ref, alt)

    logger.info(f"Dominant model filtering:")
    logger.info(f"  - Filtered {filtered_same_variant} variants (same variant as severe)")
    logger.info(f"  - Filtered {filtered_gene_in_severe} variants (gene has variant in severe)")
    logger.info(f"  - Total filtered: {len(mild)-len(result)}/{len(mild)}")
    logger.info(f"  - Remaining: {len(result)} variants")

    return result


def recessive_model_from_dict(mild, sev):
    """
    Apply recessive model filtering.

    Filters out:
    1. Homozygous variants in mild that also appear in severe (same variant)
    2. ANY variant (homo or het) in mild for genes with homozygous variants in severe
    3. Heterozygous variants in mild for genes with het variants in severe
    """
    result = {}

    # Build sets for efficient lookup
    sev_homo_set = {k for k, (_, z, _, _, _, _, _) in sev.items() if z == "homo"}

    # Build set of genes with homozygous variants in severe
    sev_genes_with_homo = set()
    for _, (__, z, annlist, ___, ____, _____, ______) in sev.items():
        if z == "homo":
            for ann in annlist:
                sev_genes_with_homo.add(ann[3])  # ann[3] is gene name

    # Build set of genes with heterozygous variants in severe
    sev_genes_with_het = set()
    for _, (__, z, annlist, ___, ____, _____, ______) in sev.items():
        if z == "het":
            for ann in annlist:
                sev_genes_with_het.add(ann[3])

    logger.info(f"Severe has {len(sev_genes_with_homo)} genes with homozygous variants")
    logger.info(f"Severe has {len(sev_genes_with_het)} genes with heterozygous variants")

    # Filter mild variants
    filtered_homo_same_variant = 0
    filtered_gene_with_severe_homo = 0
    filtered_het_in_severe_gene = 0

    for var, (patients, zygosity, annlist, chrom, pos, ref, alt) in mild.items():
        # Get genes for this variant
        genes = {ann[3] for ann in annlist}

        if zygosity == "homo":
            # Filter if exact variant is in severe homo
            if var in sev_homo_set:
                filtered_homo_same_variant += 1
                continue

            # Filter if this variant is in a gene with homozygous variant in severe
            if genes & sev_genes_with_homo:
                filtered_gene_with_severe_homo += 1
                continue

            # Keep this homo variant
            result[var] = (patients, zygosity, annlist, chrom, pos, ref, alt)

        elif zygosity == "het":
            # Filter if exact variant is in severe homo
            if var in sev_homo_set:
                filtered_homo_same_variant += 1
                continue

            # Filter if this variant is in a gene with homozygous variant in severe
            if genes & sev_genes_with_homo:
                filtered_gene_with_severe_homo += 1
                continue

            # Filter if this variant is in a gene with het variants in severe
            if genes & sev_genes_with_het:
                filtered_het_in_severe_gene += 1
                continue

            # Keep this het variant
            result[var] = (patients, zygosity, annlist, chrom, pos, ref, alt)

    logger.info(f"Recessive model filtering:")
    logger.info(f"  - Filtered {filtered_homo_same_variant} variants (same variant as severe homo)")
    logger.info(f"  - Filtered {filtered_gene_with_severe_homo} variants (gene has homo variant in severe)")
    logger.info(f"  - Filtered {filtered_het_in_severe_gene} het variants (gene has het variant in severe)")
    logger.info(f"  - Total filtered: {len(mild)-len(result)}/{len(mild)}")
    logger.info(f"  - Remaining: {len(result)} variants")

    return result


def format_ann(ann_list):
    """Format annotation list back to VCF ANN format."""
    if not ann_list or ann_list == "unknown":
        return "ANN=unknown"

    formatted_anns = []
    for ann in ann_list:
        if isinstance(ann, list):
            formatted_anns.append('|'.join(str(x) for x in ann))
        else:
            formatted_anns.append(str(ann))

    return "ANN=" + ','.join(formatted_anns)


def extract_genes(ann_list):
    """
    Extract unique gene names from annotation list.

    Args:
        ann_list: List of annotation arrays from SnpEff

    Returns:
        str: Comma-separated unique gene names, or "unknown" if none found
    """
    if not ann_list or ann_list == "unknown":
        return "unknown"

    genes = set()
    for ann in ann_list:
        if isinstance(ann, list) and len(ann) > 3:
            gene = ann[3]  # Gene name is at index 3 in SnpEff ANN
            if gene and gene != '':
                genes.add(gene)

    if not genes:
        return "unknown"

    return ','.join(sorted(genes))


@timing
def parallel_recessive_modifier_ont(discordant_sib_pairs, concordant_severe_sibs, model):
    """
    Main analysis function (no AlphaMissense / gnomAD annotation).

    Flow:
    1. Parse VCFs
    2. Apply genetic model filtering
    3. Apply multi-sib filtering
    4. Build final results table
    """
    logger.info(f"\n{'='*60}")
    logger.info(f"Starting analysis with {model} model")
    logger.info(f"Discordant pairs: {len(discordant_sib_pairs)}")
    logger.info(f"Concordant severe siblings: {len(concordant_severe_sibs)}")
    logger.info(f"Current memory usage: {get_memory_usage():.2f} GB")
    logger.info(f"{'='*60}\n")

    # Step 1: Load concordant severe variants
    logger.info("PHASE 1: Reading concordant severe VCFs...")
    logger.info("Filtering for: intron, 5'UTR, 3'UTR, exon variants only")
    concord = parallel_read_vcf_concordant(concordant_severe_sibs, model)

    # Initialize tracking structures
    res = {}
    var_to_genes_map = defaultdict(set)
    variant_sib_map = defaultdict(set)

    # Step 2: Process each discordant pair
    logger.info("\nPHASE 2: Processing discordant pairs and applying genetic model...")
    for pair_idx, discord in enumerate(discordant_sib_pairs, 1):
        mild_file = discord[0]
        sev_files = discord[1]

        logger.info(f"--- Pair {pair_idx}/{len(discordant_sib_pairs)}: {os.path.basename(mild_file)} ---")

        mild, sev = read_vcfs_discordant(mild_file, sev_files, model)

        sev_plus_concord = concord.copy()
        sev_plus_concord.update(sev)

        if model == "dominant":
            passing_dict = dominant_model_from_dict(mild, sev_plus_concord)
        else:
            passing_dict = recessive_model_from_dict(mild, sev_plus_concord)

        sample_id = os.path.basename(mild_file).split(".")[0]

        for var, (patients, zygosity, annlist, chrom, pos, ref, alt) in passing_dict.items():
            variant_sib_map[var].add(sample_id)

            if var not in res:
                res[var] = [patients, zygosity, annlist, chrom, pos, ref, alt]
            else:
                res[var][0].extend(patients)

            if model == "recessive":
                for ann in annlist:
                    gene = ann[3]
                    if gene:
                        var_to_genes_map[var].add(gene)

        gc.collect()

    logger.info(f"\n{'='*60}")
    logger.info(f"After genetic model filtering: {len(res)} variants")
    logger.info(f"Memory usage: {get_memory_usage():.2f} GB")
    logger.info(f"{'='*60}\n")

    # Step 3: Final filtering based on sib counts
    logger.info("PHASE 3: Applying multi-sib filtering...")
    final_pass = {}

    if model == "dominant":
        for var, (patients, zygosity, annlist, chrom, pos, ref, alt) in res.items():
            if len(variant_sib_map[var]) >= 2:
                final_pass[var] = [patients, zygosity, annlist, chrom, pos, ref, alt]

        logger.info(f"Dominant final filter: {len(final_pass)} variants in >=2 sibs")

    else:  # recessive
        for var, (patients, zygosity, annlist, chrom, pos, ref, alt) in res.items():
            if zygosity == "homo" and len(variant_sib_map[var]) >= 2:
                final_pass[var] = [patients, zygosity, annlist, chrom, pos, ref, alt]

        logger.info(f"Recessive homo filter: {len(final_pass)} homo variants in >=2 sibs")

        gene_sib_counts = defaultdict(lambda: defaultdict(int))

        for var, (patients, zygosity, annlist, chrom, pos, ref, alt) in res.items():
            if zygosity != "het":
                continue

            genes = var_to_genes_map.get(var)
            if not genes:
                continue

            for gene in genes:
                for sib in variant_sib_map[var]:
                    gene_sib_counts[gene][sib] += 1

        genes_to_keep = {
            gene
            for gene, sib_count_map in gene_sib_counts.items()
            if sum(1 for count in sib_count_map.values() if count >= 2) >= 2
        }

        logger.info(f"Recessive het filter: {len(genes_to_keep)} genes with compound het pattern")

        het_count = 0
        for var, (patients, zygosity, annlist, chrom, pos, ref, alt) in res.items():
            if zygosity == "het" and var_to_genes_map.get(var, set()) & genes_to_keep:
                final_pass[var] = [patients, zygosity, annlist, chrom, pos, ref, alt]
                het_count += 1

        logger.info(f"Recessive het filter: {het_count} het variants in qualifying genes")

    # Clear intermediate data
    del res, var_to_genes_map, variant_sib_map, concord
    gc.collect()

    logger.info(f"\n{'='*60}")
    logger.info(f"Variants passing all filters: {len(final_pass)}")
    logger.info(f"{'='*60}\n")

    if not final_pass:
        logger.warning("No variants passed filtering - returning empty results table")
        return pd.DataFrame(columns=[
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "zygosity", "gene", "patients", "N_sibs"
        ])

    # Step 4: Build final DataFrame
    logger.info(f"\n{'='*60}")
    logger.info(f"PHASE 4: Building final results table...")
    logger.info(f"{'='*60}\n")

    final_df = pd.DataFrame.from_dict(final_pass, orient='index').reset_index()
    final_df = final_df.rename(columns={
        'index': 'variant_id',
        0: 'patients',
        1: 'zygosity',
        2: 'annotation',
        3: 'chrom',
        4: 'pos',
        5: 'ref',
        6: 'alt',
    })

    final_df['#CHROM'] = final_df['chrom']
    final_df['POS'] = final_df['pos']
    final_df['REF'] = final_df['ref']
    final_df['ALT'] = final_df['alt']

    final_df['INFO'] = final_df['annotation'].apply(
        lambda x: format_ann(x if x else "unknown")
    )

    # Extract gene names
    final_df['gene'] = final_df['annotation'].apply(extract_genes)

    final_df['QUAL'] = 100
    final_df['FILTER'] = "PASS"
    final_df['ID'] = final_df['variant_id']
    final_df['N_sibs'] = final_df['patients'].apply(lambda x: len(set(x)))

    final_df = final_df[[
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "zygosity", "gene", "patients", "N_sibs"
    ]].reset_index(drop=True)

    # Summary statistics
    N_homo = final_df.query('zygosity == "homo"').shape[0]
    N_het = final_df.query('zygosity == "het"').shape[0]

    # Get unique genes
    all_genes = set()
    for genes_str in final_df['gene']:
        if genes_str != 'unknown':
            all_genes.update(genes_str.split(','))

    logger.info(f"\n{'='*60}")
    logger.info(f"FINAL RESULTS - {model.upper()} MODEL")
    logger.info(f"{'='*60}")
    logger.info(f"Total variants: {N_homo + N_het}")
    logger.info(f"  - Homozygous: {N_homo}")
    logger.info(f"  - Heterozygous: {N_het}")
    logger.info(f"  - Unique genes: {len(all_genes)}")
    logger.info(f"Final memory usage: {get_memory_usage():.2f} GB")
    logger.info(f"{'='*60}\n")

    return final_df


####################################
def get_args():
    """Parse command line parameters."""
    parser = argparse.ArgumentParser(
        description='Recessive Modifier Analysis (ONT long-read, no external allele-frequency/pathogenicity annotation)'
    )
    parser.add_argument("-model", type=str, required=True, choices=['dominant', 'recessive'])
    parser.add_argument("-discordant", type=str, required=True,
                        help="Path to a tab-delimited file listing discordant pairs: "
                             "mild_vcf<TAB>severe_vcf1,severe_vcf2,... (one pair per line)")
    parser.add_argument("-concordant", type=str, default=None,
                        help="Path to a file listing one concordant-severe VCF per line")
    parser.add_argument("-output", type=str, default=None)

    return parser.parse_args()


def load_discordant_pairs(path):
    """Load discordant sib pairs from a tab-delimited file: mild_vcf<TAB>severe_vcf1,severe_vcf2,..."""
    pairs = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(f"Malformed discordant pair line (expected 2 tab-separated fields): {line}")
            mild_vcf, severe_vcfs = fields
            pairs.append([mild_vcf, severe_vcfs.split(",")])
    return pairs


def load_concordant_sibs(path):
    """Load concordant severe VCF paths, one per line."""
    if not path:
        return []
    with open(path, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def main():
    """Main execution function."""
    args = get_args()

    logger.info(f"Starting analysis (no AlphaMissense/gnomAD annotation)")
    logger.info(f"Initial memory usage: {get_memory_usage():.2f} GB")

    discordant_sibs = load_discordant_pairs(args.discordant)
    concordant_sibs = load_concordant_sibs(args.concordant)

    # Run analysis
    res = parallel_recessive_modifier_ont(
        discordant_sibs,
        concordant_sibs,
        args.model,
    )

    # Save results
    output_file = args.output or f"{timestamp}_results_{args.model}.ont.vcf"
    res.to_csv(output_file, index=False, sep="\t", header=True)
    logger.info(f"\nResults saved to: {output_file}")
    logger.info(f"Peak memory usage: {get_memory_usage():.2f} GB")


if __name__ == "__main__":
    main()
