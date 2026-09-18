#!/usr/bin/env python3

import pandas as pd
import numpy as np
import time
import os
import pickle
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


class GnomADAnnotator:
    """Lazy-loading gnomAD annotator."""
    
    def __init__(self, gnomad_pickle_file=None):
        self.gnomad_file = gnomad_pickle_file
        self.af_dict = None
    
    @timing
    def load_gnomad_data(self):
        """Load gnomAD data only when needed."""
        if self.af_dict is not None:
            return  # Already loaded
        
        if not self.gnomad_file:
            logger.warning("No gnomAD file provided")
            return
        
        logger.info(f"Loading gnomAD allele frequency data from {self.gnomad_file}...")
        logger.info(f"Memory before loading: {get_memory_usage():.2f} GB")
        
        try:
            with open(self.gnomad_file, 'rb') as f:
                self.af_dict = pickle.load(f)
            
            logger.info(f"gnomAD data loaded successfully: {len(self.af_dict):,} variants")
            logger.info(f"Memory after loading: {get_memory_usage():.2f} GB")
        except Exception as e:
            logger.error(f"Failed to load gnomAD data: {e}")
            self.af_dict = {}
    
    def annotate_variant(self, chrom, pos, ref, alt):
        """Look up gnomAD allele frequency for a variant."""
        if self.af_dict is None:
            return None
        
        variant_id = f"{chrom}:{pos}:{ref}>{alt}"
        return self.af_dict.get(variant_id, None)
    
    def annotate_batch(self, variants):
        """
        Annotate a batch of variants efficiently.
        
        Args:
            variants: List of (chrom, pos, ref, alt) tuples
        
        Returns:
            dict: variant_id -> allele_frequency
        """
        if self.af_dict is None:
            return {f"{c}:{p}:{r}>{a}": None for c, p, r, a in variants}
        
        results = {}
        for chrom, pos, ref, alt in variants:
            variant_id = f"{chrom}:{pos}:{ref}>{alt}"
            results[variant_id] = self.af_dict.get(variant_id, None)
        
        return results


class AlphaMissenseAnnotator:
    """
    Lazy-loading AlphaMissense annotator with transcript-first lookup.
    
    Lookup strategy:
    1. First try: CHROM:POS:REF:ALT:TRANSCRIPT_ID (transcript-specific)
    2. Fallback: CHROM:POS:REF:ALT (variant-level, max pathogenicity)
    """
    
    def __init__(self, alphamissense_file=None):
        self.alphamissense_file = alphamissense_file
        self.transcript_lookup = None  # CHROM:POS:REF:ALT:TRANSCRIPT -> pathogenicity
        self.variant_lookup = None     # CHROM:POS:REF:ALT -> max pathogenicity
    
    @timing
    def load_alphamissense_data(self):
        """Load AlphaMissense data only when needed."""
        if self.transcript_lookup is not None:
            return  # Already loaded
        
        if not self.alphamissense_file:
            logger.warning("No AlphaMissense file provided")
            return
        
        logger.info("Loading AlphaMissense annotation data...")
        logger.info(f"Memory before loading: {get_memory_usage():.2f} GB")
        
        # Check for cached pickle
        pickle_file = f"{self.alphamissense_file}.pickle"
        load_from_pickle = False
        
        if os.path.exists(pickle_file):
            try:
                logger.info(f"Loading cached AlphaMissense data from pickle...")
                with open(pickle_file, 'rb') as f:
                    cached_data = pickle.load(f)
                    
                # Validate pickle format
                if isinstance(cached_data, dict) and 'transcript_lookup' in cached_data and 'variant_lookup' in cached_data:
                    self.transcript_lookup = cached_data['transcript_lookup']
                    self.variant_lookup = cached_data['variant_lookup']
                    load_from_pickle = True
                else:
                    logger.warning("Pickle file has incompatible format - will regenerate")
                    logger.warning(f"Deleting old pickle: {pickle_file}")
                    os.remove(pickle_file)
            except Exception as e:
                logger.warning(f"Failed to load pickle file: {e}")
                logger.warning(f"Will regenerate from TSV")
                if os.path.exists(pickle_file):
                    os.remove(pickle_file)
        
        if not load_from_pickle:
            logger.info(f"Parsing AlphaMissense file: {self.alphamissense_file}")
            
            # Read the file, skipping comment lines (lines starting with # except the header)
            # The header line starts with #CHROM so we need to handle it specially
            df = pd.read_csv(
                self.alphamissense_file,
                sep="\t",
                compression='gzip' if self.alphamissense_file.endswith('.gz') else None,
                comment=None,  # Don't auto-skip # lines
                skiprows=lambda x: x < 2  # Skip first 2 lines (copyright/license)
            )
            
            logger.info(f"Loaded {len(df)} rows from AlphaMissense file")
            logger.info(f"Columns found: {list(df.columns)}")
            
            # The first column will have the # in it, clean it up
            df.columns = [col.lstrip('#') for col in df.columns]
            
            logger.info(f"Cleaned columns: {list(df.columns)}")
            
            # Now get the columns we need
            required_cols = ['CHROM', 'POS', 'REF', 'ALT', 'transcript_id', 'am_pathogenicity']
            
            # Verify all required columns exist
            missing = [col for col in required_cols if col not in df.columns]
            
            if missing:
                raise ValueError(f"Missing required columns: {missing}. Available: {list(df.columns)}")
            
            logger.info("All required columns found!")
            
            # Select only the columns we need
            df = df[required_cols]
            
            # Create transcript-specific lookup: CHROM:POS:REF:ALT:TRANSCRIPT
            logger.info("Building transcript-specific lookup...")
            df['transcript_key'] = (
                df['CHROM'].astype(str) + ':' +
                df['POS'].astype(str) + ':' +
                df['REF'] + ':' +
                df['ALT'] + ':' +
                df['transcript_id']
            )
            self.transcript_lookup = df.set_index('transcript_key')['am_pathogenicity'].to_dict()
            
            # Create variant-level lookup: CHROM:POS:REF:ALT (max across transcripts)
            logger.info("Building variant-level fallback lookup...")
            df['variant_key'] = (
                df['CHROM'].astype(str) + ':' +
                df['POS'].astype(str) + ':' +
                df['REF'] + ':' +
                df['ALT']
            )
            self.variant_lookup = df.groupby('variant_key')['am_pathogenicity'].max().to_dict()
            
            # Save both lookups to pickle
            logger.info("Caching lookups to pickle for faster future loading...")
            with open(pickle_file, 'wb') as f:
                pickle.dump({
                    'transcript_lookup': self.transcript_lookup,
                    'variant_lookup': self.variant_lookup
                }, f)
            
            del df
            gc.collect()
        
        logger.info(f"AlphaMissense data loaded:")
        logger.info(f"  - Transcript-specific entries: {len(self.transcript_lookup):,}")
        logger.info(f"  - Unique variants (fallback): {len(self.variant_lookup):,}")
        logger.info(f"Memory after loading: {get_memory_usage():.2f} GB")
    
    def annotate_variant(self, chrom, pos, ref, alt, transcript_id=None):
        """
        Look up AlphaMissense pathogenicity for a variant.
        
        First tries transcript-specific lookup, then falls back to variant-level.
        """
        if self.transcript_lookup is None:
            return np.nan
        
        # Try transcript-specific lookup first
        if transcript_id:
            transcript_key = f"{chrom}:{pos}:{ref}:{alt}:{transcript_id}"
            if transcript_key in self.transcript_lookup:
                return self.transcript_lookup[transcript_key]
        
        # Fallback to variant-level lookup
        variant_key = f"{chrom}:{pos}:{ref}:{alt}"
        return self.variant_lookup.get(variant_key, np.nan)
    
    def annotate_batch(self, variants):
        """
        Annotate a batch of variants efficiently.
        
        Args:
            variants: List of (chrom, pos, ref, alt, transcript_id) tuples
        
        Returns:
            dict: variant_key -> pathogenicity_score
        """
        if self.transcript_lookup is None:
            return {f"{c}:{p}:{r}:{a}:{t}": np.nan for c, p, r, a, t in variants}
        
        results = {}
        transcript_hits = 0
        variant_hits = 0
        misses = 0
        
        for chrom, pos, ref, alt, transcript_id in variants:
            transcript_key = f"{chrom}:{pos}:{ref}:{alt}:{transcript_id}"
            
            # Try transcript-specific lookup first
            if transcript_id and transcript_key in self.transcript_lookup:
                results[transcript_key] = self.transcript_lookup[transcript_key]
                transcript_hits += 1
            else:
                # Fallback to variant-level lookup
                variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                if variant_key in self.variant_lookup:
                    results[transcript_key] = self.variant_lookup[variant_key]
                    variant_hits += 1
                else:
                    results[transcript_key] = np.nan
                    misses += 1
        
        logger.debug(f"Annotation stats: {transcript_hits} transcript hits, "
                    f"{variant_hits} variant fallback hits, {misses} misses")
        
        return results


def get_ann(info_str):
    """Extract and parse the ANN field from VCF INFO string."""
    for field in info_str.split(';'):
        if field.startswith('ANN='):
            ann = field.split('=', 1)[1]
            return [x.split("|") for x in ann.split(",")]
    return False


@timing
def parse_vcf_no_annotation(vcf, model):
    """
    Parse VCF WITHOUT annotation - just get variant structure.
    This is MUCH faster since we skip all annotation lookups.
    
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
            
            info = line.strip().split("\t")
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
                    # Check if any part of the annotation matches allowed types
                    if any(allowed in ann_type for allowed in allowed_annotations):
                        filtered_ann.append(a)
            
            # Skip variant if no annotations match our criteria
            if not filtered_ann:
                filtered_count += 1
                continue
            
            variant_count += 1
            
            # Store variant info WITHOUT annotation
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
    """Process a VCF file or list of VCF files WITHOUT annotation."""
    if isinstance(vcf, list):
        merged = {}
        for v in vcf:
            r = parse_vcf_no_annotation(v, model)
            for var_id, (patients, zyg, ann, chrom, pos, ref, alt) in r.items():
                if var_id in merged:
                    merged[var_id][0].extend(patients)
                else:
                    merged[var_id] = [patients, zyg, ann, chrom, pos, ref, alt]
        return merged
    else:
        return parse_vcf_no_annotation(vcf, model)


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
    """
    Read concordant severe VCF files in parallel WITHOUT annotation.
    Annotation happens later only for passing variants.
    """
    available_cpus = get_available_cpus()
    num_parallel_tasks = min(len(vcfs), 6)  # Increased back to 6 since no annotation overhead
    cpus_per_task = max(1, len(available_cpus) // num_parallel_tasks)
    
    logger.info(f"Reading {len(vcfs)} concordant VCFs (no annotation): {num_parallel_tasks} parallel tasks")
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
    """Read discordant sibling pair VCFs WITHOUT annotation."""
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
def annotate_passing_variants(final_pass, alphamissense_annotator, gnomad_annotator):
    """
    Annotate only the variants that passed all filtering.
    This is MUCH more efficient than annotating everything.
    
    Args:
        final_pass: Dict of passing variants
        alphamissense_annotator: AlphaMissenseAnnotator instance
        gnomad_annotator: GnomADAnnotator instance
    
    Returns:
        dict: final_pass with annotations added
    """
    logger.info(f"\n{'='*60}")
    logger.info(f"Annotating {len(final_pass)} passing variants...")
    logger.info(f"{'='*60}\n")
    
    # Load annotation data (lazy loading - only happens here)
    if alphamissense_annotator:
        alphamissense_annotator.load_alphamissense_data()
    
    if gnomad_annotator:
        gnomad_annotator.load_gnomad_data()
    
    # Collect all variants to annotate
    alphamissense_batch = []
    gnomad_batch = []
    
    for var_id, (patients, zygosity, annlist, chrom, pos, ref, alt) in final_pass.items():
        # Get transcript ID from first annotation (index 6 is transcript ID in SnpEff ANN)
        transcript_id = annlist[0][6] if annlist and len(annlist[0]) > 6 else None
        alphamissense_batch.append((chrom, pos, ref, alt, transcript_id))
        gnomad_batch.append((chrom, pos, ref, alt))
    
    logger.info(f"Batch annotating {len(alphamissense_batch)} variants with AlphaMissense...")
    logger.info(f"  Strategy: Try transcript-specific first, fallback to variant-level")
    logger.info(f"Batch annotating {len(gnomad_batch)} variants with gnomAD...")
    
    # Batch annotate
    alphamissense_results = {}
    gnomad_results = {}
    
    if alphamissense_annotator:
        alphamissense_results = alphamissense_annotator.annotate_batch(alphamissense_batch)
    
    if gnomad_annotator:
        gnomad_results = gnomad_annotator.annotate_batch(gnomad_batch)
    
    # Add annotations to final_pass
    annotated_pass = {}
    for var_id, (patients, zygosity, annlist, chrom, pos, ref, alt) in final_pass.items():
        # Get AlphaMissense annotation
        transcript_id = annlist[0][6] if annlist and len(annlist[0]) > 6 else None
        alpham_key = f"{chrom}:{pos}:{ref}:{alt}:{transcript_id}"
        alpham_score = alphamissense_results.get(alpham_key, np.nan)
        
        # Get gnomAD annotation
        gnomad_key = f"{chrom}:{pos}:{ref}>{alt}"
        gnomad_af = gnomad_results.get(gnomad_key, None)
        
        annotated_pass[var_id] = [patients, zygosity, annlist, alpham_score, gnomad_af]
    
    logger.info(f"Annotation complete!")
    
    return annotated_pass


@timing
def parallel_recessive_modifier_v2(discordant_sib_pairs, concordant_severe_sibs, model, 
                                   alphamissense_annotator=None, gnomad_annotator=None):
    """
    Main analysis function with LAZY ANNOTATION.
    
    Flow:
    1. Parse VCFs without annotation (fast)
    2. Apply genetic model filtering (fast)
    3. Apply multi-sib filtering (fast)
    4. Annotate ONLY passing variants (efficient)
    
    This is much more efficient than annotating everything upfront.
    """
    logger.info(f"\n{'='*60}")
    logger.info(f"Starting LAZY ANNOTATION analysis with {model} model")
    logger.info(f"Discordant pairs: {len(discordant_sib_pairs)}")
    logger.info(f"Concordant severe siblings: {len(concordant_severe_sibs)}")
    logger.info(f"Current memory usage: {get_memory_usage():.2f} GB")
    logger.info(f"{'='*60}\n")
    
    # Step 1: Load concordant severe variants WITHOUT annotation
    logger.info("PHASE 1: Reading VCFs without annotation (fast)...")
    logger.info("Filtering for: intron, 5'UTR, 3'UTR, exon variants only")
    concord = parallel_read_vcf_concordant(concordant_severe_sibs, model)
    
    # Initialize tracking structures
    res = {}
    gene_to_var_map = defaultdict(list)
    var_to_genes_map = defaultdict(set)
    variant_sib_map = defaultdict(set)
    
    # Step 2: Process each discordant pair WITHOUT annotation
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
                    gene_to_var_map[gene].append(var)
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
    del res, gene_to_var_map, var_to_genes_map, variant_sib_map, concord
    gc.collect()
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Variants passing all filters: {len(final_pass)}")
    logger.info(f"Now will annotate ONLY these {len(final_pass)} variants")
    logger.info(f"{'='*60}\n")
    
    # Step 4: Annotate ONLY passing variants
    logger.info("PHASE 4: Annotating passing variants...")
    annotated_pass = annotate_passing_variants(final_pass, alphamissense_annotator, gnomad_annotator)
    
    del final_pass
    gc.collect()
    
    # Step 5: Build final DataFrame
    logger.info(f"\n{'='*60}")
    logger.info(f"PHASE 5: Building final results table...")
    logger.info(f"{'='*60}\n")
    
    final_df = pd.DataFrame.from_dict(annotated_pass, orient='index').reset_index()
    final_df = final_df.rename(columns={
        'index': 'variant_id',
        0: 'patients',
        1: 'zygosity',
        2: 'annotation',
        3: 'alphamissense_pathogenicity',
        4: 'gnomad_AF_NFE'
    })
    
    final_df[['#CHROM', 'POS', 'info', 'zygosity']] = final_df['variant_id'].str.split('|', expand=True)
    final_df[['REF', 'ALT']] = final_df['info'].str.split('>', expand=True)
    
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
        "zygosity", "gene", "patients", "N_sibs", "alphamissense_pathogenicity", "gnomad_AF_NFE"
    ]].reset_index(drop=True)
    
    # Summary statistics
    N_homo = final_df.query('zygosity == "homo"').shape[0]
    N_het = final_df.query('zygosity == "het"').shape[0]
    N_with_alpham = final_df['alphamissense_pathogenicity'].notna().sum()
    N_with_gnomad = final_df['gnomad_AF_NFE'].notna().sum()
    
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
    logger.info(f"  - With AlphaMissense annotation: {N_with_alpham} ({N_with_alpham/(N_homo+N_het)*100:.1f}%)")
    logger.info(f"  - With gnomAD AF (NFE): {N_with_gnomad} ({N_with_gnomad/(N_homo+N_het)*100:.1f}%)")
    logger.info(f"Final memory usage: {get_memory_usage():.2f} GB")
    logger.info(f"{'='*60}\n")
    
    return final_df


####################################
def get_args():
    """Parse command line parameters."""
    parser = argparse.ArgumentParser(
        description='Recessive Modifier Analysis with Transcript-First AlphaMissense Lookup',
        epilog='AlphaMissense: First tries transcript-specific, then falls back to variant-level'
    )
    parser.add_argument("-model", type=str, required=True, choices=['dominant', 'recessive'])
    parser.add_argument("-alphamissense", type=str, default=None, 
                        help="Path to AlphaMissense TSV file with transcript_id column (can be gzipped)")
    parser.add_argument("-gnomad", type=str, default=None)
    parser.add_argument("-output", type=str, default=None)
    
    return parser.parse_args()


def main():
    """Main execution function."""
    args = get_args()
    
    logger.info(f"Starting analysis with TRANSCRIPT-FIRST AlphaMissense lookup")
    logger.info(f"Initial memory usage: {get_memory_usage():.2f} GB")
    
    # Initialize annotators (but don't load data yet!)
    alphamissense_annotator = None
    if args.alphamissense:
        logger.info("Initializing AlphaMissense annotator (lazy loading)...")
        logger.info("  Lookup strategy: transcript-specific first, variant-level fallback")
        alphamissense_annotator = AlphaMissenseAnnotator(args.alphamissense)
    
    gnomad_annotator = None
    if args.gnomad:
        logger.info("Initializing gnomAD annotator (lazy loading)...")
        gnomad_annotator = GnomADAnnotator(args.gnomad)
    
    # Define file paths
    fp = '/fs/ess/PAS0631/2024_august_recessive_modifier/SNPEFF_VCF/SNPEFF'
    
    discordant_sibs = [
        
        # s02_007
        [
            os.path.join(fp, 'S02_007_1.dv.filter.snpEff.vcf') ,
            [ os.path.join(fp, 'S02_008_1.dv.filter.snpEff.vcf'), os.path.join(fp, 'S02_009_1.dv.filter.snpEff.vcf') ]
        ],
        
        # s_65811
        [ os.path.join(fp, 's_65811_1.dv.filter.snpEff.vcf'), [os.path.join(fp, 's_65812_1.dv.filter.snpEff.vcf')] ],
        
        # s_50
        [ os.path.join(fp, 'S51_1.dv.filter.snpEff.vcf'), [os.path.join(fp, 'S50_1.dv.filter.snpEff.vcf')] ],
    
        # s04_001
        [ os.path.join(fp, 'S04_001_1.dv.filter.snpEff.vcf'), [os.path.join(fp, 'S04_002_1.dv.filter.snpEff.vcf')] ],

        # s_13673
        [ os.path.join(fp, 's_13673_1.dv.filter.snpEff.vcf'), [os.path.join(fp, 's_13674_1.dv.filter.snpEff.vcf')] ], 

        # s_113
        [ os.path.join(fp, 's_113_1.dv.filter.snpEff.vcf'), [os.path.join(fp, 's_112_1.dv.filter.snpEff.vcf')] ], 

        # s_66809
        [ os.path.join(fp, 's_66809_1.dv.filter.snpEff.vcf'), [os.path.join(fp, 's_66189_1.dv.filter.snpEff.vcf')] ]
    ]
    

    concordant_sibs = [
        os.path.join(fp, 'S45_1.dv.filter.snpEff.vcf'),
        os.path.join(fp, 'S46_1.dv.filter.snpEff.vcf'),
        os.path.join(fp, 's_79_1.dv.filter.snpEff.vcf'), 
        os.path.join(fp, 's_80_1.dv.filter.snpEff.vcf'), 
        os.path.join(fp, 'S57_1.dv.filter.snpEff.vcf'),
        os.path.join(fp, 'S58_1.dv.filter.snpEff.vcf'),
        os.path.join(fp, 'S134_2_1.dv.filter.snpEff.vcf'), 
        os.path.join(fp, 'S136_1.dv.filter.snpEff.vcf'),
        os.path.join(fp, 'S190_1.dv.filter.snpEff.vcf'),
        os.path.join(fp, 'S192_1.dv.filter.snpEff.vcf')
    ]
    
    # Run analysis
    res = parallel_recessive_modifier_v2(
        discordant_sibs,
        concordant_sibs,
        args.model,
        alphamissense_annotator,
        gnomad_annotator
    )
    
    # Save results
    output_file = args.output or f"{timestamp}_results_{args.model}.v2.vcf"
    res.to_csv(output_file, index=False, sep="\t", header=True)
    logger.info(f"\nResults saved to: {output_file}")
    logger.info(f"Peak memory usage: {get_memory_usage():.2f} GB")


if __name__ == "__main__":
    main()