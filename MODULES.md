# Module Reference — `haplotype_aware` pipeline

Documents every process actually wired into `workflows/haplotype_aware/main.nf` (directly or via a
`subworkflows/*.nf` include), grouped by pipeline stage in roughly execution order. Each entry lists
the process name, its module file, what it does, and its declared `input`/`output` tuples.

Processes that exist in `modules/` but are **not** called anywhere in `main.nf`'s workflow body are
listed separately at the end ("Included but not invoked") rather than omitted, since they're still
`include`d and easy to mistake for active pipeline steps.

Legend: **in** = process `input:` block, **out** = process `output:` block (emit names in `` `code` ``).

---

## 1. Input parsing & basecalling

### `DESIGN_INPUT` — `modules/parse_design/main.nf`
Parses `--design` (SAMPLE/ONT/ONT_UNAL_BAM tsv) into a per-row csv via `bin/parse_design.py`.
- **in**: `val(design)` — path to the design tsv
- **out**: `path("fastq.csv")` → `fastq_ch` — one row per sample/lane

### `MERGE_INPUTS` — `modules/samtools/main.nf`
Merges multi-lane fastq/bam rows for a sample into one fastq + one bam (`samtools merge` + index).
- **in**: `tuple val(sample), path(fastqs), path(bams)` — grouped lists, one sample
- **out**: `tuple val(sample), path("${sample}.fastq.gz"), path("${sample}.bam")` → `merged_inputs_ch`

### `DORADO_BASECALL` — `modules/dorado/main.nf`
Basecalls a pod5 directory with `dorado basecaller` (used only for `--pod5_design` runs). Model string
(`--dorado_basecall_model`, default `sup,5mCG_5hmCG`) carries the mod-calling tag that embeds MM/ML
methylation tags.
- **in**: `tuple val(sampleID), val(pod5_dir)`
- **out**: `tuple val(sampleID), path("*.unal.bam")` → `dorado_basecall_ch`

### `DORADO_TRIM` — `modules/dorado/main.nf`
Trims sequencing adapters from the unaligned bam with `dorado trim`.
- **in**: `tuple val(sampleID), val(fastq), val(unal_bam)`
- **out**: `tuple val(sampleID), path("*dorado.trim.bam")` → `dorado_trim_output_ch`

### `SAMTOOLS_CONVERT_BAM_TO_FASTQ` — `modules/samtools/main.nf`
Converts the trimmed unaligned bam to fastq (`samtools fastq -T MM,ML,MN` preserves methylation tags
as fastq comments for `MINIMAP2_ALIGN -y` to restore later).
- **in**: `tuple val(sampleID), val(unal_bam)`
- **out**: `tuple val(sampleID), val(unal_bam), path("*.fastq.gz")` → `samtool_convert_unalbam2fastq`
  — this becomes the pipeline's canonical `ont_ch` shape: `(sampleID, unal_bam, fastq)`

### `NANOPLOT_RAW` — `modules/nanoplot/main.nf`
Raw-read QC (`NanoPlot`) on the pre-HERRO fastq. Report-only, no downstream consumers.
- **in**: `tuple val(sampleID), val(bam), val(fastq)`
- **out**: `path("*")` — NanoPlot report directory

### `INDEX_REFERENCE` — `modules/samtools/main.nf`
`samtools faidx` on the reference genome.
- **in**: `path(ref)`
- **out**: `tuple path(ref), path("${ref}.fai")` → `ref_indexed_ch` — passed everywhere downstream as
  `samtools_fai_index`

### `MINIMAP2_INDEX` — `modules/minimap2/main.nf`
Builds (or reuses, if already cached at `index_dir`) the minimap2 `.mmi` index for the reference.
- **in**: `path genome`, `val index_dir`
- **out**: `path("*.mmi")` → `minimap_index`

---

## 2. Error correction (HERRO) — `subworkflows/herro.nf`

GPU deep-learning ONT self-correction ([lbcb-sci/herro](https://github.com/lbcb-sci/herro)) — **not**
the same tool as the earlier "HERO" (`kangxiongbin/HERO`) hybrid Illumina+ONT corrector this replaced.
herro corrects ONT reads via all-vs-all self-alignment + neural-net inference; it never touches
Illumina data, so unlike the old path there's no per-sample gating — every sample in `ont_ch` runs
through the same three steps.

### `HERRO_PREPROCESS` — `modules/herro/main.nf`
Adapter-trims/splits/length-filters the ONT reads (`herro`'s `scripts/preprocess.sh`, wrapping
porechop + duplex_tools + seqkit). `--herro_preprocess_parts` chunks the job to bound memory on large
read sets.
- **in**: `tuple val(sampleID), val(unal_bam), val(fastq)`
- **out**: `tuple val(sampleID), path("${sampleID}.herro_preprocessed.fastq.gz")` → `preprocessed_ch`

### `HERRO_ALIGN_BATCHES` — `modules/herro/main.nf`
All-vs-all self-alignment of the preprocessed reads (`minimap2`) via herro's
`scripts/create_batched_alignments.sh`, batched for the GPU inference step.
- **in**: `tuple val(sampleID), path(preprocessed_fastq)`
- **out**: `tuple val(sampleID), path(preprocessed_fastq), path("${sampleID}_alignment_batches")` → `batches_ch`

### `HERRO_INFERENCE` — `modules/herro/main.nf`
Neural-net correction (`herro inference`, requires GPU + `--herro_model` weights). Output is FASTA
re-consensus, gzipped under a `.fastq.gz` name so downstream extension-based basename logic keeps
working (hifiasm/minimap2 auto-detect FASTA vs FASTQ by content, not extension).
- **in**: `tuple val(sampleID), path(preprocessed_fastq), path(alignment_batches)`
- **out**: `tuple val(sampleID), path("*.herro.corrected.fastq.gz")` → `herro_ch`

### `SAMTOOLS_FASTQ_TO_BAM` — `modules/samtools/main.nf`
Reimports the herro-corrected fastq into an unaligned bam, reusing the original `@RG` line (needed
later for `dorado polish`'s model auto-resolution).
- **in**: `tuple val(sampleID), path(fastq), path(orig_bam)`
- **out**: `tuple val(sampleID), path("*.unal.bam"), path(fastq)` → `fastq_to_bam_ch`

---

## 3. Trio / haplotype family structure — `subworkflows/trio.nf`, gated by `--run_trio_analysis` + `--family_json`

Runs only for parent samples identified via `family_json`; splits them out of the main sample channel.

### `YAK_COUNT` — `modules/yak/main.nf`
Builds a k-mer database (`yak count -k31 -b37`) from a parent's own ONT reads, for hifiasm trio-binning.
- **in**: `tuple val(sampleID), val(fastqs)` — `fastqs` may be one ONT fastq or an Illumina R1/R2 pair
- **out**: `tuple val(sampleID), path("*.yak")` → `yak_ch`

Parents also run `minimap2` (§5) + `CLAIR3` (§5) to produce a genotype-only VCF for pedigree phasing —
same processes the rest of the pipeline uses, just applied to parent-only reads and not carried further.

---

## 4. Assembly — `subworkflows/ont.nf` (`hifasim`/`hifasim_trio`/`gfatools`/`dorado`), gated by `--assemble_genome`

### `HIFASIM` — `modules/hifasim/main.nf`
De novo ONT assembly (`hifiasim --ont`) for non-trio samples.
- **in**: `tuple val(sampleID), val(bam), val(reads)`
- **out**: `tuple val(sampleID), val(bam), val(reads), path("*hap1.p_ctg.gfa"), path("*hap2.p_ctg.gfa")` → `hifasim_asm`

### `HIFASIM_TRIO` — `modules/hifasim/main.nf`
Trio-binned assembly (`hifiasm --ont -1 <paternal.yak> -2 <maternal.yak>`) for children with both
parents' yak databases resolved.
- **in**: `tuple val(sampleID), val(bam), val(reads), val(pat_yak), val(mat_yak)`
- **out**: `tuple val(sampleID), val(bam), val(reads), path("*hap1.p_ctg.gfa"), path("*hap2.p_ctg.gfa")` → `hifasim_trio_asm`

### `GFA_CONVERT` — `modules/gfatools/main.nf`
Converts each haplotype's GFA to FASTA (`gfatools gfa2fa`) and concatenates into one combined fasta
with `hap1_`/`hap2_` prefixed headers.
- **in**: `tuple val(sampleID), val(bam), val(reads), path(hap1), path(hap2)` — GFAs
- **out**: `tuple val(sampleID), val(bam), val(reads), path("*.haps.combined.fa")` → `combined_fa`

### `GFA_FAIDX` — `modules/gfatools/main.nf`
Indexes the combined fasta and re-splits it back into per-haplotype fastas + `.fai`.
- **in**: `tuple val(sampleID), val(bam), val(reads), path(combined_fa)`
- **out**: `tuple val(sampleID), val(bam), val(reads), path(combined_fa), path("${combined_fa}.fai")` → `fasta_asm`;
  `tuple val(sampleID), val(bam), val(reads), path("*hap1.fa"), path("*hap2.fa")` → `haplotype_fasta_asm`

### `DORADO` — `modules/dorado/main.nf`
Polishes the assembly (`dorado polish`) by realigning the original reads to it, gated by `--dorado`
(default `true`).
- **in**: `tuple val(sampleID), val(unal_bam), val(reads), val(hapfasta), val(hapfai)`
- **out**: `tuple val(sampleID), val(reads), path("*hap1.doradopolish.fa"), path("*hap2.doradopolish.fa")` → `dorado_output_ch`

### `GENOME_STATS` — `modules/genome_stats/main.nf`
Read/base counts vs. assembled genome size and implied coverage — a QC tsv, no downstream consumers.
- **in**: `tuple val(sampleID), val(reads), val(hap1_fa), val(hap2_fa)`
- **out**: `tuple val(sampleID), path("*.genome_stats.tsv")` → `genome_stats_ch`

### `MAP_ASM_TO_REF` — `modules/minimap2/main.nf`
Maps each haplotype assembly back to the reference (`minimap2 -ax asm20`), for downstream
per-chromosome coverage QC.
- **in**: `path ref_mmi`, `tuple val(sampleID), val(reads), val(hap1_fa), val(hap2_fa)`
- **out**: `tuple val(sampleID), path(hap1.ref.bam+bai), path(hap2.ref.bam+bai)` → `hap_bams`

### `CHROM_COVERAGE` — `modules/genome_stats/main.nf`
Per-chromosome `samtools coverage` on both haplotype-vs-reference bams from `MAP_ASM_TO_REF`.
- **in**: `tuple val(sampleID), path(hap1_bam), path(hap1_bai), path(hap2_bam), path(hap2_bai)`
- **out**: `tuple val(sampleID), path("*.chrom_coverage.tsv")` → `chrom_coverage_ch`

### `DIPCALL` — `modules/dipcall/main.nf`
Assembly-based small-variant + structural calling (`run-dipcall`) from the two polished haplotypes.
- **in**: `tuple val(sampleID), val(reads), val(hap1), val(hap2)`, `tuple path(ref), path(ref_fai)`
- **out**: `tuple val(sampleID), path("*.dip.vcf.gz"), path("*.dip.vcf.gz.tbi")` → `dipcall_ch`;
  plus `hap1_bam`/`hap2_bam` emits

### `HAPDIFF` — `modules/hapdiff/main.nf`
Assembly-based SV calling (`hapdiff.py`) against the reference, phased + unphased VCFs.
- **in**: `tuple val(sampleID), val(reads), val(hap1), val(hap2)`, `tuple path(ref), path(ref_fai)`
- **out**: `path("*.hapdiff.phased.vcf.gz")` → `hapdiff_ch`; `path("*.hapdiff.unphased.vcf.gz")` → `hapdiff_unphased_ch`

---

## 5. Alignment-based variant calling — `subworkflows/ont.nf`, gated by `--alignment_based_variant_calling` (default `true`)

### `MINIMAP2_ALIGN` — `modules/minimap2/main.nf`
Aligns corrected ONT reads to the reference (`minimap2 -x map-ont -y --MD`); `-y` restores the MM/ML
methylation tags `SAMTOOLS_CONVERT_BAM_TO_FASTQ` carried through as fastq comments.
- **in**: `path index`, `tuple val(sampleID), val(fastq)`
- **out**: `tuple val(sampleID), path("*.sorted.bam"), path("*.sorted.bam.bai")` → `bam_ch`

### `CLAIR3` — `modules/clair3/main.nf`
ONT small-variant calling (`run_clair3.sh`, GPU, `r941_prom_sup_g5014` model).
- **in**: `tuple val(ref), val(ref_fai)`, `tuple val(sampleID), val(bam), val(bai)`
- **out**: `tuple val(sampleID), path("*.clair3.vcf.gz"), path("*.clair3.vcf.gz.tbi")` → `clair3_ch`

### `LONGPHASE` — `modules/longphase/main.nf`
First-pass phasing + haplotagging from clair3 SNP/indel calls alone (`longphase phase` + `haplotag`).
- **in**: `tuple val(ref_fa), val(ref_fai)`, `tuple val(sampleID), val(vcf), val(tbi), val(bam), val(bai)`
- **out**: haplotagged bam → `longphase_ch`; phased vcf → `longphase_vcf_ch`; combined → `longphase_full_ch`

### `SNIFFLES` — `modules/sniffles/main.nf`
Structural-variant calling from the longphase-haplotagged bam.
- **in**: `tuple val(sampleID), val(bam), val(bai)`
- **out**: `tuple val(sampleID), path("*.sniffles.vcf.gz"), path("*.sniffles.vcf.gz.tbi")` → `sniffles_ch`

### `LONGPHASE_SV` — `modules/longphase/main.nf`
Second-pass phasing that additionally incorporates the sniffles SV calls (`--sv-file`), re-haplotagging
the bam. This is the "final" haplotagged bam most downstream tools (straglr, spectre, paraphase,
modkit) consume.
- **in**: `tuple val(ref_fa), val(ref_fai)`, `tuple val(sampleID), val(snp_vcf), val(snp_tbi), val(sv_vcf), val(sv_tbi), val(bam), val(bai)`
- **out**: haplotagged bam → `longphase_sv_ch`; phased vcf → `longphase_sv_vcf_ch`; combined → `longphase_sv_full_ch`

### `STRAGLR` — `modules/straglr/main.nf`
Repeat-expansion genotyping (`straglr.py`) on the final haplotagged bam.
- **in**: `tuple val(sampleID), val(bam), val(bai)`, `tuple val(ref_fa), val(ref_fai)`
- **out**: `tuple val(sampleID), path("*.straglr.tsv"), path("*.straglr.bed")` → `straglr_ch`

### `SPECTRE` — `modules/spectre/main.nf`
CNV calling: `mosdepth` coverage windows fed into `spectre CNVCaller`.
- **in**: `tuple val(sampleID), val(bam), val(bai)`, `tuple val(ref_fa), val(ref_fai)`
- **out**: `tuple val(sampleID), path("*.spectre.*")` → `spectre_cnv_ch`; mosdepth regions bed

### `PARAPHASE` — `modules/paraphase/main.nf`
Segmental-duplication region reconstruction (e.g. SMN1/SMN2) via `paraphase`, restricted to
`--paraphase_genes` if set.
- **in**: `tuple val(sampleID), val(bam), val(bai)`, `tuple val(ref_fa), val(ref_fai)`
- **out**: `tuple val(sampleID), path("*.paraphase.json"), path("*.paraphase.bam+bai"), path("*_paraphase_vcfs")` → `paraphase_ch`

### `MODKIT_PILEUP` — `modules/modkit/main.nf`
Per-haplotype methylation pileup (`modkit pileup --partition-tag HP`) on the final haplotagged bam.
Silently produces empty output if MM/ML tags were ever dropped upstream (see `SAMTOOLS_CONVERT_BAM_TO_FASTQ`/`MINIMAP2_ALIGN` notes above).
- **in**: `tuple val(sampleID), val(bam), val(bai)`, `tuple val(ref), val(ref_fai)`
- **out**: `tuple val(sampleID), path("${sampleID}*.bed.gz")` → `modkit_ch`

---

## 6. Variant consolidation

### `CONSOLIDATE_VARIANTS` — `modules/consolidate_variants/main.nf`
Merges per-caller VCFs (clair3+longphase, sniffles, dipcall, hapdiff) into one master per-sample VCF
via `bin/consolidate_variants.py`, tagging each variant with which caller(s) found it.
- **in**: `tuple val(sampleID), val(clair3_longphase_vcf), val(sniffles_vcf), val(dipcall_vcf), val(hapdiff_vcf)`
  — assembly-caller slots are `"NO_FILE"` placeholders when `--assemble_genome` is off
- **out**: `tuple val(sampleID), path("*.master.vcf")` → `consolidated_ch`

---

## 7. Trio-aware pedigree phasing (children only)

### `WHATSHAP_PHASE_TRIO` — `modules/whatshap/main.nf`
Additive refinement of the child's small-variant phasing using Mendelian constraints from both
parents' genotype VCFs (`whatshap phase --ped`). Runs alongside, not instead of, `LONGPHASE_SV` —
a separate, additionally-published output, not fed back into consolidation.
- **in**: `tuple val(ref), val(ref_fai)`, `tuple val(child), val(father), val(mother), val(child_vcf), val(child_bam), val(child_bai), val(father_vcf), val(mother_vcf)`
- **out**: `tuple val(child), path("*.trio.whatshapphase.vcf.gz"), path("*.trio.whatshapphase.vcf.gz.tbi")` → `whatshap_phase_trio_ch`

---

## 8. Variant annotation — gated by `--run_alphagenome` (scaffold, default `false`)

### `ALPHAGENOME` — `modules/alphagenome/main.nf`
Variant-effect annotation of the master VCF via `bin/run_alphagenome.py`.
- **in**: `tuple val(sampleID), val(vcf)`, `tuple path(ref), path(ref_fai)`
- **out**: `tuple val(sampleID), path("*.alphagenome.tsv")` → `alphagenome_ch`

---

## Included but not invoked

These are `include`d in `main.nf` (via `subworkflows/illumina.nf`) but never actually called in the
workflow body — present in case a future Illumina-only branch is wired in, not part of any current run:

- `ATRIA` — `modules/atria/main.nf` (lost its only caller when `hero_correction` was replaced by
  `herro_correction`, which doesn't touch Illumina reads at all)
- `BWA_INDEX` / `BWA_MEM` — `modules/bwa/main.nf`
- `FILTER_BAM` / `INDEX_BAM` — `modules/samtools/main.nf` (via `illumina.nf`'s `filter_bam` workflow)
- `DEEPVARIANT_CALL_VARIANTS` / `SPLIT_DV_VCF` — `modules/deepvariant/main.nf`
- `BAM_TO_BW` / `BAM_COMPARE` / `BW_COMPARE` / `MERGE_BW` — `modules/deeptools/main.nf`

`WHATSHAP_PHASE` / `WHATSHAP_HAPLOTAG` (`modules/whatshap/main.nf`) are similarly unused — only
`WHATSHAP_PHASE_TRIO` is called.
