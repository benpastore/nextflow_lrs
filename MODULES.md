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

## 4. Assembly — `subworkflows/assembly.nf`, gated by `--assemble_genome` (default `false` — currently deactivated)

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

## 5. Alignment-based variant calling — `subworkflows/ont.nf` (always on)

How each caller here ends up with (or without) haplotype information, since it differs per tool:

| Caller | Haplotype info | How |
|---|---|---|
| Clair3 | Real, native phase | Phased downstream by `LONGPHASE`/`LONGPHASE_SV` (`GT` with `\|` + `PS`), not by Clair3 itself |
| Sniffles | Real, native phase | `--phase` reads the HP/PS tags already on its input bam, writes `PHASE=` in `INFO` |
| Straglr | Post-hoc annotation | No native mode — `HAPLOTYPE_ANNOTATE_STRAGLR` tags calls from HP-tagged reads afterward |
| Spectre | Post-hoc annotation | Same as Straglr — `HAPLOTYPE_ANNOTATE_SPECTRE` |
| Dipcall / Hapdiff | Inherent (assembly-based) | Already phase-resolved by construction — called from the haplotype-split assembly itself, see §4 |

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
Structural-variant calling from the longphase-haplotagged bam, with `--phase` (reads the HP/PS tags
already on that bam and writes a per-SV `PHASE=` field in `INFO` — Sniffles doesn't need re-haplotagging
of its own, unlike Straglr/Spectre below). Emits a raw `.vcf`, not `.vcf.gz` — see `INDEX_SNIFFLES_VCF`.
- **in**: `tuple val(sampleID), val(bam), val(bai)`
- **out**: `tuple val(sampleID), path("*.sniffles.vcf")` → `sniffles_raw_ch`

### `INDEX_SNIFFLES_VCF` — `modules/sniffles/main.nf`
Bgzip + tabix-index Sniffles' raw VCF, under the `bcftools` container/label — the `sniffles` biocontainer
image ships only the `sniffles` binary itself, no htslib CLI tools, so compression can't happen inline
in `SNIFFLES` the way most other modules do it.
- **in**: `tuple val(sampleID), path(vcf)`
- **out**: `tuple val(sampleID), path("*.sniffles.vcf.gz"), path("*.sniffles.vcf.gz.tbi")` → `sniffles_ch`

### `LONGPHASE_SV` — `modules/longphase/main.nf`
Second-pass phasing that additionally incorporates the sniffles SV calls (`--sv-file`), re-haplotagging
the bam. This is the "final" haplotagged bam most downstream tools (straglr, spectre, paraphase,
modkit) consume.
- **in**: `tuple val(ref_fa), val(ref_fai)`, `tuple val(sampleID), val(snp_vcf), val(snp_tbi), val(sv_vcf), val(sv_tbi), val(bam), val(bai)`
- **out**: haplotagged bam → `longphase_sv_ch`; phased vcf → `longphase_sv_vcf_ch`; combined → `longphase_sv_full_ch`

### `STRAGLR` — `modules/straglr/main.nf`
Repeat-expansion genotyping (`straglr.py`) on the final haplotagged bam. No native haplotype-aware
genotyping mode (only a private ONT fork of Straglr has one) — see `HAPLOTYPE_ANNOTATE_STRAGLR` below.
- **in**: `tuple val(sampleID), val(bam), val(bai)`, `tuple val(ref_fa), val(ref_fai)`
- **out**: `tuple val(sampleID), path("*.straglr.tsv"), path("*.straglr.bed")` → `straglr_ch`

### `HAPLOTYPE_ANNOTATE_STRAGLR` — `modules/haplotype_annotate/main.nf`
Post-hoc haplotype tagging for Straglr calls, via `bin/annotate_haplotype.py`: tallies `HP:1`/`HP:2`
read counts (`samtools view -c -d HP:...`) in a flank window (`--haplotype_flank_bp`, default 1000bp)
around each locus's start/end, and appends an `HP1`/`HP2`/`AMBIGUOUS` column — rather than re-genotyping
on a haplotype-split bam, which would halve read depth and hurt sensitivity.
- **in**: `tuple val(sampleID), path(tsv), path(bed), path(bam), path(bai)` — straglr's own tsv/bed joined
  with the final haplotagged bam (`longphase_sv`'s `longphase_sv_ch`)
- **out**: `tuple val(sampleID), path("*.straglr.haplotagged.tsv")` → `straglr_haplotagged_ch`

### `SPECTRE` — `modules/spectre/main.nf`
CNV calling: `mosdepth` coverage windows fed into `spectre CNVCaller`. Normalizes whatever
`spectre CNVCaller` names its own output files (varies by version/compression state) into three
predictable emits.
- **in**: `tuple val(sampleID), val(bam), val(bai)`, `tuple val(ref_fa), val(ref_fai)`
- **out**: `tuple val(sampleID), path("*.spectre.vcf.gz"), path("*.spectre.vcf.gz.tbi")` → `spectre_vcf_ch`;
  `tuple val(sampleID), path("*.spectre.bed.gz"), path("*.spectre.bed.gz.tbi")` → `spectre_bed_ch`;
  `tuple val(sampleID), path("*.spectre.spc")` → `spectre_spc_ch`; plus mosdepth regions bed

### `HAPLOTYPE_ANNOTATE_SPECTRE` — `modules/haplotype_annotate/main.nf`
Same rationale/mechanism as `HAPLOTYPE_ANNOTATE_STRAGLR` above, applied to Spectre's CNV segments (flank
around each segment's start/end, not the whole span — CNV calls can be 100kb-2Mb, and boundary-flanking
reads are what's actually informative about haplotype).
- **in**: `tuple val(sampleID), path(bed_gz), path(bed_gz_tbi), path(bam), path(bai)` — `spectre_bed_ch`
  joined with the final haplotagged bam
- **out**: `tuple val(sampleID), path("*.spectre.haplotagged.bed.gz"), path("*.spectre.haplotagged.bed.gz.tbi")` → `spectre_haplotagged_ch`

### `CNV_COMPARE` — `modules/cnv_repeat_compare/main.nf`
Cohort-wide, gated by `--family_json` + `--run_cnv_compare`. `bin/cnv_compare.py` turns every child's
Spectre VCF into a per-gene copy number (genes from `--cnv_repeat_gtf`, default `--snpeff_gtf`; no call = CN 2),
then compares each mild child with its severe siblings (family.json mild/severe grouping).
- **in**: family.json, a `sampleID<TAB>vcf` manifest + every `*.spectre.vcf.gz`, the GTF
- **out** (`10_modifier_analysis/cnv`): `cohort.cnv_gene_summary.tsv` (per gene: `n_mild_higher` = mild children
  with more copies than every severe sibling, `n_mild_lower`, per-pair counts, sign test, concordant-severe background
  gains/losses), `cohort.cnv_mild_vs_severe.tsv` (one row per differing mild/severe pair x gene),
  `cohort.cnv_gene_cn.tsv` (gene x child CN matrix)

### `REPEAT_COMPARE` — `modules/cnv_repeat_compare/main.nf`
Cohort-wide, gated by `--family_json` + `--run_repeat_compare`. `bin/repeat_compare.py` merges every child's
Straglr BED into cohort loci by position only (overlapping or within `--repeat_merge_dist`; Straglr's per-sample
consensus motif varies between children for long VNTR motifs), takes each child's longest allele size in bp expressed
in copies of the locus's most common motif, annotates each locus with its most direct gene feature (CDS > 5'UTR > 3'UTR > ncRNA exon >
promoter > intron, else nearest gene) and optional `--repeat_known_loci_bed` names, then compares mild vs severe.
- **in**: family.json, a `sampleID<TAB>bed` manifest + every `*.straglr.bed`, the GTF (or `NO_FILE`), known loci (or `NO_FILE`)
- **out** (`10_modifier_analysis/repeat_expansion`): `cohort.repeat_locus_summary.tsv` (within-family `n_mild_longer`/`n_mild_shorter`,
  `n_mild_specific_expanded` = mild children expanded where no severe sibling is, group medians, Mann-Whitney + BH q),
  `cohort.repeat_gene_summary.tsv`, `cohort.repeat_pattern_summary.tsv` (mild-longer vs -shorter by feature class,
  motif period, motif, known locus; per-child expansion burden mild vs siblings), `cohort.repeat_mild_vs_severe.tsv`,
  `cohort.repeat_sample_burden.tsv`, `cohort.repeat_loci.tsv` (catalog with every child's alleles)

### `PARAPHASE` — `modules/paraphase/main.nf`
Segmental-duplication region reconstruction (e.g. SMN1/SMN2) via `paraphase`, restricted to
`--paraphase_genes` if set.
- **in**: `tuple val(sampleID), val(bam), val(bai)`, `tuple val(ref_fa), val(ref_fai)`
- **out**: `tuple val(sampleID), path("*.paraphase.json"), path("*.paraphase.bam+bai"), path("*_paraphase_vcfs")` → `paraphase_ch`

### `MODKIT_PILEUP` — `modules/modkit/main.nf`
Per-haplotype methylation pileup (`modkit pileup --phased`, modkit >= 0.6) on the final haplotagged bam.
Silently produces empty output if MM/ML tags were ever dropped upstream (see `SAMTOOLS_CONVERT_BAM_TO_FASTQ`/`MINIMAP2_ALIGN` notes above).
- **in**: `tuple val(sampleID), val(bam), val(bai)`, `tuple val(ref), val(ref_fai)`
- **out**: `tuple val(sampleID), path("${sampleID}*.bed.gz")` → `modkit_ch`

---

## 6. Variant consolidation

### `CONSOLIDATE_VARIANTS` — `modules/consolidate_variants/main.nf`
Merges the **alignment-based** per-caller VCFs (phased small variants in the `clair3_longphase` column — WhatsHap
family-phased with parents, else longphase_sv — and Sniffles SVs) into one master per-sample VCF via
`bin/consolidate_variants.py`, tagging each variant with which caller(s) found it. Assembly-based calls
(dipcall / hapdiff, `--assemble_genome`) are kept separate and not merged in.
With parents (`--run_trio_analysis`), every site also gets **`INFO/ORIGIN`** = paternal / maternal / both / unknown —
the parental haplotype carrying the ALT — plus `ORIGIN_SOURCE` and a `DE_NOVO` flag:
small variants from the parents' Clair3 genotypes (`mendelian`), else the child's WhatsHap phase block whose HP1 was
assigned to a parent by its Mendelian-informative SNVs (`phase_block`); SVs from `SV_TRIO_ORIGIN`, matched by Sniffles
ID (`sv_trio_origin`). This is the one label small variants and SVs share (their PS values are not comparable).
- **in**: `tuple val(sampleID), val(clair3_longphase_vcf), val(sniffles_vcf), val(father_vcf), val(mother_vcf), val(sv_origin)`
  — parent/SV-origin slots are `"NO_FILE"` without trios
- **out**: `tuple val(sampleID), path("*.master.vcf")` → `consolidated_ch`; trio mode also
  `<sample>.phase_block_origin.tsv` (per phase block: informative SNVs, HP1's parent)

---

## 7. Joint family pedigree phasing

### `WHATSHAP_PHASE_FAMILY` — `modules/whatshap/main.nf`
ONE `whatshap phase --ped` run per family: every child of a father/mother pair + both parents, then split back
into one VCF per child. Genotypes: the children's `LONGPHASE_SV` VCFs + the parents' Clair3 VCFs, merged into
one multi-sample VCF (whatshap only phases the genotypes of its first positional VCF). Reads: every child's
phased BAM and — with `whatshap_family_parent_reads` (default `true`) — the parents' BAMs from `trio_family`;
each gets an `@RG SM:<sample>` when it lacks one (pedigree mode assigns reads by SM). Run per chromosome,
`whatshap_family_jobs` (6) at a time, then concatenated. Phasing all members together means the parents are
phased too, every child is phased against the same parental haplotypes, and Mendelian conflicts are judged
against all children. Logs het-SNV phasing rates per member. With `--run_trio_analysis` each child's VCF (ALT
sites only) **replaces** `LONGPHASE_SV`'s as its phased small-variant VCF for `CONSOLIDATE_VARIANTS` and
the hap-ibd steps. SVs aren't included (whatshap can't phase them) — their parent of origin comes from
`SV_TRIO_ORIGIN`, and both meet in `CONSOLIDATE_VARIANTS`' `INFO/ORIGIN`.
- **in**: `tuple val(ref), val(ref_fai)`, `tuple val(family), val(father), val(mother), val(father_vcf), val(mother_vcf), val(father_bam), val(mother_bam), val(children), val(child_vcfs), val(child_bams)`
- **out** (`06_variants/whatshap_family_phased`): `<child>.trio.whatshapphase.vcf.gz` (+ `.tbi`) per child →
  `child_vcfs_ch` (split to `(child, vcf, tbi)` in `main.nf`); `<father>_<mother>.family.whatshapphase.vcf.gz`
  (every member, parents phased) → `family_vcf_ch`

### `SV_TRIO_ORIGIN` — `modules/sv_trio_origin/main.nf`
Parent of origin for each child Sniffles SV (`bin/sv_trio_origin.py`). Parents get their own `SNIFFLES` run in
`subworkflows/trio.nf`. Two lines of evidence: **transmission** (SV also in father's / mother's Sniffles calls —
same type with INS~DUP compatible, start within `sv_merge_dist`, size ratio ≥ 0.5) and
**haplotype** (the SV's Sniffles `PHASE` hp/block, with each LONGPHASE round-1 block's HP1 assigned to a parent from
the child's phased het SNVs where the parents' genotypes differ; ≥ 3 informative SNVs, ≥ 80% agreement).
- **in**: `tuple val(child), val(child_sv), val(child_snv), val(father), val(mother), val(father_snv), val(mother_snv), val(father_sv), val(mother_sv)`
- **out** (`06_variants/sv_trio_origin`): `<child>.sv_trio_origin.tsv` (per SV: transmission, haplotype_origin, origin,
  origin_source = transmission+haplotype / transmission / haplotype / conflict), `<child>.sv_trio_blocks.tsv`
  (per block: informative SNVs, HP1 parent)

## 7b. Annotation, sibling sharing, methylation — `subworkflows/annotation.nf`, `subworkflows/haplotypes.nf`

Flow: `WHATSHAP_PHASE_FAMILY` + `SV_TRIO_ORIGIN` → `CONSOLIDATE_VARIANTS` (INFO/ORIGIN) → `SNPEFF` → `ANNOTATE_AF`;
with parents also `WHATSHAP_PHASE_FAMILY` family VCF → `ORIENT_PHASE` → `HAP_IBD` → `IBD_SEGMENTS` →
`MILD_HAPLOTYPE_VARIANTS` (+ `METHYLATION_COMPARE`). No genetic (recessive / dominant) model is run;
the end products are the annotated VCFs, the sibling sharing segments and the CNV / repeat / methylation comparisons.

### `ANNOTATE_AF` — `modules/af_annotate/main.nf`
`bin/annotate_af.py`, once for the whole cohort, right after `SNPEFF`. Writes population allele frequency into
every child's snpEff VCF: `INFO/gnomAD_AF_NFE` (`params.gnomad_chm13_vcf`, Ensembl's projection of gnomAD onto
CHM13) and `INFO/KG_AF_EUR` (1000 Genomes EUR re-called on CHM13, `kg_af` subworkflow; off unless
`kg_chm13_vcfs` / `kg_af_vcf` is set). The two fill each other's gaps — the gnomAD projection misses sites where
CHM13 carries the minor allele. Annotation only, nothing filtered; a `##af_sources` header records which
sources were used.
- **in**: snpEff VCF manifest + VCFs, 1000G sites VCF (or `NO_FILE`)
- **out** (`06_variants/snpeff_af`): `<child>.snpeff.af.vcf.gz`

### `ORIENT_PHASE` — `modules/haplotype/main.nf`
Per family, `bin/orient_family_phase.py` on `WHATSHAP_PHASE_FAMILY`'s family VCF. hap-ibd reads GT order as two
chromosome-long haplotypes and ignores phase blocks, so each child's WhatsHap blocks are flipped to
**paternal|maternal**: a block's HP1 is assigned to a parent by voting over its Mendelian-informative het SNVs
(parents' genotypes from the same VCF; ≥ 3 votes, ≥ 80% agreement). Children only, biallelic PASS SNVs; a site is
dropped if any child's het is unoriented or below `hapibd_min_gq`. chrX/Y/M excluded (`hap_exclude_chroms`).
- **out** (`07_haplotypes/oriented`): `<family>.oriented.vcf.gz`, `<family>.orient_blocks.tsv`

### `HAP_IBD` — `modules/haplotype/main.nf`
Per family, hap-ibd (`params.hapibd_jar`, `hapibd` label, Java container) on the oriented VCF with a uniform
1 cM/Mb map generated from the reference `.fai`, `min-mac=1` (a variant only one child carries is still evidence
their haplotypes differ). Segment-length and gap settings: `hapibd_min_seed` / `_min_extend` / `_min_output` /
`_max_gap`. Haplotype 1 = paternal, 2 = maternal for every child.
- **out** (`07_haplotypes/hapibd`): `<family>.hapibd.ibd.gz`, `.log`

### `IBD_SEGMENTS` — `modules/haplotype/main.nf`
Cohort-wide, `bin/hapibd_segments.py`. Per mild × severe pair: mild hap 1 IBD with severe hap 1 → paternal copy
shared, hap 2 with hap 2 → maternal copy shared; the rest of each chromosome (first to last marker) is
not_shared. A not_shared stretch shorter than `hapibd_bridge_cm` with shared on both sides is bridged to shared
(a real one needs two nearby crossovers; hap-ibd segments are routinely broken by genotype / switch errors).
Cross-haplotype segments (hap 1 with hap 2) are counted, not used.
- **out** (`07_haplotypes`): `cohort.sib_parent_segments.tsv` (per parent: shared / not_shared, bridged, n_sites,
  n_mismatch), `cohort.sib_ibd_segments.tsv` (IBD0/1/2), `cohort.sib_ibd_summary.tsv` (genome fractions — expect
  ~25/50/25% — mismatch rate inside shared segments, cross-haplotype and bridged counts)

### `METHYLATION_COMPARE` — `modules/methylation_compare/main.nf` (`--run_methylation` + trio)
`bin/methylation_compare.py`. Mild vs severe promoter 5mC (−2 kb/+500 bp around each TSS, merged per gene), per
pair: combined (all reads) and per parental haplotype. modkit's HP1/HP2 are longphase labels, consistent only
within a phase block, so each longphase_sv block's HP1 is assigned to a parent by voting over informative het
SNVs (one parent carries the ALT, the other hom-ref). Each side is labelled `shared` / `not_shared` from
`cohort.sib_parent_segments.tsv` — a difference on a not_shared haplotype may be sequence-driven (cis); on a
shared one it is epigenetic on identical sequence. Flag at |Δ| ≥ `meth_min_delta` (0.2).
- **in**: modkit `hp1/hp2/combined` bedMethyl + longphase_sv VCF per child, parents' Clair3 VCFs, sib segments, GTF
- **out** (`10_modifier_analysis/methylation`): `cohort.promoter_methylation.tsv` (+ `.flagged.tsv`),
  `cohort.methylation_phase_blocks.tsv` (block → parent), `cohort.methylation_imprinting_check.tsv`
  (SNRPN/PEG3/MEST/… should be methylated on the maternal side, H19/MEG3 paternal — `consistent=no` across the
  board means the HP → parent mapping is inverted)

### `MODIFIER_VARIANTS` — `modules/modifier/main.nf` (trio + `--run_snpeff`, `--run_modifier_analysis`)
`bin/modifier_variants.py`, one task for the whole cohort. Children are grouped into families (`<father>_<mother>`
in `family.json`); in each *discordant* family (≥ 1 mild + ≥ 1 severe) the mild child is tested against its
severe sib(s), from the annotated master VCFs. Severe children of *concordant* families (severe only) are pooled
as extra exclusion controls: every rule below that excludes on the severe sib also excludes on any concordant
severe child, but they are never written to any output (the summary counts what they excluded). With no
concordant family in `family.json` that step is skipped. Plain genotype rules, no statistical model, no AF filter. A variant is *genic*
when snpEff puts it in a protein-coding gene (CDS-affecting, splice, exon, intron or UTR term); otherwise it is
*regulatory* and linked to the nearest protein-coding gene on each side (GTF, strand ignored). Small variants
match exactly; SVs by type + POS within `modifier_sv_match_dist` (500 bp). A gene passes as:
- `homozygous` — mild homozygous for a genic variant the severe sib is not homozygous for (het / absent is fine)
- `compound_het` — ≥ 1 het on the paternal and ≥ 1 on the maternal copy (`INFO/ORIGIN`), ≥ 1 of them absent from
  the severe sib; no-ORIGIN hets fall back to opposite haplotypes of one phase block. SVs count like any variant
- `regulatory` — a non-genic variant the severe sib doesn't carry with the same zygosity (both flanking genes)
- **in**: annotated master VCFs of the family, `cohort.sib_ibd_segments.tsv`, GTF (`cnv_repeat_gtf` / `snpeff_gtf`)
- **out** (`10_modifier_analysis/variants/families/<family>`): `<mild>.passing_genes.tsv`,
  `<mild>.mild.vcf.gz` and `<severe>.severe.vcf.gz` (original records + `INFO/MOD_*`, `INFO/IBD` /
  `IBD_PAT` / `IBD_MAT` at the variant, AF kept), `<mild>.modifier_variants.tsv.gz` (gene × variant),
  `<mild>.gene_variants.tsv.gz` (every mild **and** severe variant in each passing gene, with per-sib zygosity
  and `carriers` = `mild_only` / `severe_only` / `shared` / `shared_diff`), `<mild>.modifier_summary.tsv`
  (counts per step)

### `GENE_RECURRENCE` — `modules/modifier/main.nf`
`bin/gene_recurrence.py`. Counts the mild children each gene passes in (total and per rule). For genes in
≥ `modifier_min_mild` (2) mild children, every variant in the gene carried by those mild children or their severe
sib(s) goes into two deduplicated AlphaGenome / AlphaMissense input VCFs, `mild` and `severe` (a variant carried on
both sides is in both), with `INFO/MILD_SAMPLES` / `SEVERE_SAMPLES` (sample:zygosity), `SAMPLES` (this side),
`CARRIERS` (family:class), `GENE`, `SOURCE_MODEL` (rule supported, or `context`) and `REASON`. Symbolic SVs left
out; optional `modifier_candidate_max_af`, off by default.
- **out** (`10_modifier_analysis/variants`): `cohort.passing_genes.tsv`, `cohort.gene_recurrence.tsv`,
  `cohort.candidate_gene_variants.tsv.gz`, `candidate_variants.{mild,severe}.chm13.vcf`

### AlphaGenome / AlphaMissense — `modules/alphagenome/main.nf`, `modules/alphamissense/main.nf` (`modifier_annotation`, `subworkflows/annotation.nf`)
On each of `candidate_variants.{mild,severe}.chm13.vcf` (one chain per side), lifted CHM13 → hg38: `LIFTOVER_ALPHAGENOME_VARIANTS` → `ALPHAGENOME` (chunked,
`--run_alphagenome`) → `MERGE_ALPHAGENOME_CHUNKS`, and/or `ALPHAGENOME_API` (`--run_alphagenome_api`);
`LIFTOVER_ALPHAMISSENSE_VARIANTS` → `ALPHAMISSENSE` (when `alphamissense_tsv` is set). Per-variant outputs only,
no report; every row carries `MILD_SAMPLES` / `SEVERE_SAMPLES` / `CARRIERS`.
- **out**: `10_modifier_analysis/variants/{alphagenome,alphagenome_api,alphamissense}/{mild,severe}`

---

## Included but not invoked

**Kept, not called.**
- `MILD_HAPLOTYPE_VARIANTS` — `modules/haplotype/main.nf` (`bin/mild_haplotype_variants.py`): every mild-child
  variant on a parental copy the severe sib didn't inherit (by `INFO/ORIGIN` + `cohort.sib_parent_segments.tsv`).
  For the haplotype-aware step to be added after the modifier search.
- `VARIANT_REPORT` — `modules/variant_report/main.nf` (`bin/build_variant_report.py`): the old combined per-variant
  report over the recessive / dominant TSVs. Not compatible with the modifier search outputs.

**Illumina.**

These modules exist (`subworkflows/illumina.nf`) but are no longer `include`d in `main.nf` — present in case a future Illumina-only branch is wired in, not part of any current run:

- `ATRIA` — `modules/atria/main.nf` (lost its only caller when `hero_correction` was replaced by
  `herro_correction`, which doesn't touch Illumina reads at all)
- `BWA_INDEX` / `BWA_MEM` — `modules/bwa/main.nf`
- `FILTER_BAM` / `INDEX_BAM` — `modules/samtools/main.nf` (via `illumina.nf`'s `filter_bam` workflow)
- `DEEPVARIANT_CALL_VARIANTS` / `SPLIT_DV_VCF` — `modules/deepvariant/main.nf`
- `BAM_TO_BW` / `BAM_COMPARE` / `BW_COMPARE` / `MERGE_BW` — `modules/deeptools/main.nf`

`WHATSHAP_PHASE` / `WHATSHAP_HAPLOTAG` (`modules/whatshap/main.nf`) are similarly unused — only
`WHATSHAP_PHASE_FAMILY` is called.
