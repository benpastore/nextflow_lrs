#!/usr/bin/env python3
"""
Write population allele frequencies into every child's consolidated+snpEff
VCF, once, right after snpEff -- so every later step reads AF from
INFO instead of re-scanning the population files itself.

  INFO/gnomAD_AF_NFE   gnomAD non-Finnish European AF (-gnomad_vcf, field AF_nfe):
                       Ensembl's projection of gnomAD onto CHM13. Large, but it
                       misses sites where CHM13 itself carries the minor allele.
  INFO/KG_AF_EUR       1000 Genomes EUR AF (-kg_af_vcf, field AF): 1KGP re-called
                       natively on CHM13 (modules/kg_af), which fills those gaps.

A field is written only where the source has the exact allele
(chrom/pos/ref/alt); absent = not in that source (or the source wasn't
given -- see the ##af_sources header line). Annotation only: nothing is
filtered here.

All samples' sites are collected first, so each population file is
streamed exactly once for the whole cohort (GnomadNFEAnnotator, shared
with the modifier scripts).

    annotate_af.py -manifest vcf_manifest.tsv -gnomad_vcf gnomad.vcf.gz \\
        -kg_af_vcf 1kgp.EUR.sites.vcf.gz

Writes <sampleID>.snpeff.af.vcf.gz next to the inputs' names.
Written for Python 3.6 (the rnaseq container).
"""
import argparse
import gzip
import os

from _modifier_common import GnomadNFEAnnotator, load_vcf_manifest, logger, opener


def present(path):
    return bool(path) and os.path.basename(path) != "NO_FILE" and os.path.exists(path)


def collect_sites(paths):
    """Every (chrom, pos, ref, alt) across the cohort's VCFs, deduplicated."""
    sites = set()
    for path in paths:
        with opener(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.split("\t", 5)
                if len(c) < 5:
                    continue
                sites.add((c[0], c[1], c[3], c[4]))
    logger.info(f"{len(sites):,} unique sites across {len(paths)} VCF(s)")
    return sites


def fmt_af(v):
    return f"{v:.6g}"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-manifest", required=True, help="sampleID<TAB>snpEff VCF, VCFs staged alongside")
    p.add_argument("-gnomad_vcf", default="", help="gnomAD CHM13 VCF (INFO/AF_nfe); empty = skip")
    p.add_argument("-kg_af_vcf", default="", help="1000 Genomes CHM13 EUR sites VCF (INFO/AF); empty = skip")
    args = p.parse_args()

    manifest = load_vcf_manifest(args.manifest)
    sites = collect_sites(list(manifest.values()))

    sources = []   # (INFO key, description, {variant_id: AF})
    for key, path, field, label, desc in (
            ("gnomAD_AF_NFE", args.gnomad_vcf, "AF_nfe", "gnomAD",
             "gnomAD non-Finnish European allele frequency (CHM13 projection)"),
            ("KG_AF_EUR", args.kg_af_vcf, "AF", "1000 Genomes EUR",
             "1000 Genomes EUR allele frequency (1KGP re-called on CHM13)")):
        if present(path):
            afs = GnomadNFEAnnotator(path, field=field, label=label).annotate_batch(list(sites))
            sources.append((key, desc, afs, os.path.basename(path)))
        else:
            logger.warning(f"{label}: no file given -- INFO/{key} not written")
            sources.append((key, desc, None, "none"))

    header = [f'##INFO=<ID={key},Number=1,Type=Float,Description="{desc}">\n'
              for key, desc, afs, _ in sources if afs is not None]
    header.append("##af_sources=" + ",".join(f"{key}:{name}" for key, _, _, name in sources) + "\n")

    for sample, path in sorted(manifest.items()):
        out_path = f"{sample}.snpeff.af.vcf.gz"
        n, hits = 0, {key: 0 for key, _, _, _ in sources}
        with opener(path) as fin, gzip.open(out_path, "wt", compresslevel=4) as fout:
            for line in fin:
                if line.startswith("##"):
                    fout.write(line)
                    continue
                if line.startswith("#CHROM"):
                    fout.writelines(header)
                    fout.write(line)
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 8:
                    fout.write(line)
                    continue
                n += 1
                vid = f"{c[0]}:{c[1]}:{c[3]}>{c[4]}"
                add = []
                for key, _, afs, _ in sources:
                    if afs is not None and afs.get(vid) is not None:
                        add.append(f"{key}={fmt_af(afs[vid])}")
                        hits[key] += 1
                if add:
                    c[7] = ";".join(add) if c[7] in (".", "") else c[7] + ";" + ";".join(add)
                fout.write("\t".join(c) + "\n")
        logger.info(f"{sample}: {n:,} records -> {out_path}; "
                    + ", ".join(f"{k} on {v:,}" for k, v in hits.items()))


if __name__ == "__main__":
    main()
