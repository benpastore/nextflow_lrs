#!/usr/bin/env python3
"""
Consolidate per-caller VCFs (alignment-based and assembly-based) for a single
sample into one master VCF. Records, for every variant site, which callers
detected it, how many callers detected it, and keeps each caller's own
genotype (including phase / haplotype block, via GT + PS) as a separate
sample column so haplotype information is not lost in the merge.

Small variants (SNP/indel) are matched by exact (chrom, pos, ref, alt).
Structural variants (INFO/SVTYPE present) are matched by (chrom, svtype)
plus proximity clustering, since breakpoint coordinates rarely agree exactly
across callers. This is a lightweight positional clustering, not a full
breakend-aware SV merge (see --sv-merge-dist).

Two ways to run this:

  1. Explicit, single sample (what the Nextflow module uses):
       consolidate_variants.py --sample S1 \\
           --caller clair3_longphase:alignment:S1.longphase_sv.vcf.gz \\
           --caller sniffles:alignment:S1.sniffles.vcf.gz \\
           --caller dipcall:assembly:S1.dip.vcf.gz \\
           --caller hapdiff:assembly:S1.hapdiff.phased.vcf.gz \\
           --output S1.master.vcf

  2. Auto-discovery over a results/variants-style directory tree, for
     running standalone outside the pipeline (every sample found gets its
     own output file):
       consolidate_variants.py --input-dir /path/to/results --outdir out/

     This recursively scans --input-dir for filenames matching the
     naming convention this pipeline's variant-calling modules use
     (*.longphase_sv.vcf.gz, *.longphase.vcf.gz, *.clair3.vcf.gz,
     *.sniffles.vcf.gz, *.dip.vcf.gz, *.hapdiff.phased.vcf.gz), infers
     the sample ID from each filename, and groups files by sample ID.
     For the small-variant slot it prefers the most-refined file found
     per sample (longphase_sv > longphase > clair3) rather than using
     more than one, since those all trace back to the same clair3 calls.
"""

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict


def opener(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_caller_arg(arg):
    parts = arg.split(":", 2)
    if len(parts) != 3:
        sys.exit(f"--caller must be NAME:TYPE:PATH, got: {arg}")
    name, ctype, path = parts
    if ctype not in ("alignment", "assembly"):
        sys.exit(f"--caller TYPE must be 'alignment' or 'assembly', got: {ctype} ({arg})")
    return name, ctype, path


SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
END_RE = re.compile(r"(?:^|;)END=(-?[0-9]+)")


def parse_vcf(path, caller):
    """Yield one dict per ALT allele in the VCF."""
    records = []
    with opener(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, _vid, ref, alt, qual, filt, info = f[:8]
            fmt = f[8] if len(f) > 8 else ""
            sample = f[9] if len(f) > 9 else ""
            gt, ps = "", ""
            if fmt and sample:
                keys = fmt.split(":")
                vals = sample.split(":")
                d = dict(zip(keys, vals))
                gt = d.get("GT", "")
                ps = d.get("PS", "")

            m = SVTYPE_RE.search(info)
            svtype = m.group(1) if m else None
            m = END_RE.search(info)
            end = int(m.group(1)) if m else None

            for a in alt.split(","):
                if a in (".", ""):
                    continue
                records.append(
                    {
                        "chrom": chrom,
                        "pos": int(pos),
                        "ref": ref,
                        "alt": a,
                        "qual": qual,
                        "filter": filt,
                        "caller": caller,
                        "gt": gt,
                        "ps": ps,
                        "svtype": svtype,
                        "end": end,
                    }
                )
    return records


def cluster_svs(sv_records, merge_dist):
    """Greedy proximity clustering of SV records: same chrom + svtype,
    start positions within merge_dist of the running cluster start."""
    sv_records = sorted(
        sv_records, key=lambda r: (r["chrom"], r["svtype"] or "", r["pos"])
    )
    clusters = []
    current = None
    for rec in sv_records:
        if (
            current is not None
            and rec["chrom"] == current["chrom"]
            and (rec["svtype"] or "") == (current["svtype"] or "")
            and rec["pos"] - current["cluster_pos"] <= merge_dist
        ):
            current["members"].append(rec)
            current["cluster_pos"] = rec["pos"]
        else:
            if current is not None:
                clusters.append(current)
            current = {
                "chrom": rec["chrom"],
                "svtype": rec["svtype"],
                "pos": rec["pos"],
                "cluster_pos": rec["pos"],
                "members": [rec],
            }
    if current is not None:
        clusters.append(current)
    return clusters


def build_sites(all_records, merge_dist):
    small = {}
    sv = []
    for rec in all_records:
        if rec["svtype"]:
            sv.append(rec)
        else:
            key = (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
            small.setdefault(key, []).append(rec)

    sites = []
    for (chrom, pos, ref, alt), members in small.items():
        sites.append(
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "svtype": None,
                "end": None,
                "members": members,
            }
        )

    for cluster in cluster_svs(sv, merge_dist):
        members = cluster["members"]
        ends = [m["end"] for m in members if m["end"] is not None]
        sites.append(
            {
                "chrom": cluster["chrom"],
                "pos": min(m["pos"] for m in members),
                "ref": "N",
                "alt": f"<{cluster['svtype'] or 'SV'}>",
                "svtype": cluster["svtype"],
                "end": max(ends) if ends else None,
                "members": members,
            }
        )

    sites.sort(key=lambda s: (s["chrom"], s["pos"]))
    return sites


# (filename-suffix regex, caller name, caller type) -- checked in order, first
# match per sample wins. Grouped so the small-variant family (all downstream
# of clair3) only ever contributes one caller per sample.
SMALL_VARIANT_PATTERNS = [
    (re.compile(r"\.longphase_sv\.vcf\.gz$"), "clair3_longphase", "alignment"),
    (re.compile(r"\.longphase\.vcf\.gz$"), "clair3_longphase", "alignment"),
    (re.compile(r"\.clair3\.vcf\.gz$"), "clair3", "alignment"),
]
STANDALONE_PATTERNS = [
    (re.compile(r"\.sniffles\.vcf\.gz$"), "sniffles", "alignment"),
    (re.compile(r"\.dip\.vcf\.gz$"), "dipcall", "assembly"),
    (re.compile(r"\.hapdiff\.phased\.vcf\.gz$"), "hapdiff", "assembly"),
]
CALLER_PRIORITY = ["clair3_longphase", "clair3", "sniffles", "dipcall", "hapdiff"]


def discover_samples(input_dir):
    """Walk input_dir and group recognized per-caller VCFs by sample ID.

    Returns {sample_id: {caller_name: (type, path)}}.
    """
    all_files = []
    for root, _, files in os.walk(input_dir):
        for fname in files:
            all_files.append((fname, os.path.join(root, fname)))

    samples = defaultdict(dict)

    small_variant_done = set()
    for pattern, caller_name, ctype in SMALL_VARIANT_PATTERNS:
        for fname, path in all_files:
            m = pattern.search(fname)
            if not m or fname[: m.start()] in small_variant_done:
                continue
            sample = fname[: m.start()]
            samples[sample][caller_name] = (ctype, path)
            small_variant_done.add(sample)

    for pattern, caller_name, ctype in STANDALONE_PATTERNS:
        for fname, path in all_files:
            m = pattern.search(fname)
            if not m:
                continue
            sample = fname[: m.start()]
            samples[sample][caller_name] = (ctype, path)

    return samples


def write_master_vcf(output_path, sample, callers, caller_type, sites):
    out = opener(output_path, "wt") if output_path.endswith(".gz") else open(
        output_path, "wt"
    )
    with out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##source=consolidate_variants.py (sample={sample})\n")
        out.write(
            f"##consolidate_variants_callers={','.join(f'{c}:{caller_type[c]}' for c in callers)}\n"
        )
        out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant (SV clusters: max END across supporting callers)">\n')
        out.write('##INFO=<ID=NUM_CALLERS,Number=1,Type=Integer,Description="Number of variant callers that detected this variant">\n')
        out.write('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Names of variant callers that detected this variant">\n')
        out.write('##INFO=<ID=NUM_ALIGNMENT_CALLERS,Number=1,Type=Integer,Description="Number of alignment-based callers that detected this variant">\n')
        out.write('##INFO=<ID=NUM_ASSEMBLY_CALLERS,Number=1,Type=Integer,Description="Number of assembly-based callers that detected this variant">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype as reported by this caller (phased with | when the caller phased it)">\n')
        out.write('##FORMAT=<ID=PS,Number=1,Type=String,Description="Phase set ID as reported by this caller, if phased">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(callers) + "\n")

        for i, site in enumerate(sites):
            calls = {}
            for m in site["members"]:
                # keep the first call seen per caller at this site
                calls.setdefault(m["caller"], m)

            detected = sorted(calls.keys(), key=callers.index)
            n_align = sum(1 for c in detected if caller_type[c] == "alignment")
            n_asm = sum(1 for c in detected if caller_type[c] == "assembly")

            if site["svtype"]:
                vid = f"{site['chrom']}_{site['pos']}_{site['svtype']}_{i}"
            else:
                vid = f"{site['chrom']}_{site['pos']}_{site['ref']}_{site['alt']}"

            info_fields = [
                f"NUM_CALLERS={len(detected)}",
                f"CALLERS={','.join(detected)}",
                f"NUM_ALIGNMENT_CALLERS={n_align}",
                f"NUM_ASSEMBLY_CALLERS={n_asm}",
            ]
            if site["svtype"]:
                info_fields.insert(0, f"SVTYPE={site['svtype']}")
                if site["end"] is not None:
                    info_fields.insert(1, f"END={site['end']}")

            row = [
                site["chrom"],
                str(site["pos"]),
                vid,
                site["ref"],
                site["alt"],
                ".",
                ".",
                ";".join(info_fields),
                "GT:PS",
            ]
            for c in callers:
                m = calls.get(c)
                if m is None:
                    row.append("./.:.")
                else:
                    gt = m["gt"] or "./."
                    ps = m["ps"] or "."
                    row.append(f"{gt}:{ps}")
            out.write("\t".join(row) + "\n")

    sys.stderr.write(
        f"[consolidate_variants] {sample}: {len(sites)} sites from "
        f"{len(callers)} callers ({', '.join(callers)}) -> {output_path}\n"
    )


def run_explicit(args):
    callers = []
    caller_type = {}
    all_records = []
    for arg in args.caller:
        name, ctype, path = parse_caller_arg(arg)
        callers.append(name)
        caller_type[name] = ctype
        if path == "NO_FILE":
            continue
        try:
            recs = parse_vcf(path, name)
        except (FileNotFoundError, OSError):
            recs = []
        if not recs:
            continue
        all_records.extend(recs)

    sites = build_sites(all_records, args.sv_merge_dist)
    write_master_vcf(args.output, args.sample, callers, caller_type, sites)


def run_discovery(args):
    samples = discover_samples(args.input_dir)
    if not samples:
        sys.exit(f"No recognizable caller VCFs found under {args.input_dir}")

    os.makedirs(args.outdir, exist_ok=True)

    for sample, found in sorted(samples.items()):
        callers = sorted(found.keys(), key=CALLER_PRIORITY.index)
        caller_type = {name: found[name][0] for name in callers}

        all_records = []
        for name in callers:
            _ctype, path = found[name]
            try:
                all_records.extend(parse_vcf(path, name))
            except (FileNotFoundError, OSError) as e:
                sys.stderr.write(f"[consolidate_variants] {sample}: skipping {name} ({path}): {e}\n")

        missing = [c for c in CALLER_PRIORITY if c not in callers and not (
            c in ("clair3_longphase", "clair3") and any(x in callers for x in ("clair3_longphase", "clair3"))
        )]
        if missing:
            sys.stderr.write(f"[consolidate_variants] {sample}: no VCF found for {', '.join(missing)}\n")

        sites = build_sites(all_records, args.sv_merge_dist)
        output_path = os.path.join(args.outdir, f"{sample}.master.vcf")
        write_master_vcf(output_path, sample, callers, caller_type, sites)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--caller",
        action="append",
        help="NAME:TYPE:PATH, TYPE is 'alignment' or 'assembly'; repeatable. "
        "PATH of 'NO_FILE' (or a missing/empty file) is skipped. "
        "Single-sample mode; requires --sample and --output.",
    )
    mode.add_argument(
        "--input-dir",
        help="Recursively scan this directory for per-caller VCFs and "
        "auto-consolidate every sample found. Requires --outdir.",
    )
    p.add_argument("--sample", help="Sample ID (required with --caller)")
    p.add_argument("--output", help="Output VCF path, .vcf or .vcf.gz (required with --caller)")
    p.add_argument("--outdir", help="Output directory, one <sample>.master.vcf per sample (required with --input-dir)")
    p.add_argument(
        "--sv-merge-dist",
        type=int,
        default=500,
        help="Max bp between SV start positions (same chrom+SVTYPE) to be "
        "merged into one site (default: 500)",
    )
    args = p.parse_args()

    if args.caller:
        if not (args.sample and args.output):
            p.error("--caller mode requires --sample and --output")
        run_explicit(args)
    else:
        if not args.outdir:
            p.error("--input-dir mode requires --outdir")
        run_discovery(args)


if __name__ == "__main__":
    main()
