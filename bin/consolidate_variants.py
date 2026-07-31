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
"""
import argparse
import gzip
import re
import sys


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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True, help="Sample ID")
    p.add_argument(
        "--caller",
        action="append",
        required=True,
        help="NAME:TYPE:PATH, TYPE is 'alignment' or 'assembly'; repeatable. "
        "PATH of 'NO_FILE' (or a missing/empty file) is skipped.",
    )
    p.add_argument("--output", required=True, help="Output VCF path (.vcf or .vcf.gz)")
    p.add_argument(
        "--sv-merge-dist",
        type=int,
        default=500,
        help="Max bp between SV start positions (same chrom+SVTYPE) to be "
        "merged into one site (default: 500)",
    )
    args = p.parse_args()

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

    out = opener(args.output, "wt") if args.output.endswith(".gz") else open(
        args.output, "wt"
    )
    with out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##source=consolidate_variants.py (sample={args.sample})\n")
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
        f"[consolidate_variants] {args.sample}: {len(sites)} sites from "
        f"{len(callers)} callers ({', '.join(callers)}) -> {args.output}\n"
    )


if __name__ == "__main__":
    main()
