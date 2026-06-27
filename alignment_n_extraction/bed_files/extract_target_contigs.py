#!/usr/bin/env python3
"""
==============================================================================
 extract_target_contigs.py
==============================================================================
 Reads {species}_ASCT_mRNA.bed files to identify which contigs/scaffolds
 contain ASCT gene copies, then extracts those contigs (with configurable
 flanking) from the corresponding genome FASTAs.

 Pairing uses species_map from config.yaml:
   BED: Rphi_ASCT_mRNA.bed  →  species code "Rphi"
   config.yaml: ASM2657151v2 → Rphi
   FASTA: any file containing "ASM2657151v2" in its name

 Usage:
   python extract_target_contigs.py --config config.yaml
   python extract_target_contigs.py --config config.yaml --catalog-only
   python extract_target_contigs.py --config config.yaml --full-contig
==============================================================================
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import yaml

FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


# ─────────────────────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────────────────────

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


# ─────────────────────────────────────────────────────────────────────────────
# PARSE _ASCT_mRNA.bed
# ─────────────────────────────────────────────────────────────────────────────

def find_asct_beds(bed_dir: str) -> dict:
    """Find {species}_ASCT_mRNA.bed files. Returns {species_code: path}."""
    found = {}
    for f in Path(bed_dir).iterdir():
        if f.name.endswith("_ASCT_mRNA.bed"):
            species = f.name.replace("_ASCT_mRNA.bed", "")
            found[species] = str(f)
    return dict(sorted(found.items()))


def parse_asct_bed(bed_path: str) -> list:
    """
    Parse one _ASCT_mRNA.bed file.
    Format: scaffold  start  end  name  cum_bitscore  strand  qstart  qend  qcov  n_exons
    """
    copies = []
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                parts = line.split()
            if len(parts) < 4:
                continue
            copy = {
                "scaffold": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
                "name": parts[3],
                "bitscore": float(parts[4]) if len(parts) > 4 else 0.0,
                "strand": parts[5] if len(parts) > 5 else "+",
                "qstart": int(parts[6]) if len(parts) > 6 else 0,
                "qend": int(parts[7]) if len(parts) > 7 else 0,
                "qcov": float(parts[8]) if len(parts) > 8 else 0.0,
                "n_exons": int(parts[9]) if len(parts) > 9 else 0,
            }
            copies.append(copy)
    return copies


# ─────────────────────────────────────────────────────────────────────────────
# FASTA PAIRING & EXTRACTION
# ─────────────────────────────────────────────────────────────────────────────

def find_fasta_for_species(species_code: str, fasta_dir: str, code_to_assembly: dict):
    """Find FASTA whose filename contains the assembly substring for this species."""
    assembly_substr = code_to_assembly.get(species_code)
    if assembly_substr is None:
        return None
    for f in Path(fasta_dir).iterdir():
        if f.is_file() and f.name.endswith(FASTA_EXTENSIONS):
            if assembly_substr in f.name:
                return str(f)
    return None


def index_fasta_lengths(fasta_path: str) -> dict:
    """Get {seq_id: length} for all sequences in a FASTA."""
    lengths = {}
    current = None
    length = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = length
                current = line[1:].strip().split()[0]
                length = 0
            else:
                length += len(line.strip())
    if current is not None:
        lengths[current] = length
    return lengths


def extract_full_contigs(fasta_path: str, contig_names: set, out_path: str) -> dict:
    """Extract whole contigs from FASTA. Returns {contig: length}."""
    extracted = {}
    writing = False
    current = None
    current_len = 0

    with open(fasta_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                if writing and current:
                    extracted[current] = current_len
                seq_id = line[1:].strip().split()[0]
                if seq_id in contig_names:
                    writing = True
                    current = seq_id
                    current_len = 0
                    fout.write(line)
                else:
                    writing = False
                    current = None
            elif writing:
                fout.write(line)
                current_len += len(line.strip())
        if writing and current:
            extracted[current] = current_len
    return extracted


def extract_flanked_regions(fasta_path: str, regions: list, out_path: str):
    """
    Extract ±flank regions via samtools faidx.
    regions: list of (contig, start, end, label)
    Falls back to full contig if samtools unavailable.
    """
    # Index if needed
    fai = fasta_path + ".fai"
    if not os.path.exists(fai):
        try:
            subprocess.run(["samtools", "faidx", fasta_path],
                           check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("      [WARN] samtools not available, extracting full contigs")
            contig_names = {r[0] for r in regions}
            extract_full_contigs(fasta_path, contig_names, out_path)
            return

    with open(out_path, "w") as fout:
        for contig, start, end, label in regions:
            spec = f"{contig}:{start}-{end}"
            try:
                result = subprocess.run(
                    ["samtools", "faidx", fasta_path, spec],
                    capture_output=True, text=True, check=True
                )
                lines = result.stdout.strip().split("\n")
                fout.write(f">{contig}:{start}-{end} {label}\n")
                seq = "".join(lines[1:])
                for i in range(0, len(seq), 80):
                    fout.write(seq[i:i+80] + "\n")
                print(f"      {label}: {contig}:{start:,}-{end:,} ({len(seq):,} bp)")
            except subprocess.CalledProcessError:
                print(f"      [WARN] Failed: {spec}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def merge_windows(windows):
    """Merge overlapping (start, end) intervals."""
    if not windows:
        return []
    windows = sorted(windows)
    merged = [list(windows[0])]
    for start, end in windows[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def process_species(species, bed_path, fasta_path, flank, out_dir, mode, contig_lengths=None):
    """Parse BED, group by contig, extract ±flank around EACH gene copy."""
    print(f"\n  [{species}]")
    copies = parse_asct_bed(bed_path)
    print(f"    {len(copies)} copies from {Path(bed_path).name}")

    by_contig = defaultdict(list)
    for c in copies:
        by_contig[c["scaffold"]].append(c)

    for contig, cc in sorted(by_contig.items()):
        names = [c["name"] for c in cc]
        lo = min(c["start"] for c in cc)
        hi = max(c["end"] for c in cc)
        print(f"    {contig}: {len(cc)} copies ({', '.join(names)}) "
              f"span {lo:,}-{hi:,}")

    if fasta_path:
        if mode == "contig":
            out_fa = out_dir / f"{species}.target_contigs.fa"
            extracted = extract_full_contigs(fasta_path, set(by_contig.keys()), str(out_fa))
            print(f"    → {out_fa.name} ({len(extracted)} contigs)")
        else:
            lengths = contig_lengths or index_fasta_lengths(fasta_path)
            regions = []
            for contig, cc in sorted(by_contig.items()):
                contig_len = lengths.get(contig, 10**9)
                # Create one window per copy, then merge overlapping
                windows = []
                for c in cc:
                    w_start = max(0, c["start"] - flank)
                    w_end = min(contig_len, c["end"] + flank)
                    windows.append((w_start, w_end))
                merged = merge_windows(windows)
                for i, (r_start, r_end) in enumerate(merged):
                    if len(merged) == 1:
                        label = f"{species}_{contig}"
                    else:
                        label = f"{species}_{contig}_r{i+1}"
                    regions.append((contig, r_start, r_end, label))
                    span = r_end - r_start
                    # Report what copies fall in this window
                    copies_in = [c["name"] for c in cc
                                 if c["start"] >= r_start and c["end"] <= r_end]
                    print(f"    → region {label}: {contig}:{r_start:,}-{r_end:,} "
                          f"({span:,} bp, {len(copies_in)} copies)")

            out_fa = out_dir / f"{species}.target_regions.fa"
            extract_flanked_regions(fasta_path, regions, str(out_fa))

    return {
        "species": species,
        "copies": copies,
        "contigs": dict(by_contig),
        "n_copies": len(copies),
        "n_contigs": len(by_contig),
    }


def main():
    parser = argparse.ArgumentParser(description="Extract ASCT target contigs")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--catalog-only", action="store_true",
                        help="Just catalog copies, skip FASTA extraction")
    parser.add_argument("--full-contig", action="store_true",
                        help="Override config: extract full contigs")
    args = parser.parse_args()

    cfg = load_config(args.config)
    bed_dir = cfg["bed_dir"]
    fasta_dir = cfg.get("fasta_dir")
    out_dir = Path(cfg["output_dir"]) / "00_extracted"
    out_dir.mkdir(parents=True, exist_ok=True)
    flank = cfg.get("gene_flank_bp", cfg.get("flank_bp", 50_000))
    mode = "contig" if args.full_contig else cfg.get("extract_mode", "flank")

    # Build reverse map: species_code → assembly_substring
    species_map = cfg.get("species_map", {})
    code_to_assembly = {v: k for k, v in species_map.items()}

    print("=" * 70)
    print("EXTRACT TARGET CONTIGS")
    print("=" * 70)

    # Find BED files
    asct_beds = find_asct_beds(bed_dir)
    print(f"\n  Found {len(asct_beds)} species BED files")

    # Pair with FASTAs
    fasta_map = {}
    if fasta_dir and not args.catalog_only:
        for sp in asct_beds:
            fasta = find_fasta_for_species(sp, fasta_dir, code_to_assembly)
            if fasta:
                fasta_map[sp] = fasta
                print(f"    {sp} → {Path(fasta).name}")
            else:
                print(f"    {sp} → NOT FOUND (key: {code_to_assembly.get(sp, '???')})")

    # Process
    results = []
    for sp, bed_path in asct_beds.items():
        fasta = fasta_map.get(sp) if not args.catalog_only else None
        result = process_species(sp, bed_path, fasta, flank, out_dir, mode)
        results.append(result)

    # ── Write copy catalog ──
    catalog_path = out_dir / "copy_catalog.tsv"
    with open(catalog_path, "w") as f:
        f.write("species\tcopy_name\tscaffold\tstart\tend\tstrand\t"
                "bitscore\tqcov\tn_exons\n")
        for r in results:
            for c in r["copies"]:
                f.write(f"{r['species']}\t{c['name']}\t{c['scaffold']}\t"
                        f"{c['start']}\t{c['end']}\t{c['strand']}\t"
                        f"{c['bitscore']}\t{c['qcov']}\t{c['n_exons']}\n")
    print(f"\n  Catalog: {catalog_path}")

    # ── Write copy counts ──
    counts_path = out_dir / "copy_counts.tsv"
    with open(counts_path, "w") as f:
        f.write("species\tn_copies\tn_contigs\tcopy_names\tscaffolds\n")
        for r in sorted(results, key=lambda x: -x["n_copies"]):
            names = ", ".join(c["name"] for c in r["copies"])
            scaffolds = ", ".join(sorted(r["contigs"].keys()))
            f.write(f"{r['species']}\t{r['n_copies']}\t{r['n_contigs']}\t"
                    f"{names}\t{scaffolds}\n")
    print(f"  Counts:  {counts_path}")

    # ── Print summary ──
    max_copies = max(r["n_copies"] for r in results)
    print(f"\n  {'Species':<8} {'Copies':>6}  {'':>1}")
    print(f"  {'─'*8} {'─'*6}  {'─'*30}")
    for r in sorted(results, key=lambda x: -x["n_copies"]):
        bar = "█" * r["n_copies"] + "░" * (max_copies - r["n_copies"])
        print(f"  {r['species']:<8} {r['n_copies']:>6}  {bar}")

    print(f"\n  Next: python 00_preprocess.py --config config.yaml")


if __name__ == "__main__":
    main()
