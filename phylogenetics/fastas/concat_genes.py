#!/usr/bin/env python3
"""
Concatenate gene copy FASTA files per species in genomic order and orientation.

Uses BED files to determine:
  - Order: sorted by scaffold (alphabetical), then start position
  - Orientation: minus-strand sequences are reverse-complemented

Output: one FASTA per species with a single concatenated sequence.
"""

import glob
import os
import sys
from pathlib import Path

COMP = str.maketrans("ACGTacgtRYSWKMBDHVryswkmbdhv",
                      "TGCAtgcaYRWSMKVHDByrwsmkvhdb")

def reverse_complement(seq):
    return seq.translate(COMP)[::-1]

def read_fasta(path):
    """Return the sequence from a single-record FASTA file."""
    lines = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            lines.append(line.strip())
    return "".join(lines)

def parse_bed(path):
    """Return list of dicts with scaffold, start, name, strand."""
    entries = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            entries.append({
                "scaffold": cols[0],
                "start": int(cols[1]),
                "end": int(cols[2]),
                "name": cols[3],
                "strand": cols[5],
            })
    return entries

def main():
    workdir = Path(__file__).parent

    species_with_fastas = set()
    for fasta in workdir.glob("*_ASCT*_gene.fasta"):
        species_with_fastas.add(fasta.name.split("_ASCT")[0])

    for species in sorted(species_with_fastas):
        bed_path = workdir / f"{species}_ASCT_mRNA.bed"
        if not bed_path.exists():
            print(f"WARNING: no BED file for {species}, skipping", file=sys.stderr)
            continue

        entries = parse_bed(bed_path)
        entries.sort(key=lambda e: (e["scaffold"], e["start"]))

        print(f"{species}: {len(entries)} copies, order = ", end="")
        for e in entries:
            print(f"{e['name']}({e['strand']}) ", end="")
        print()

        concat_seq = []
        for e in entries:
            fasta_path = workdir / f"{e['name']}_gene.fasta"
            if not fasta_path.exists():
                print(f"  WARNING: {fasta_path.name} not found, skipping", file=sys.stderr)
                continue
            seq = read_fasta(fasta_path)
            if e["strand"] == "-":
                print(f"  Reverse-complementing {e['name']}")
                seq = reverse_complement(seq)
            concat_seq.append(seq)

        out_path = workdir / f"{species}_ASCT_concat.fasta"
        with open(out_path, "w") as f:
            full_seq = "".join(concat_seq)
            f.write(f">{species}_ASCT_concat\n")
            for i in range(0, len(full_seq), 80):
                f.write(full_seq[i:i+80] + "\n")

        print(f"  -> {out_path.name} ({len(''.join(concat_seq))} bp)")

if __name__ == "__main__":
    main()
