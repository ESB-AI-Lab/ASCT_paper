#!/usr/bin/env python3
"""
concat_exon_fastas.py

Concatenate individual exon FASTA sequences into a single gene FASTA per copy.
Exon order is determined by genomic position:
  - Plus strand (+): ascending genomic position (low to high)
  - Minus strand (-): descending genomic position (high to low)

Strand information is read from the corresponding mRNA BED file.

Usage:
    python concat_exon_fastas.py \
        --fasta-dir ./ \
        --bed-dir bed_files/ \
        --outdir gene_fastas/

Expected input:
    - Exon FASTA files: {copy}_exons.fasta
    - mRNA BED files: {species}_ASCT_mRNA.bed with strand in column 6

Output:
    - One FASTA per copy: {copy}_gene.fasta
    - Combined FASTA: all_ASCT_genes.fasta
"""

import argparse
import glob
import os
import re
import sys


def parse_mRNA_beds(bed_dir):
    """
    Parse all mRNA BED files to build a map of copy_name -> strand.
    """
    strand_map = {}
    bed_files = glob.glob(os.path.join(bed_dir, "*_ASCT_mRNA.bed"))

    for bed_file in bed_files:
        with open(bed_file, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) >= 6:
                    copy_name = fields[3]
                    strand = fields[5]
                    strand_map[copy_name] = strand

    return strand_map


def parse_position_from_header(header):
    """
    Extract genomic start position from a FASTA header.
    Handles common bedtools getfasta formats:
      >scaffold:start-end
      >scaffold:start-end(+)
      >scaffold:start-end(-)
      >copy_name/exon_N::scaffold:start-end(strand)
    """
    match = re.search(r'(\S+):(\d+)-(\d+)', header)
    if match:
        scaffold = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        return (scaffold, start, end)

    return None


def parse_exon_fasta(filepath):
    """
    Parse a multi-sequence FASTA file.
    Returns list of (header, sequence, parsed_position) tuples.
    """
    records = []
    current_header = None
    current_seq = []

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    seq = "".join(current_seq)
                    pos = parse_position_from_header(current_header)
                    records.append((current_header, seq, pos))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

    if current_header is not None:
        seq = "".join(current_seq)
        pos = parse_position_from_header(current_header)
        records.append((current_header, seq, pos))

    return records


def sort_exons(records, strand):
    """
    Sort exon records by genomic position.
    Plus strand: ascending (low to high)
    Minus strand: descending (high to low)
    """
    with_pos = [(r, r[2]) for r in records if r[2] is not None]
    without_pos = [r for r in records if r[2] is None]

    if without_pos:
        sys.stderr.write(
            "Warning: %d exon(s) could not be sorted by position. "
            "Appending at end.\n" % len(without_pos)
        )

    if strand == "-":
        with_pos.sort(key=lambda x: x[1][1], reverse=True)
    else:
        with_pos.sort(key=lambda x: x[1][1])

    sorted_records = [r[0] for r in with_pos] + without_pos
    return sorted_records


def concat_and_write(records, copy_name, strand, outdir):
    """
    Concatenate exon sequences and write a single gene FASTA.
    Header contains only the copy name.
    """
    sorted_records = sort_exons(records, strand)

    gene_seq = "".join(r[1] for r in sorted_records)

    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "%s_gene.fasta" % copy_name)

    with open(outpath, "w") as fh:
        fh.write(">%s\n" % copy_name)
        for i in range(0, len(gene_seq), 80):
            fh.write(gene_seq[i:i + 80] + "\n")

    return outpath, gene_seq, sorted_records


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Concatenate exon FASTA sequences into single gene FASTAs per "
            "copy. Respects strand orientation for exon ordering."
        )
    )
    parser.add_argument(
        "--fasta-dir", "-f", required=True,
        help="Directory containing {copy}_exons.fasta files."
    )
    parser.add_argument(
        "--bed-dir", "-b", required=True,
        help="Directory containing {species}_ASCT_mRNA.bed files for strand info."
    )
    parser.add_argument(
        "--outdir", "-o", default="gene_fastas",
        help="Output directory for gene FASTAs (default: gene_fastas/)."
    )
    parser.add_argument(
        "--copy", "-c", default=None,
        help="Process only this copy (e.g., Mcal_copy1). Default: all."
    )

    args = parser.parse_args()

    strand_map = parse_mRNA_beds(args.bed_dir)

    if not strand_map:
        sys.stderr.write(
            "Error: no strand information found in BED files in '%s'.\n"
            % args.bed_dir
        )
        sys.exit(1)

    print("Loaded strand info for %d copies from mRNA BED files." % len(strand_map))

    if args.copy:
        pattern = os.path.join(args.fasta_dir, "%s_exons.fasta" % args.copy)
    else:
        pattern = os.path.join(args.fasta_dir, "*_exons.fasta")

    fasta_files = sorted(glob.glob(pattern))

    if not fasta_files:
        sys.stderr.write(
            "Error: no exon FASTA files found matching '%s'.\n" % pattern
        )
        sys.exit(1)

    print("Found %d exon FASTA files to process.\n" % len(fasta_files))

    all_genes = []

    for fasta_path in fasta_files:
        basename = os.path.basename(fasta_path)
        copy_name = basename.replace("_exons.fasta", "")

        if copy_name not in strand_map:
            sys.stderr.write(
                "Warning: no strand info for '%s' in mRNA BED files. "
                "Skipping.\n" % copy_name
            )
            continue

        strand = strand_map[copy_name]
        records = parse_exon_fasta(fasta_path)

        if not records:
            sys.stderr.write(
                "Warning: no sequences found in '%s'. Skipping.\n"
                % fasta_path
            )
            continue

        outpath, gene_seq, sorted_records = concat_and_write(
            records, copy_name, strand, args.outdir
        )

        all_genes.append((copy_name, strand, gene_seq, sorted_records))

        exon_order = []
        for r in sorted_records:
            if r[2]:
                exon_order.append("%d-%d" % (r[2][1], r[2][2]))
            else:
                exon_order.append("unknown")

        print("  %s (strand: %s)" % (copy_name, strand))
        print("    Exons: %d" % len(sorted_records))
        print("    Gene length: %d bp" % len(gene_seq))
        print("    Exon order: %s" % " -> ".join(exon_order))
        print("    Written to: %s" % outpath)
        print("")

    if all_genes:
        combined_path = os.path.join(args.outdir, "all_ASCT_genes.fasta")
        with open(combined_path, "w") as fh:
            for copy_name, strand, gene_seq, _ in all_genes:
                fh.write(">%s\n" % copy_name)
                for i in range(0, len(gene_seq), 80):
                    fh.write(gene_seq[i:i + 80] + "\n")

        print("=" * 60)
        print("Combined FASTA: %s (%d genes)" % (
            combined_path, len(all_genes)
        ))
        print("=" * 60)


if __name__ == "__main__":
    main()
