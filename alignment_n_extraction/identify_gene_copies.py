#!/usr/bin/env python3
import sys
from collections import defaultdict

def group_gene_copies(bedfile, max_gap=20000):
    # Read all lines
    lines = []
    with open(bedfile) as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            strand = fields[6]
            lines.append((chrom, start, end, strand, line.strip()))
    
    # Group by chromosome first
    chrom_groups = defaultdict(list)
    for chrom, start, end, strand, original_line in lines:
        chrom_groups[chrom].append((start, end, strand, original_line))
    
    results = []
    global_copy_counter = 0
    
    # Process each chromosome
    for chrom in sorted(chrom_groups.keys()):
        # Group by strand within chromosome
        strand_groups = defaultdict(list)
        for start, end, strand, original_line in chrom_groups[chrom]:
            strand_groups[strand].append((start, end, original_line))
        
        # Get all loci (strand+position combos) and sort by genomic position
        all_loci = []
        for strand in strand_groups:
            exons = strand_groups[strand]
            exons.sort(key=lambda x: x[0])
            
            # Group exons into loci based on gaps
            current_locus = []
            for start, end, original_line in exons:
                if current_locus and start - current_locus[-1][1] > max_gap:
                    # Close current locus, start new one
                    all_loci.append((current_locus[0][0], strand, current_locus))
                    current_locus = []
                current_locus.append((start, end, original_line))
            
            if current_locus:
                all_loci.append((current_locus[0][0], strand, current_locus))
        
        # Sort loci by genomic position
        all_loci.sort(key=lambda x: x[0])
        
        # Assign copy numbers in genomic order
        for start_pos, strand, locus_exons in all_loci:
            global_copy_counter += 1
            
            for start, end, original_line in locus_exons:
                fields = original_line.split('\t')
                base_name = fields[4]
                fields[4] = f"{base_name}_copy{global_copy_counter}"
                results.append('\t'.join(fields))
    
    return results

if __name__ == "__main__":
    bedfile = sys.argv[1]
    max_gap = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    
    results = group_gene_copies(bedfile, max_gap)
    
    # Sort output for readability
    def sort_key(line):
        fields = line.split('\t')
        return (fields[0], int(fields[1]))
    
    for line in sorted(results, key=sort_key):
        print(line)
