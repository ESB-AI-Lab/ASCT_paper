#!/bin/bash

for outfn in $(ls -t *.outfmt6 | head -n 17); do 
   species=$(grep $(echo $outfn | awk -F "_" '{print $7}') generate_asct_bed.py | awk '{print $3}' | sed 's/"//g' | sed 's/,//g')
   python3.9 generate_asct_bed.py --input $outfn --species $species --outdir bed_files/
done > extract_gene_copies.txt
