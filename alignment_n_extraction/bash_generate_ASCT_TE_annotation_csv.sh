# Process all RepeatMasker .out files into a single CSV
echo "species,paralog,contig,start,end,strand,repeat_name,repeat_class,sw_score,perc_div,perc_del,perc_ins" > ASCT_TE_landscape.csv

for out in rm_output/*.out; do
    # Extract species and paralog from filename
    # e.g., Mcal_ASCT1.3_region.fasta.out -> Mcal, ASCT1.3
    fname=$(basename "$out" | sed 's/_region.*//') 
    species=$(echo "$fname" | cut -d'_' -f1)
    paralog=$(echo "$fname" | cut -d'_' -f2-)
    
    awk -v sp="$species" -v para="$paralog" '
    NR>3 && NF>=15 {
        # Standard RepeatMasker .out columns
        # $1=score $2=%div $3=%del $4=%ins $5=query $6=start $7=end $8=(left) 
        # $9=strand $10=repeat_name $11=repeat_class/family
        strand = ($9 == "C") ? "-" : "+"
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", \
            sp, para, $5, $6, $7, strand, $10, $11, $1, $2, $3, $4
    }' "$out" >> ASCT_TE_landscape.csv
done

# Check output
head -n 5 ASCT_TE_landscape.csv
wc -l ASCT_TE_landscape.csv
