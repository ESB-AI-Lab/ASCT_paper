#!/bin/bash

#Mgal
GENOME=GCF_965363235.1_xbMytGall1.hap1.1_genomic.fasta
# ASCT2b1 on NC_134840.1 (57544849-57551566)
samtools faidx $GENOME NC_134840.1:57490000-57600000 > Mgal_ASCT2b1_region.fasta

# ASCT2a2 on NC_134850.1 (60701917-60709118)
samtools faidx $GENOME NC_134850.1:60650000-60760000 > Mgal_ASCT2a2_region.fasta

# ASCT1.3 on NC_134850.1 (21977828-21983936)
samtools faidx $GENOME NC_134850.1:21920000-22040000 > Mgal_ASCT1.3_region.fasta

#Mcal
GENOME=GCF_021869535.1_xbMytCali1.0.p_genomic.fasta
# ASCT2a2 on NW_026262593.1 (12323973-12333562)
samtools faidx $GENOME NW_026262593.1:12270000-12390000 > Mcal_ASCT2a2_region.fasta

# ASCT2b1 on NW_026262690.1 (65540866-65548530)
samtools faidx $GENOME NW_026262690.1:65490000-65600000 > Mcal_ASCT2b1_region.fasta

# ASCT1.3 on NW_026262593.1 (47660639-47670281)
samtools faidx $GENOME NW_026262593.1:47610000-47720000 > Mcal_ASCT1.3_region.fasta

# ============================================================
# Mcor
# ============================================================
GENOME=GCA_017311375.1_Mcoruscus_HiC_genomic.fna

# ASCT2b1 on CM029599.1 (57049002-57055636)
samtools faidx $GENOME CM029599.1:56999000-57106000 > Mcor_ASCT2b1_region.fasta

# ASCT1.3a on CM029607.1 (37886514-37900527)
samtools faidx $GENOME CM029607.1:37836000-37951000 > Mcor_ASCT1.3a_region.fasta

# ASCT1.3b on JAASAO010001760.1 (2940-11875) - small scaffold, pad what we can
samtools faidx $GENOME JAASAO010001760.1:1-62000 > Mcor_ASCT1.3b_region.fasta

# ============================================================
# Medi
# ============================================================
GENOME=GCA_019925275.1_PEIMed_genomic.fna

# ASCT2a1 on JAHUZM010000350.1 (29918-35992) - small scaffold, pad what we can
samtools faidx $GENOME JAHUZM010000350.1:1-86000 > Medi_ASCT2a1_region.fasta

# ASCT2b2 on CM034350.1 (60958994-60964221)
samtools faidx $GENOME CM034350.1:60909000-61015000 > Medi_ASCT2b2_region.fasta

# ASCT2b3 on CM034350.1 (57659492-57665833)
samtools faidx $GENOME CM034350.1:57609000-57716000 > Medi_ASCT2b3_region.fasta

# ASCT1.4 on CM034359.1 (29064098-29070814)
samtools faidx $GENOME CM034359.1:29014000-29121000 > Medi_ASCT1.4_region.fasta

# ASCT2a5 on CM034359.1 (73136389-73140720)
samtools faidx $GENOME CM034359.1:73086000-73191000 > Medi_ASCT2a5_region.fasta

# ============================================================
# Pvir - ASCT2b1 and ASCT2a2 are ~13kb apart, extract as separate regions
# ============================================================
GENOME=GCA_037379345.1_ASM3737934v1_genomic.fasta

# ASCT2b1 on CM074522.1 (8171380-8175471)
samtools faidx $GENOME CM074522.1:8121000-8226000 > Pvir_ASCT2b1_region.fasta

# ASCT2a2 on CM074522.1 (8188527-8193638)
samtools faidx $GENOME CM074522.1:8138000-8244000 > Pvir_ASCT2a2_region.fasta

# ASCT1.3 on CM074522.1 (14292138-14302028)
samtools faidx $GENOME CM074522.1:14242000-14352000 > Pvir_ASCT1.3_region.fasta

# ============================================================
# Asen - ASCT2a1 and ASCT2b2 are ~12kb apart, extract as separate regions
# ============================================================
GENOME=GCA_963971305.1_xbArcSenh1.1_genomic.fasta

# ASCT2a1 on OZ020316.1 (42032990-42043786)
samtools faidx $GENOME OZ020316.1:41983000-42094000 > Asen_ASCT2a1_region.fasta

# ASCT2b2 on OZ020316.1 (42055951-42062003)
samtools faidx $GENOME OZ020316.1:42006000-42112000 > Asen_ASCT2b2_region.fasta

# ASCT1.3 on OZ020316.1 (38094223-38103125)
samtools faidx $GENOME OZ020316.1:38044000-38153000 > Asen_ASCT1.3_region.fasta

# ============================================================
# Run RepeatMasker on all new regions
# ============================================================
mkdir -p rm_output

for fasta in Mcor_*_region.fasta Medi_*_region.fasta Pvir_*_region.fasta Asen_*_region.fasta; do
    singularity exec \
        --bind $(pwd):/home/edwin \
        --bind /gpool/shared/repbase:/gpool/shared/repbase \
        /gpool/bin/singularityimages/repeatmodeler.sif \
        RepeatMasker -lib /gpool/shared/repbase/repbase_master.fasta \
        -gff -xsmall -dir ./rm_output "$fasta"
done

# ============================================================
# After RepeatMasker finishes, append new results to existing CSV
# ============================================================
for out in rm_output/Mcor_*.out rm_output/Medi_*.out rm_output/Pvir_*.out rm_output/Asen_*.out; do
    fname=$(basename "$out" | sed 's/_region.*//') 
    species=$(echo "$fname" | cut -d'_' -f1)
    paralog=$(echo "$fname" | cut -d'_' -f2-)
    
    awk -v sp="$species" -v para="$paralog" '
    NR>3 && NF>=15 {
        strand = ($9 == "C") ? "-" : "+"
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", \
            sp, para, $5, $6, $7, strand, $10, $11, $1, $2, $3, $4
    }' "$out" >> ASCT_TE_landscape.csv
done

# Quick summary of new results
for out in rm_output/Mcor_*.out rm_output/Medi_*.out rm_output/Pvir_*.out rm_output/Asen_*.out; do
    echo "=== $(basename $out .fasta.out) ==="
    awk 'NR>3 {print $10, $11}' "$out" | sort | uniq -c | sort -rn | head -5
    echo ""
done
