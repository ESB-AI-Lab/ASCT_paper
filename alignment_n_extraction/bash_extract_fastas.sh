#!/bin/bash

GENOME=GCF_021869535.1_xbMytCali1.0.p_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCF_965363235.1_xbMytGall1.hap1.1_genomic.fasta
for bed in $(ls bed_files/*Mgal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_963971305.1_xbArcSenh1.1_genomic.fasta
for bed in $(ls bed_files/*Asen*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_037379345.1_ASM3737934v1_genomic.fasta
for bed in $(ls bed_files/*Pvir*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_019925275.1_PEIMed_genomic.fna
for bed in $(ls bed_files/*Medi*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_017311375.1_Mcoruscus_HiC_genomic.fna
for bed in $(ls bed_files/*Mcor*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_948099385.2_Fh240107Braker_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_048127485.1_Solemya_velum_v1_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCF_932274485.2_xgPatVulg1.2_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCF_023055435.1_xgHalRufe1.0.p_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCF_003073045.1_ASM307304v1_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCF_006345805.1_ASM634580v1_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_047652355.1_ASM4765235v1_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCF_032854445.1_CUHK_Ljap_v2_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_903813345.2_Owenia_chromosome_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_978046135.1_Genome_assembly_of_Meretrix_meretrix_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_029931535.1_MarmarV2_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCA_030141615.1_Upi_BIV9798_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCF_021730395.1_MADL_Memer_1_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCF_026571515.1_ASM2657151v2_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCF_902652985.1_xPecMax1.1_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

GENOME=GCF_963853765.1_xbMagGiga1.1_genomic.fasta
for bed in $(ls bed_files/*Mcal*exons.bed); do bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; done

for outfn in $(ls -t *.outfmt6 | head -n 14); do 
    species=$(grep $(echo $outfn | awk -F "_" '{print $7}') generate_asct_bed.py | awk '{print $3}' | sed 's/"//g' | sed 's/,//g');
    GENOME=$(ls $(echo $outfn | awk -F "_" '{print $5"_"$6}')*.fasta);
    for bed in $(ls bed_files/*${species}*exons.bed); do 
        bedtools getfasta -s -fi $GENOME -bed $bed -fo bed_files/$(basename $bed .bed).fasta; 
    done; 
done

