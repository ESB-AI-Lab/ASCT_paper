#!/bin/bash
#SBATCH -J kwasiblastdb		# jobname
#SBATCH -o kwasiblastdb.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e kwasiblastdb.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1			# start and stop of the array start-end
###SBATCH -n 1			# -n, --ntasks=INT Maximum number of tasks. Use for requesting a whole node. env var SLURM_NTASKS
#SBATCH -c 8			# -c, --cpus-per-task=INT The # of cpus/task. env var for threads is SLURM_CPUS_PER_TASK
#SBATCH -p p1			# queue (partition) -- normal, development, largemem, etc.
#SBATCH --mem-per-cpu=1750
###SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
###SBATCH --mail-user=solarese@uci.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes

#import blastn
###blast 2.2.31+
PATH=/gpool/bin/ncbi-blast-2.2.31+/bin:$PATH

bestn=1

JOBFILE=exon_queries.txt
#MYDB=GHIK01.1_sl.fasta
#MYDB=GAEM01.1.fasta
#MYDB=GHII01.1.fasta
#MYLABEL="GAEM01TSA"
#MYLABEL="GHII01"
QUERY=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)
#QUERY=gracey2ndlibsequence.fasta
#QUERY=Mytilus_californianus_andy_gracey_11k_sl.fasta
#QUERY=Mytilus_californianus_ESTs_45k_sl.fasta

# Retrotransposon search
#MYDB=Mytilus_californianus_ESTs_45k_sl.fasta
#MYLABEL="m_californianus_ests"
#QUERY=retroviridae.fasta
#MYDB=retroviridae.fasta
#MYLABEL="retroviridae"
#QUERY=Mytilus_californianus_ESTs_45k_sl.fasta
#QUERY=GHII01.1.fasta
#MYDB=GHII01.1.fasta


# Viral Search
#MYDB=viral.1.1.genomic.fasta
#MYLABEL="viral_db1"
#MYDB=viral.2.1.genomic.fasta
#MYLABEL="viral_db2"
#QUERY=viral.1.1.genomic.fasta
#QUERY=viral.2.1.genomic.fasta

# worm genes
#QUERY=Tb11.02.0290.gene.fasta
#QUERY=Tb11.02.0290.pep.fasta
#QUERY=Tb11.02.0290.rna.fasta
# Mytilus Gali Transcriptome
#MYDB=GHII01.1.fasta
#MYDB=GAEM01.1.fasta
#MYDB=GHIK01.1_sl.fasta

# Genome Databases

#QUERY=Mytilus_californianus_ESTs_45k_sl.fasta
#MYDB=Mytilus_californianus_ESTs_45k_sl.fasta
#MYLABEL="m_californianus_ests"
#QUERY=GHII01.1.fasta
# Mytilus Gali Transcriptome
#MYDB=GHII01.1.fasta

#MYLABEL="m_californianus_v1p"
#MYDB=GCA_021869535.1_xbMytCali1.0.p_genomic.fna
#MYLABEL="m_californianus_v1p"
#MYDB=GCF_021869535.1_xbMytCali1.0.p_genomic.fasta
#MYLABEL="m_coruscus_v1"
#MYDB=GCA_017311375.1_Mcoruscus_HiC_genomic.fna
#MYLABEL="m_edulis_v1"
#MYDB=GCA_019925275.1_PEIMed_genomic.fna
#MYLABEL="m_edulis_v2"
#MYDB=GCA_019925275.2_PEIMed_v2_genomic.fna
#MYLABEL="m_californianus_v1a"
#MYDB=GCA_021869935.1_xbMytCali1.0.a_genomic.fna
#MYLABEL="m_galloprovincialis_v1"
#MYDB=GCF_965363235.1_xbMytGall1.hap1.1_genomic.fasta
#MYLABEL="perna_viridis_v1"
#MYDB=GCA_037379345.1_ASM3737934v1_genomic.fasta
#MYLABEL="arcuatula_senhousia_v1"
#MYDB=GCA_963971305.1_xbArcSenh1.1_genomic.fasta
MYLABEL="Fasciola_hepatica"
MYDB="GCA_948099385.2_Fh240107Braker_genomic.fasta"
MYLABEL="Solemya_velum"
MYDB="GCA_048127485.1_Solemya_velum_v1_genomic.fasta"
MYLABEL="Patella_vulgata"
MYDB="GCF_932274485.2_xgPatVulg1.2_genomic.fasta"
MYLABEL="Haliotis_rufescens"
MYDB="GCF_023055435.1_xgHalRufe1.0.p_genomic.fasta"
MYLABEL="Pomacea_canaliculata"
MYDB="GCF_003073045.1_ASM307304v1_genomic.fasta"
MYLABEL="Octopus_sinensis"
MYDB="GCF_006345805.1_ASM634580v1_genomic.fasta"
MYLABEL="Nautilus_pompilius"
MYDB="GCA_047652355.1_ASM4765235v1_genomic.fasta"
MYLABEL="Liolophura_sinensis"
MYDB="GCF_032854445.1_CUHK_Ljap_v2_genomic.fasta"
MYLABEL="Owenia_fusiformis"
MYDB="GCA_903813345.2_Owenia_chromosome_genomic.fasta"
MYLABEL="Meretrix_meretrix"
MYDB="GCA_978046135.1_Genome_assembly_of_Meretrix_meretrix_genomic.fasta"
#
#MYLABEL="Margaritifera_margaritifera"
#MYDB="GCA_029931535.1_MarmarV2_genomic.fasta"
#MYLABEL="Unio_pictorum"
#MYDB="GCA_030141615.1_Upi_BIV9798_genomic.fasta"
#MYLABEL="Mercenaria_mercenaria"
#MYDB="GCF_021730395.1_MADL_Memer_1_genomic.fasta"
#MYLABEL="Ruditapes_philippinarum"
#MYDB="GCF_026571515.1_ASM2657151v2_genomic.fasta"
#MYLABEL="Pecten_maximus"
#MYDB="GCF_902652985.1_xPecMax1.1_genomic.fasta"
#MYLABEL="Magallana_gigas"
#MYDB="GCF_963853765.1_xbMagGiga1.1_genomic.fasta"
#MYLABEL="Argopecten_irradians"
#MYDB="GCF_041381155.1_Ai_NY_genomic.fasta"
#MYLABEL=""
#MYDB=""

makeblastdb -in ${MYDB} -input_type fasta -dbtype nucl -parse_seqids -out $(basename ${MYDB} .fasta) -title ${MYLABEL}
#blastn -query ${QUERY} -db $(basename ${MYDB} .fasta)  -max_target_seqs ${bestn} -outfmt 6 -evalue 1e-5 -num_threads ${SLURM_CPUS_PER_TASK} > $(basename ${QUERY} .fasta)_$(basename ${MYDB} .fasta)_e5_blastn_best${bestn}.outfmt6
#blastn -query ${QUERY} -db $(basename ${MYDB} .fasta)  -max_target_seqs ${bestn} -outfmt 6 -evalue 1 -num_threads ${SLURM_CPUS_PER_TASK} > $(basename ${QUERY} .fasta)_$(basename ${MYDB} .fasta)_e1_blastn_best${bestn}.outfmt6
#tblastx -query ${QUERY} -db $(basename ${MYDB} .fasta)  -max_target_seqs ${bestn} -outfmt 6 -evalue 1 -num_threads ${SLURM_CPUS_PER_TASK} > $(basename ${QUERY} .fasta)_$(basename ${MYDB} .fasta)_tblastxe1_blastn_best${bestn}.outfmt6
