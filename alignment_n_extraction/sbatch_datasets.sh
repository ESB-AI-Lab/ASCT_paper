#!/bin/bash
#SBATCH -J jjbdatadl		# jobname
#SBATCH -o jjbdatadl.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e jjddatadl.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 523-524			# start and stop of the array start-end
###SBATCH -n 1			# -n, --ntasks=INT Maximum number of tasks. Use for requesting a whole node. env var SLURM_NTASKS
#SBATCH -c 8			# -c, --cpus-per-task=INT The # of cpus/task. env var for threads is SLURM_CPUS_PER_TASK
#SBATCH -p p1			# queue (partition) -- normal, development, largemem, etc.
#SBATCH --mem-per-cpu=1750
###SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
###SBATCH --mail-user=solarese@uci.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes

#import jq
PATH=/gpool/bin/jq:$PATH

#import ncbi datasets
###datasets 18.14.0
PATH=/gpool/bin/ncbi:$PATH

# example call
# datasets download genome accession GCA_019925275.1 --include gff3,rna,cds,protein,genome,seq-report

#JOBFILE=dbs_dl.txt
JOBFILE=dbs_new.txt
#JOBFILE=dbs.missing.txt

#GCFID=$(datasets summary genome accession $GCAID --as-json-lines | jq -r '.paired_accession')
#GCAID=$1
GCAID=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)

echo $GCAID

datasets download genome accession $GCAID --include gtf,gff3,rna,cds,protein,genome,seq-report --filename ${GCAID}_ncbi_dataset.zip

GCFID=$(datasets summary genome accession $GCAID --as-json-lines | jq -r '.paired_accession')

echo $GCFID

if [[ -n "$GCFID" && "$GCFID" != "null" ]]; then
    datasets download genome accession $GCFID --include gtf,gff3,rna,cds,protein,genome,seq-report --filename ${GCFID}_ncbi_dataset.zip
else
    echo "No RefSeq accession found for $GCAID"
fi
