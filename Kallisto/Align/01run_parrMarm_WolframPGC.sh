#!/bin/bash
#!
#! Example SLURM job script for Gurdon Institute Cluster
#! Last updated: Sat Apr 18 13:05:53 BST 2015
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J PGCmeth
#! How many (MPI) tasks will there be in total? (<= nodes*16)
#SBATCH -n 6
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL
#SBATCH --mem=100G
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#! Job farming
#SBATCH --array=1-319%5



cd /mnt/scratch/gurdon/cap76/Thorsten/hPGCs/

FILES=( $(ls *_1.fastq.gz | sed -e 's/\_1.fastq.gz$//' ) )

mkdir fastqc_results2 
mkdir fastqc_trimmed_results2 
mkdir STAR_results2 

cp /mnt/scratch/gurdon/cap76/Thorsten/MissingTranscriptomics/wigToBigWig STAR_results2
cp /mnt/scratch/gurdon/cap76/Thorsten/MissingTranscriptomics/wigToBigWig fastqc_trimmed_results2

T=8 
T=8
wdir='/mnt/scratch/gurdon/cap76/Thorsten/hPGCs/'
FASTA='/mnt/scratch/gurdon/cap76/Wolfram/RNASeq/GRCh38.primary_assembly.genome.fa'
GTF='/mnt/scratch/gurdon/cap76/Wolfram/RNASeq/gencode.v31.annotation.gtf'
INDICES='/mnt/scratch/gurdon/cap76/Wolfram/RNASeq/GRCh38'
CHRSIZES='/mnt/scratch/gurdon/cap76/Wolfram/RNASeq/chr.size'

mystring="fastqc -o ${wdir}fastqc_results2 -t $T ${FILES[$SLURM_ARRAY_TASK_ID-1]}_1.fastq.gz ${FILES[$SLURM_ARRAY_TASK_ID-1]}_2.fastq.gz ; echo 'FASTQC has finished'  >> ${wdir}o.smartseq2 ; ~/trim_galore --paired  -o ${wdir}fastqc_trimmed_results2 ${FILES[$SLURM_ARRAY_TASK_ID-1]}_1.fastq.gz ${FILES[$SLURM_ARRAY_TASK_ID-1]}_2.fastq.gz ; echo 'trim_galore has finished' >> ${wdir}o.smartseq2 ; cd ${wdir}fastqc_trimmed_results2 ; cat ${FILES[$SLURM_ARRAY_TASK_ID-1]}_*.fq.gz_trimming_report.txt | grep -A2 'processed' > QC_fastq.${FILES[$SLURM_ARRAY_TASK_ID-1]}.txt ; cd ${wdir}fastqc_trimmed_results2"

echo $mystring

eval $mystring
