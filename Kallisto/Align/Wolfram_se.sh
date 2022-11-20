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
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#SBATCH --mem=100G
#! Job farming
#SBATCH --array=1-2%10


cd /mnt/scratch/gurdon/cap76/Wolfram/RNASeq

m=( $(ls SRR153729*.fastq.gz | sed -e 's/\.fastq.gz$//' ) )
FILES=( $(ls SRR153729*.fastq.gz | sed -e 's/\.fastq.gz$//' ) )
a=($(echo $t | tr ',' "\n"))
mkdir fastqc_results3
mkdir fastqc_trimmed_results3
mkdir STAR_results3

cp /mnt/scratch/gurdon/cap76/Thorsten/MissingTranscriptomics/wigToBigWig STAR_results2
cp /mnt/scratch/gurdon/cap76/Thorsten/MissingTranscriptomics/wigToBigWig fastqc_trimmed_results2

T=8
wdir='/mnt/scratch/gurdon/cap76/Wolfram/RNASeq/'
FASTA='/mnt/scratch/gurdon/cap76/Wolfram/RNASeq/GRCh38.primary_assembly.genome.fa'
GTF='/mnt/scratch/gurdon/cap76/Wolfram/RNASeq/gencode.v31.annotation.gtf'
INDICES='/mnt/scratch/gurdon/cap76/Wolfram/RNASeq/GRCh38'
CHRSIZES='/mnt/scratch/gurdon/cap76/Wolfram/RNASeq/chr.size'
mystring="fastqc -o ${wdir}fastqc_results3 -t $T ${FILES[$SLURM_ARRAY_TASK_ID-1]}.fastq.gz ; echo 'FASTQC has finished'  >> ${wdir}o.smartseq2 ; ~/trim_galore --fastqc -o ${wdir}fastqc_trimmed_results3 ${FILES[$SLURM_ARRAY_TASK_ID-1]}.fastq.gz ; echo 'trim_galore has finished' >> ${wdir}o.smartseq2 ; cd ${wdir}fastqc_trimmed_results3 ; cat ${FILES[$SLURM_ARRAY_TASK_ID-1]}_*.fastq.gz_trimming_report.txt | grep -A2 'processed' > QC_fastq.${FILES[$SLURM_ARRAY_TASK_ID-1]}.txt; STAR --runThreadN $T --runMode alignReads --genomeDir $INDICES --readFilesIn ${FILES[$SLURM_ARRAY_TASK_ID-1]}_trimmed.fq.gz  --readFilesCommand gunzip -c --outFileNamePrefix ${wdir}STAR_results3/${FILES[$SLURM_ARRAY_TASK_ID-1]} --sjdbGTFfile $GTF --sjdbOverhang 149 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outWigType wiggle --outWigStrand Unstranded --outWigNorm RPM ; echo 'STAR has finished' >> ${wdir}o.smartseq2 ; cd .. ; cd STAR_results3 ; cat  ${FILES[$SLURM_ARRAY_TASK_ID-1]}Log.final.out | grep 'Uniquely mapped reads' > QC_uniquelyMappedReads.${FILES[$SLURM_ARRAY_TASK_ID-1]}.txt ; samtools index ${FILES[$SLURM_ARRAY_TASK_ID-1]}Aligned.sortedByCoord.out.bam ; ./wigToBigWig ${FILES[$SLURM_ARRAY_TASK_ID-1]}Signal.UniqueMultiple.str1.out.wig $CHRSIZES ${sample}_uniqueMultiple.bw ; ./wigToBigWig ${FILES[$SLURM_ARRAY_TASK_ID-1]}Signal.Unique.str1.out.wig $CHRSIZES ${FILES[$SLURM_ARRAY_TASK_ID-1]}_unique.bw ; cd .. ; echo 'Smart-SEQ2 analysis completed...' >> ${wdir}o.smartseq2 ; echo ${FILES[$SLURM_ARRAY_TASK_ID-1]} >> ${wdir}o.smartseq2 ; echo '--------------------------------' >> ${wdir}o.smartseq2"
echo $mystring

eval $mystring
