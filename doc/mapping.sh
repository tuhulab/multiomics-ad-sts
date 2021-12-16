#!/bin/sh
module load gcc/9.4.0 star/2.7.9a parallel/20210722
cd /home/projects/ku_00015/people/tuhu/multiomics-ad-sts/data
fn_path=../_tmp/fn.txt
# cat ${fn_path} | parallel --tempdir ../_tmp "STAR --runMode alignReads --genomeDir STAR_Index/ --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix BAM/{} --readFilesIn fastq/{}_R2_001.fastq.gz --runThreadN 20 --outTmpDir /home/projects/ku_00015/people/tuhu/multiomics-ad-sts/_tmp/STAR/{}"
STAR --genomeLoad LoadAndExit --genomeDir STAR_Index/
for i in $(cat ${fn_path}); do STAR --genomeDir STAR_Index/ --genomeLoad LoadAndKeep\
--runMode alignReads \
--outFilterMultimapNmax 1 \
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--outFileNamePrefix BAM/$i. \
--readFilesIn fastq/$i\_R2_001.fastq.gz \
--runThreadN 20 ; done
STAR --genomeLoad Remove --genomeDir STAR_Index/