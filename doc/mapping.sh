#!/bin/sh
module load gcc/9.4.0 star/2.7.9a parallel/20210722
cd /home/projects/ku_00015/people/tuhu/multiomics-ad-sts/data
fn_path=../_tmp/fn.txt
for i in $(cat ${fn_path}); do STAR --genomeDir STAR_Index/ \
--runMode alignReads \
--outFilterMultimapNmax 1 \
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--outFileNamePrefix BAM/$i. \
--readFilesIn fastq/$i\_R2_001.fastq.gz \
--runThreadN 20 ; done
