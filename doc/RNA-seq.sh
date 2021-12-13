#!/bin/sh
module load jdk/17 perl/5.30.2 fastqc/0.11.9
cd /home/projects/ku_00015/people/tuhu/multiomics-ad-sts
mkdir -p data/fastqc
fastqc data/fastq/*.fastq.gz -o