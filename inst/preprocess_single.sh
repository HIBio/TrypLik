#!/bin/bash

set -e -u -o pipefail

samtools index $1
samtools view -b -h $1 chr16:1200000-1300000 > $1.tryptase.cram
samtools index $1.tryptase.cram

bwa index consensus.fa

filtered=$1.tryptase.cram

samtools index $filtered
samtools sort -n $filtered -o ${filtered/.tryptase.cram/.tryptase.sorted.cram}
samtools fastq -1 tempR1.fastq -2 tempR2.fastq -0 /dev/null -s /dev/null -n ${filtered/.tryptase.cram/.tryptase.sorted.cram}
bwa mem consensus.fa tempR1.fastq tempR2.fastq > out.sam
samtools sort -n out.sam -o out.sorted.sam
samtools fixmate -m out.sorted.sam fixmate.sam
samtools sort -o positionsort.sam fixmate.sam
samtools markdup -r positionsort.sam temp.sam

exit 0
