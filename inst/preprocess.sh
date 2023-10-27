#!/bin/bash

set -e -u -o pipefail

for x in *.cram; do
  samtools index $x
  samtools view -b -h $x chr16:1200000-1300000 > $x.tryptase.cram
  samtools index $x.tryptase.cram
done

bwa index consensus.fa

for x in *tryptase.cram; do
  samtools index $x
  samtools sort -n $x -o ${x/.cram/.sorted.cram}
  samtools fastq -1 tempR1.fastq -2 tempR2.fastq -0 /dev/null -s /dev/null -n ${x/.cram/.sorted.cram}
  bwa mem consensus.fa tempR1.fastq tempR2.fastq > out.sam
  samtools sort -n out.sam -o out.sorted.sam
  samtools fixmate -m out.sorted.sam fixmate.sam
  samtools sort -o positionsort.sam fixmate.sam
  samtools markdup -r positionsort.sam temp.sam
  ./count_tryptase.sh >> tryptase.out
done

cat tryptase.out | awk '{printf "./Tryplik %d %d %d %d %d | sort -nk2 | tail -n 1\n", $6, $7, $8, $9, $10}' > in1

chmod +x in1

./in1 > tryptaseMLE

exit 0
