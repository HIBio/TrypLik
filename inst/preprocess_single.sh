#!/bin/bash

# set -e -u -o pipefail

nthreads=${4-8}

SAMTOOLS=$(which samtools)
BWA="/home/dnanexus/out/out/bwa/bwa"

echo "Using samtools: $SAMTOOLS"
echo "Using bwa: $BWA"

## output dir
outdir="${2-.}"
echo "Using $outdir as output dir"

consensusdir="${3-.}"
echo "consensusdir is $consensusdir"

# if [[ -f "$consensusdir/consensus.fa.amb" ]]; then
#   echo "consensus file already indexed"
# else
#   echo "Processing consensus..."
#   $BWA index $consensusdir/consensus.fa
# fi

## index assumed to be in the same dir as .cram
if [[ -f "$1.crai" ]]; then
  echo ".crai index file already exists."
else
  $SAMTOOLS index "$1"
fi

## temp files and output
basename="$(basename $1)"
filtered="$outdir/$basename.tryptase.cram"
sorted="${filtered/.tryptase.cram/.tryptase.sorted.cram}"
fastq1="$outdir/$basename.tempR1.fastq"
fastq2="$outdir/$basename.tempR2.fastq"
out="$outdir/$basename.out.sam"
outsorted="$outdir/$basename.out.sorted.sam"
fixmate="$outdir/$basename.fixmate.sam"
possort="$outdir/$basename.positionsort.sam"
output="$outdir/${basename/.cram/.sam}"

# not indexing $filtered as it takes time
$SAMTOOLS view -u --fast -b -h "$1" chr16:1200000-1300000 > $filtered
$SAMTOOLS sort -u -@$nthreads -m1G -n $filtered -o $sorted
$SAMTOOLS fastq -1 $fastq1 -2 $fastq2 -0 /dev/null -s /dev/null -n $sorted
$BWA mem -t$nthreads $consensusdir/consensus.fa $fastq1 $fastq2 > $out
$SAMTOOLS sort -u -@$nthreads -m1G -n $out -o $outsorted
$SAMTOOLS fixmate -m $outsorted $fixmate
$SAMTOOLS sort -u -@$nthreads -m1G -o $possort $fixmate
$SAMTOOLS markdup -r $possort $output
echo "Done with samtools for $1"

## clean up temporary files:
echo "Cleanup"
rm $filtered
rm $sorted
rm $fastq1
rm $fastq2
rm $out
rm $outsorted
rm $fixmate
rm $possort

echo "Completed."

exit 0
