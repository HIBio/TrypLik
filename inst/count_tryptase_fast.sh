#!/bin/bash

# set -e -u -o pipefail ## can't use this because ripgrep returns 1 if no matches

rg=$(which rg)
echo "Using $rg"

infile="$1"
echo "Processing $1"

## output dir
outdir="${2-.}"
echo "Using $outdir as output dir"

basename="$(basename $1)"
temp1="$outdir/$basename.temp1"
temp1a="$outdir/$basename.temp1a"
temp1b="$outdir/$basename.temp1b"
temp1c="$outdir/$basename.temp1c"
temp1d="$outdir/$basename.temp1d"
temp1e="$outdir/$basename.temp1e"
temp1f="$outdir/$basename.temp1f"
temp1g="$outdir/$basename.temp1g"
temp1h="$outdir/$basename.temp1h"
temp1i="$outdir/$basename.temp1i"
temp2="$outdir/$basename.temp2"

outfile="$outdir/${basename/.sam/.counts.txt}"

## rg added --include-zero in v12 but DNANexus uses v11.0.2 so as a
## workaround, use rg -c pattern file since no matches
## returns an exit code of 1

{ $rg -c 'CTGCAGC[A]AG[C]GGG[T]ATCGT[C]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGCCCTGGCAGGTGAGCCTGAGAGTCC[G]CG[A]CC[G]' $infile || echo 0; } > $temp1a
{ $rg -c 'CTGCAGC[G]AG[T]GGG[C]ATCGT[C]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGCCCTGGCAGGTGAGCCTGAGAGTCC[A]CG[G]CC[C]' $infile || echo 0; } > $temp1b
{ $rg -c 'CTGCAGC[G]AG[T]GGG[C]ATCGT[T]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGCCCTGGCAGGTGAGCCTGAGAGTCC[A]CG[G]CC[C]' $infile || echo 0; } > $temp1c
{ $rg -c 'CTGCAGC[G]AG[T]GGG[C]ATCGT[T]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGCCCTGGCAGGTGAGCCTGAGAGTCC[G]CG[A]CC[G]' $infile || echo 0; } > $temp1d

echo 0 > $temp1
$rg -c -F 'CCAGTCCAGGCCCTGCAGCAAGCGGGTATC' $infile >> $temp1
$rg -c 'C[AG]GGGCCTGGAGGGGTGGGCAAGGGCTGGA' $infile >> $temp1
$rg -c -F 'GACGTCAAGGATCTGGCCACCCTCAGGGTG' $infile >> $temp1
$rg -c -F 'AGTTCTACATCATCCAGACTGGAGCGGATA' $infile >> $temp1
$rg -c -F 'CCGCGTCCACACGGTCATGCTGCCCCCTGC' $infile >> $temp1
$rg -c -F 'TGGGGACAGTGGGAGGTGGGGCCAGGGTCT' $infile >> $temp1
$rg -c -F 'TTGCCCGGCCCCCTCCTCAGGCTGCACCCT' $infile >> $temp1
$rg -c 'TGCAGAGCCCCTC[CT]CACCGCCATTTCCCCT' $infile >> $temp1
$rg -c -F 'GGAGACGACGTCCGCATCATCCGTGACGAC' $infile >> $temp1
$rg -c -F 'CCAGGGCGACTCTGGAGGGCCCCTGGTGTG' $infile >> $temp1
$rg -c -F 'TACAGGCGGGCGTGGTCAGCTGGGACGAGG' $infile >> $temp1
paste -sd+ "$temp1" | bc > $temp1g
# awk '{ sum += $1 } END { print sum }' temp1 > temp1g

echo 0 > $temp2
$rg -c -F 'CCAGGCCAGGCCCTGCAGCGAGTGGGCATC' $infile >> $temp2
$rg -c -F 'CGGGGCCTGGAGGGGTGGGGAAGGGCTGGA' $infile >> $temp2
$rg -c -F 'GACGTCAAGGATCTGGCCGCCCTCAGGGTG' $infile >> $temp2
$rg -c -F 'AGTTCTACACCGCCCAGATCGGAGCGGACA' $infile >> $temp2
$rg -c -F 'CCACGTCCACACGGTCACCCTGCCCCCTGC' $infile >> $temp2
$rg -c -F 'TGGGGACAGTGGAGGTGGGGCCAGGGTCTT' $infile >> $temp2
$rg -c -F 'TTGCCCGGCCCCCTCCTGAGGCTGCACCCT' $infile >> $temp2
$rg -c -F 'TGCAGAGCGCCTCCCACCGCCATTTCCTCT' $infile >> $temp2
$rg -c -F 'GGAGACGACGTCCGCATCGTCCGTGACGAC' $infile >> $temp2
$rg -c -F 'CCAGGGCGACTCCGGAGGGCCCCTGGTGTG' $infile >> $temp2
$rg -c -F 'TGCAGGCGGGCGTGGTCAGCTGGGGCGAGG' $infile >> $temp2
paste -sd+ "$temp2" | bc > $temp1h
# awk '{ sum += $1 } END { print sum }' temp2 > temp1h

{ $rg -c -F 'CTCAGAGACCTTCCCCCCGGGGATGCCGT' $infile || echo 0; } > $temp1e
{ $rg -c -F 'CTCAGAGACCTTCCCCCCCGGGGATGCCG' $infile || echo 0; } > $temp1f

echo 0 > $temp1
$rg -c -F 'CCAGGCCAGGCCCTGCAGCAAACGGGCATT' $infile >> $temp1
$rg -c -F 'TGGGGCCTGGAGGGGTGGGCAAGGGCTGGA' $infile >> $temp1
$rg -c -F 'GACATCAAGGATCTGGCCGCCCTCAGGGTG' $infile >> $temp1
$rg -c -F 'AGTTCTACATCATCCAGACCGGGGCGGACA' $infile >> $temp1
$rg -c 'CCACATCCACAC[GC]GTCACGCTGCCCCCTGC' $infile >> $temp1
$rg -c -F 'GGGGACAGCGGGAGGCCGGGCCAGGTGGGC' $infile >> $temp1
$rg -c -F 'CCCAGCCGGCCCCAGACCCGGCTCCACGCC' $infile >> $temp1
$rg -c -F 'CCCAGTGCACCTGCCGCCGCCATACCCGCT' $infile >> $temp1
$rg -c -F 'GGCCACAGCTTTCAAATCGTCCGCGATGAC' $infile >> $temp1
$rg -c -F 'CCAGGGTGACTCTGGAGGGCCCCTGGTCTG' $infile >> $temp1
$rg -c -F 'TGCAGGCGGGCGTGGTCAGCTGGGAGGAGA' $infile >> $temp1
paste -sd+ "$temp1" | bc > $temp1i
# awk '{ sum += $1 } END { print sum }' temp1 > temp1i

# use `paste -s -` on mac
echo "Generating output for $basename"
cat $outdir/$basename.temp1[abcdefghi] | paste -s - > $outfile
cat $outfile
rm $outdir/$basename.temp[12]*
rm $infile

exit 0
