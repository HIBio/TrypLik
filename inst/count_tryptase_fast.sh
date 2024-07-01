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

outfile="$outdir/${basename/.sam/.counts}"

$rg -c --include-zero 'CTGCAGC[A]AG[C]GGG[T]ATCGT[C]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGCCCTGGCAGGTGAGCCTGAGAGTCC[G]CG[A]CC[G]' $infile > $temp1a
$rg -c --include-zero 'CTGCAGC[G]AG[T]GGG[C]ATCGT[C]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGCCCTGGCAGGTGAGCCTGAGAGTCC[A]CG[G]CC[C]' $infile > $temp1b
$rg -c --include-zero 'CTGCAGC[G]AG[T]GGG[C]ATCGT[T]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGCCCTGGCAGGTGAGCCTGAGAGTCC[A]CG[G]CC[C]' $infile > $temp1c
$rg -c --include-zero 'CTGCAGC[G]AG[T]GGG[C]ATCGT[T]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGCCCTGGCAGGTGAGCCTGAGAGTCC[G]CG[A]CC[G]' $infile > $temp1d

$rg -c --include-zero -F 'CCAGTCCAGGCCCTGCAGCAAGCGGGTATC' $infile > $temp1
$rg -c --include-zero 'C[AG]GGGCCTGGAGGGGTGGGCAAGGGCTGGA' $infile >> $temp1
$rg -c --include-zero -F 'GACGTCAAGGATCTGGCCACCCTCAGGGTG' $infile >> $temp1
$rg -c --include-zero -F 'AGTTCTACATCATCCAGACTGGAGCGGATA' $infile >> $temp1
$rg -c --include-zero -F 'CCGCGTCCACACGGTCATGCTGCCCCCTGC' $infile >> $temp1
$rg -c --include-zero -F 'TGGGGACAGTGGGAGGTGGGGCCAGGGTCT' $infile >> $temp1
$rg -c --include-zero -F 'TTGCCCGGCCCCCTCCTCAGGCTGCACCCT' $infile >> $temp1
$rg -c --include-zero 'TGCAGAGCCCCTC[CT]CACCGCCATTTCCCCT' $infile >> $temp1
$rg -c --include-zero -F 'GGAGACGACGTCCGCATCATCCGTGACGAC' $infile >> $temp1
$rg -c --include-zero -F 'CCAGGGCGACTCTGGAGGGCCCCTGGTGTG' $infile >> $temp1
$rg -c --include-zero -F 'TACAGGCGGGCGTGGTCAGCTGGGACGAGG' $infile >> $temp1
paste -sd+ $temp1 | bc > $temp1g
# awk '{ sum += $1 } END { print sum }' temp1 > temp1g

$rg -c --include-zero -F 'CCAGGCCAGGCCCTGCAGCGAGTGGGCATC' $infile > $temp2
$rg -c --include-zero -F 'CGGGGCCTGGAGGGGTGGGGAAGGGCTGGA' $infile >> $temp2
$rg -c --include-zero -F 'GACGTCAAGGATCTGGCCGCCCTCAGGGTG' $infile >> $temp2
$rg -c --include-zero -F 'AGTTCTACACCGCCCAGATCGGAGCGGACA' $infile >> $temp2
$rg -c --include-zero -F 'CCACGTCCACACGGTCACCCTGCCCCCTGC' $infile >> $temp2
$rg -c --include-zero -F 'TGGGGACAGTGGAGGTGGGGCCAGGGTCTT' $infile >> $temp2
$rg -c --include-zero -F 'TTGCCCGGCCCCCTCCTGAGGCTGCACCCT' $infile >> $temp2
$rg -c --include-zero -F 'TGCAGAGCGCCTCCCACCGCCATTTCCTCT' $infile >> $temp2
$rg -c --include-zero -F 'GGAGACGACGTCCGCATCGTCCGTGACGAC' $infile >> $temp2
$rg -c --include-zero -F 'CCAGGGCGACTCCGGAGGGCCCCTGGTGTG' $infile >> $temp2
$rg -c --include-zero -F 'TGCAGGCGGGCGTGGTCAGCTGGGGCGAGG' $infile >> $temp2
paste -sd+ $temp2 | bc > $temp1h
# awk '{ sum += $1 } END { print sum }' temp2 > temp1h

$rg -c --include-zero -F 'CTCAGAGACCTTCCCCCCGGGGATGCCGT' $infile > $temp1e
$rg -c --include-zero -F 'CTCAGAGACCTTCCCCCCCGGGGATGCCG' $infile > $temp1f

$rg -c --include-zero -F 'CCAGGCCAGGCCCTGCAGCAAACGGGCATT' $infile  > $temp1
$rg -c --include-zero -F 'TGGGGCCTGGAGGGGTGGGCAAGGGCTGGA' $infile  >> $temp1
$rg -c --include-zero -F 'GACATCAAGGATCTGGCCGCCCTCAGGGTG' $infile  >> $temp1
$rg -c --include-zero -F 'AGTTCTACATCATCCAGACCGGGGCGGACA' $infile  >> $temp1
$rg -c --include-zero 'CCACATCCACAC[GC]GTCACGCTGCCCCCTGC' $infile >> $temp1
$rg -c --include-zero -F 'GGGGACAGCGGGAGGCCGGGCCAGGTGGGC' $infile  >> $temp1
$rg -c --include-zero -F 'CCCAGCCGGCCCCAGACCCGGCTCCACGCC' $infile  >> $temp1
$rg -c --include-zero -F 'CCCAGTGCACCTGCCGCCGCCATACCCGCT' $infile  >> $temp1
$rg -c --include-zero -F 'GGCCACAGCTTTCAAATCGTCCGCGATGAC' $infile  >> $temp1
$rg -c --include-zero -F 'CCAGGGTGACTCTGGAGGGCCCCTGGTCTG' $infile  >> $temp1
$rg -c --include-zero -F 'TGCAGGCGGGCGTGGTCAGCTGGGAGGAGA' $infile  >> $temp1
paste -sd+ $temp1 | bc > $temp1i
# awk '{ sum += $1 } END { print sum }' temp1 > temp1i

# use `paste -s -` on mac
cat $outdir/$basename.temp1[abcdefghi] | paste -s - > $outfile
cat $outfile
rm $outdir/$basename.temp[12]*

exit 0
