#!/bin/bash

set -e -u -o pipefail -o noclobber

cat temp.sam | grep CTGCAGC[A]AG[C]GGG[T]ATCGT[C]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGC
CCTGGCAGGTGAGCCTGAGAGTCC[G]CG[A]CC[G] | wc | awk '{print $1}' > temp1a
cat temp.sam | grep CTGCAGC[G]AG[T]GGG[C]ATCGT[C]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGC
CCTGGCAGGTGAGCCTGAGAGTCC[A]CG[G]CC[C] | wc | awk '{print $1}' > temp1b
cat temp.sam | grep CTGCAGC[G]AG[T]GGG[C]ATCGT[T]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGC
CCTGGCAGGTGAGCCTGAGAGTCC[A]CG[G]CC[C] | wc | awk '{print $1}' > temp1c
cat temp.sam | grep CTGCAGC[G]AG[T]GGG[C]ATCGT[T]GGGGGTCAGGAGGCCCCCAGGAGCAAGTGGC
CCTGGCAGGTGAGCCTGAGAGTCC[G]CG[A]CC[G] | wc | awk '{print $1}' > temp1d

cat temp.sam | grep CCAGTCCAGGCCCTGCAGCAAGCGGGTATC > temp1
cat temp.sam | grep C[AG]GGGCCTGGAGGGGTGGGCAAGGGCTGGA >> temp1
cat temp.sam | grep GACGTCAAGGATCTGGCCACCCTCAGGGTG >> temp1
cat temp.sam | grep AGTTCTACATCATCCAGACTGGAGCGGATA >> temp1
cat temp.sam | grep CCGCGTCCACACGGTCATGCTGCCCCCTGC >> temp1
cat temp.sam | grep TGGGGACAGTGGGAGGTGGGGCCAGGGTCT >> temp1
cat temp.sam | grep TTGCCCGGCCCCCTCCTCAGGCTGCACCCT >> temp1
cat temp.sam | grep TGCAGAGCCCCTC[CT]CACCGCCATTTCCCCT >> temp1
cat temp.sam | grep GGAGACGACGTCCGCATCATCCGTGACGAC >> temp1
cat temp.sam | grep CCAGGGCGACTCTGGAGGGCCCCTGGTGTG >> temp1
cat temp.sam | grep TACAGGCGGGCGTGGTCAGCTGGGACGAGG >> temp1
cat temp1 | wc | awk '{print $1}' > temp1g

cat temp.sam | grep CCAGGCCAGGCCCTGCAGCGAGTGGGCATC > temp2
cat temp.sam | grep CGGGGCCTGGAGGGGTGGGGAAGGGCTGGA >> temp2
cat temp.sam | grep GACGTCAAGGATCTGGCCGCCCTCAGGGTG >> temp2
cat temp.sam | grep AGTTCTACACCGCCCAGATCGGAGCGGACA >> temp2
cat temp.sam | grep CCACGTCCACACGGTCACCCTGCCCCCTGC >> temp2
cat temp.sam | grep TGGGGACAGTGGAGGTGGGGCCAGGGTCTT >> temp2
cat temp.sam | grep TTGCCCGGCCCCCTCCTGAGGCTGCACCCT >> temp2
cat temp.sam | grep TGCAGAGCGCCTCCCACCGCCATTTCCTCT >> temp2
cat temp.sam | grep GGAGACGACGTCCGCATCGTCCGTGACGAC >> temp2
cat temp.sam | grep CCAGGGCGACTCCGGAGGGCCCCTGGTGTG >> temp2
cat temp.sam | grep TGCAGGCGGGCGTGGTCAGCTGGGGCGAGG >> temp2
cat temp2 | wc | awk '{print $1}' > temp1h

cat temp.sam | grep CTCAGAGACCTTCCCCCCGGGGATGCCGT | wc | awk '{print $1}' > temp1e
cat temp.sam | grep CTCAGAGACCTTCCCCCCCGGGGATGCCG | wc | awk '{print $1}' > temp1f

cat temp.sam | grep CCAGGCCAGGCCCTGCAGCAAACGGGCATT > temp1
cat temp.sam | grep TGGGGCCTGGAGGGGTGGGCAAGGGCTGGA >> temp1
cat temp.sam | grep GACATCAAGGATCTGGCCGCCCTCAGGGTG >> temp1
cat temp.sam | grep AGTTCTACATCATCCAGACCGGGGCGGACA >> temp1
cat temp.sam | grep CCACATCCACAC[GC]GTCACGCTGCCCCCTGC >> temp1
cat temp.sam | grep GGGGACAGCGGGAGGCCGGGCCAGGTGGGC >> temp1
cat temp.sam | grep CCCAGCCGGCCCCAGACCCGGCTCCACGCC >> temp1
cat temp.sam | grep CCCAGTGCACCTGCCGCCGCCATACCCGCT >> temp1
cat temp.sam | grep GGCCACAGCTTTCAAATCGTCCGCGATGAC >> temp1
cat temp.sam | grep CCAGGGTGACTCTGGAGGGCCCCTGGTCTG >> temp1
cat temp.sam | grep TGCAGGCGGGCGTGGTCAGCTGGGAGGAGA >> temp1

cat temp1 | wc | awk '{print $1}' > temp1i

# use `paste -s -` on mac
cat temp1[abcdefghi] | paste -s
rm temp[12]*

exit 0
