#!/bin/bash
#gunzip hg19 tarballs for TrackSig
for i in `ls ./annotation/hg19/* | grep "chr\d\(\d\)\?.fa.gz"`;
do 
	gunzip $i 
done
gunzip ./annotation/hg19/chrX.fa
gunzip ./annotation/hg19/chrY.fa

