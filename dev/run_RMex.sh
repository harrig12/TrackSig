#!/bin/bash
# run track sig README example from clean install 

cd ~/Documents/BCB430/TrackSig/

#python src/make_corrected_vaf.py --vcf data/example.vcf --output data/example_vaf.txt
python src/make_corrected_vaf.py --vcf data/example.vcf --purity data/example_purity.txt --output data/example_vaf.txt

# uncomment to download and unzip hg19 
# rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ ./annotation/hg19/ 
# bash gunzip_hg19.c 

src/make_counts.sh data/example.vcf data/example_vaf.txt
Rscript src/compute_mutational_signatures.R
