#!/bin/bash
# run track sig README example from clean install 

cd ~/Documents/Cait-TrackSig

#python src/make_corrected_vaf.py --vcf data/example.vcf --output data/example_vaf.txt
python src/make_corrected_vaf.py --vcf data/example.vcf --purity data/example_purity.txt --output data/example_vaf.txt

# uncomment to download and unzip hg19 
# rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ ./annotation/hg19/ 
# bash dev/gunzip_hg19.sh 

src/make_counts.sh data/example.vcf data/example_vaf.txt
Rscript src/compute_mutational_signatures.R
