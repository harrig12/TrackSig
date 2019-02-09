#!/bin/bash

# vcfToTsv.sh
# format vcf file for pyclone input
vcfFile=$1
cnaFile=$2
outFile=data/pyclone_sim940.tsv

# extract mutation_id
cut -f1,2 $vcfFile | tr "\t" "_" > tmp_col1
# rename in header
sed -i '' '1s/.*/mutation_id/' tmp_col1

# extract ref_counts
cut -f1 -d";" $vcfFile | cut -f2 -d"=" > tmp_col2
# rename in header
sed -i '' '1s/.*/ref_counts/' tmp_col2

# extract var_counts
cut -f3 -d"=" $vcfFile > tmp_col3
# rename in header
sed -i '' '1s/.*/var_counts/' tmp_col3

# extract normal_cn
cut -f3 -d"=" $vcfFile > tmp_col4
sed -i '' 's/.*/2/' tmp_col4
# rename in header
sed -i '' '1s/.*/normal_cn/' tmp_col4

# extract minor_cn
cut -f3 -d"=" $vcfFile > tmp_col5
sed -i '' 's/.*/0/' tmp_col5
# rename in header
sed -i '' '1s/.*/minor_cn/' tmp_col5

# extract major_cn
cut -f4 $cnaFile > tmp_col6
# rename in header
sed -i '' '1s/.*/major_cn/' tmp_col6

# combine columns
paste tmp_col1 tmp_col2 tmp_col3 tmp_col4 tmp_col5 tmp_col6 > $outFile

# clean up
rm tmp*

# [END]
