#!/bin/bash

# mutTypesToDS.sh
# format signature count tables for DeconstructSigs

fileToFormat=$1

cut -f4 $fileToFormat > alt
cut -f5 $fileToFormat > ref
cut -f6 $fileToFormat > tri

# === for compatibility with alex signatures ===
#paste -d"_" alt ref tri > outFile

#mv outFile data/vaf_DSformat_alex.txt

# === for DS internal reference signatures ====
paste -d">" <(paste -d[ <(cut -c1 tri) <(cut -c1 ref)) <(paste -d] <(cut -c2 tri) <(cut -c3 tri)) > outFile

# complement purines only
sed -i .bk '/^..A/s/C/g/g' outFile
sed -i .bk '/^..A/s/G/c/g' outFile
sed -i .bk '/^..A/s/T/a/g' outFile
sed -i .bk '/^..A/s/A/t/g' outFile
sed -i .bk '/^..G/s/C/g/g' outFile
sed -i .bk '/^..G/s/T/a/g' outFile
sed -i .bk '/^..G/s/A/t/g' outFile
sed -i .bk '/^..G/s/G/c/g' outFile
sed -i .bk '/./s/a/A/g' outFile
sed -i .bk '/./s/g/G/g' outFile
sed -i .bk '/./s/c/C/g' outFile
sed -i .bk '/./s/t/T/g' outFile

mv outFile data/vaf_DSformat_cosmic.txt

# === clean up ==============================
rm alt
rm ref
rm tri
rm *.bk



#[END]
