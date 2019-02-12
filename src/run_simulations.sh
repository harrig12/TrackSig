#!/bin/bash
# run TrackSig simulations from rpackage
# skips perl code

outdir=$1
sim_dir=$2
simulation_name=$3

if [[ -z $outdir ]]; then
	echo "Please provide an output directory ... exiting"
	exit
fi

if [ ! -d "$outdir" ]; then
  mkdir -p $outdir
fi

if [[ -z $simulation_name ]]; then
	echo "Please provide a simulation name ... exiting"
	exit
fi


# !!!!!!!!! TODO: also provide cna file




# make vaf
python src/make_corrected_vaf.py --vcf $sim_dir/"$simulation_name".vcf --output $sim_dir/"$simulation_name"_vaf.txt

# take out tri header
tail -n +2 $sim_dir/"$simulation_name"_tri.txt > tmp

# take out tri chr (interferes with sort) and sort
sed 's/^...//' tmp | sort -n > tmp_tri

# sort and extract phis
sort -n $sim_dir/"$simulation_name"_vaf.txt > tmp_vaf

# put together mutation_types file
sort -k 3 -r <(paste <(cut -f1,2 tmp_tri) <(cut -d= -f2 tmp_vaf) <(cut -f3,4,5 tmp_tri)) | cat > "$simulation_name".mut_types.txt

# restore chr prefix
ex -sc '%s/^/chr/|wq' "$simulation_name".mut_types.txt

# relocate to outdir
mv "$simulation_name".mut_types.txt $outdir

# clean up
rm tmp*

# make counts
src/sim_make_counts.sh $sim_dir/"$simulation_name".vcf $sim_dir/"$simulation_name"_vaf.txt

# compute mutational signtures
#Rscript src/compute_mutational_signatures.R

# [END]
