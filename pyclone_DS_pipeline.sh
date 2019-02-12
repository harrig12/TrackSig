#!/bin/bash
# pyclone + deconstructSigs pipeline

pipelineDir=pyclone_DS
dataDir=data

if [ ! -d "$pipelineDir" ]; then
	mkdir -p $pipelineDir
fi

# for each simulation, find clusters with pyclone
for sim in data/*.vcf; do
    simName=

    # run pyclone on simulations
    # format for pyclone input
    bash vcfToTsv.sh data/$simName.vcf data/$simName_cna.txt $pipelineDir/pyclone_$simName.tsv

    # analysis
    PyClone build_mutations_file --in_file $pipelineDir/$simName.tsv --out_file $pipelineDir/$simName.yaml --prior total_copy_number
    PyClone setup_analysis --working_dir $pipelineDir --in_files $pipelineDir/$simName.tsv --init_method connected --prior total_copy_number
    PyClone run_analysis --config_file $pipelineDir/config.yaml

done

# for each cluster found by pyclone, find signatures with deconstructSigs

# assign mut_types clusters


# sort mut_types
# sort loci table
# merge





# [END]
