# TrackSig package workflow

# load CNA and purity dataframe (not loaded with VCF for parallelization memory saving)
# could be done as a single annotation load.... one function to load each file
# loads the following - all shared between all VCF's, all optional (but not necessarily independent)
# cna, purity, tumortypes, signatures (alex, cosmic), trinucleotide, sigactivities

load_annotation <- function(cna_file = "", purity_file = "", tumortype_file = "", signature_file = "", trinucleotide_file = "", active_signatures_file = ""){



}

# load VCF

#



# [END]
