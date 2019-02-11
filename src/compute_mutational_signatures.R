# AUTHOR: Yulia Rubanova

# Run TrackSig with the data specified in header.R
source("src/header.R")

group = 0
EXAMPLES_PER_GROUP <- 500

#save_data_for_samples()
suppressMessages(compute_signatures_for_all_examples())
if (compute_bootstrap) {compute_errorbars_for_all_examples()}

