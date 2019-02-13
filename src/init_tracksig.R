# AUTHOR: Yulia Rubanova
# 2018

# suppressPackageStartupMessages(library("argparse"))
# parser <- ArgumentParser()

# parser$add_argument("--tumortypes", default="data/tumortypes.txt", help="file with cancer types of each sample")
# parser$add_argument("--purity", default="annotation/example_purity.txt", help="")
# parser$add_argument("--signatures", default="annotation/alexSignatures.txt", help="file with signatures definitions")
# parser$add_argument("--trinucleotide", default="annotation/trinucleotide.txt", help="")
# parser$add_argument("--active", default="annotation/active_in_samples.txt", help="specifies active signatures in TCGA cancer types")
# parser$add_argument("--pcawg", type=bool)
# args <- parser$parse_args()
# group <- args$group
# tumortype_file <- args$tumortypes
# purity_file = args$purity
# signature_file = args$signatures
# trinucleotide_file = args$trinucleotide_file
# active_signatures_file = args$active

options( expressions = 5e5 )

library(reshape2)
library(ggplot2)
library(NMF)

nmf.options(grid.patch=TRUE)

group = 0
# Select fitting all signatures or only signatures for the particular cancer type
sig_amount <- "onlyKnownSignatures" # recommended
#sig_amount <- "all" # not recommended, time-consuming

# if the signatures are specified per cancer type or per sample
cancer_type_signatures = FALSE

# if signatures trajectories need to be computed on bootstrapped signatures as well
# bootstrapping provides the uncertainty estimations on the trajectories
# warning: by default, mutations are bootstrapped 30 times and the script will run 30 time longer
compute_bootstrap = FALSE

sliding_window = FALSE
noise_sig = NULL

simulated_data = FALSE

postfix = ""

# specifies the changepoint detection algorithm.
changepoint_method = "PELT"

# folders with mutation counts, mutation order and bootstrapped mutations
# don't need to be changed unless different folder were specified in make_counts.sh
DIR_COUNTS = "data/counts/"
mutation_order = "data/mut_order/"
BOOTSTRAP_COUNTS = "data/bootstrap/"

# specifies active signatures in each sample. Contains the active signatures for the example
# active_signatures_file = "annotation/active_in_samples.txt"

SAVED_SAMPLES_DIR = "saved_data/"

PLOT_FULL_NAME = F
mutation_assignments = ""

if (!file.exists(DIR_RESULTS))
{
  dir.create(DIR_RESULTS, recursive=T)
}

src_files <- setdiff(grep(".*R$", list.files(paste0( "src"),full.names = T), value = T),
                     c(paste0( "src/compute_mutational_signatures.R"),
                       paste0( "src/header.R"),
                       paste0( "src/init_tracksig.R")))
for (file in src_files)
{
  source(file)
}

if (pcawg_format) {
  list[alex, tumortypes, active_signatures, active_signatures.our_samples] <-
      load_annotation_pcawg(tumortype_file, signature_file, active_signatures_file, trinucleotide_file)
} else {
  list[alex, tumortypes, active_signatures, active_signatures.our_samples] <-
    load_annotation(tumortype_file, signature_file, active_signatures_file)
}
