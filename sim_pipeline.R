#################
# TrackSig Simulations
#################
#library(TrackSig)

# generate simulations
#setwd("~/Desktop/TrackSig_simulations/")
getwd()

#set up
outdir <- "data"
#simulations <- c("sim940", "sim5000", "sim9500")
#sim_muts <- list(c(500, 300, 140), c(2500, 1500, 1000), c(5000, 3000, 1500))

mut_per_sim = 5000

#annotation files
sim_activities_file <- "annotation/sim_active_in_sample.txt"
sim_purity_file <- "annotation/sim_purity.txt"
sim_tumortype_file <-"annotation/sim_tumortypes.txt"
signature_file = "annotation/sigProfiler_SBS_signatures.txt"
trinucleotide_file = "annotation/trinucleotide.txt"

# create data directory if missing
# if (dir.exists(outdir) == FALSE) {
#   dir.create(outdir)
# }
dir.create(outdir, showWarnings=FALSE)

source("src/ccf_simulations.R")

####################################################
# Create simulation vcfs
print(sprintf("Generating simulations ..."))
simulations = create_simulation_set(outdir = outdir, 
  mut_per_sim = mut_per_sim,
  sim_activity_file = sim_activities_file,
  sim_purity_file = sim_purity_file,
  sim_tumortype_file = sim_tumortype_file,
  signature_file = signature_file,
  trinucleotide_file = trinucleotide_file)

####################################################
# Run simulations through TrackSig
# library("TrackSig")

# # config (same variables as in in header.R)
# TrackSig.options(purity_file = sim_purity_file,
#                  signature_file = signature_file,
#                  trinucleotide_file = trinucleotide_file,
#                  active_signatures_file = sim_activities_file,
#                  tumortype_file = sim_tumortype_file,
#                  sig_amount = "onlyKnownSignatures",
#                  compute_bootstrap = FALSE,
#                  cancer_type_signatures = FALSE)

for (sim_i in 1:length(simulations)){
  # tracksig - make counts
  print(sprintf("%s", simulations[sim_i]))
  system(sprintf("src/run_simulations.sh data/mut_types/ data/%s %s", simulations[sim_i], simulations[sim_i]))
}

pcawg_format = T
tumortype_file <- sim_tumortype_file
purity_file <- sim_purity_file
active_signatures_file <- sim_activities_file
source("src/init_tracksig.R")

# source(sprintf("src/init_tracksig.R --tumortypes %s --purity %s  --signatures %s --trinucleotide %s --active %s", 
#   sim_tumortype_file, sim_purity_file, signature_file, trinucleotide_file, sim_activities_file))

# tracksig - compute mutational signatures
compute_signatures_for_all_examples(countsDir = "data/counts/", bootstrapDir = "data/bootstrap/")

# tracksig - bootstrap not functional yet
#compute_errorbars_for_all_examples()




#################
# deconstructSig Simulations
#################
#library(deconstructSigs)


#vcfData <- read.delim(sprintf("data/%s.vcf", simulation_name))
#
## format for deconstructSigs according to https://github.com/raerose01/deconstructSigs
#vcfData <- cbind(1, vcfData)
#colnames(vcfData)[1] <- "ids"
#
#
#sigs <- mut.to.sigs.input(vcfData[1:100,],
#                          sample.id = "ids",
#                          chr = "X.CHROM",
#                          pos = "POS",
#                          ref = "REF",
#                          alt = "ALT")
#

# fileName <- "data/vaf_DSformat_cosmic.txt"

# DSformatSigs <- read.delim(fileName)

# # bad habits
# DSformatSigs <- as.data.frame(t(as.matrix(table(DSformatSigs))))
# DSformatSigs[1,] <- lapply(DSformatSigs[1,], as.numeric)
# rownames(DSformatSigs) <- "sim5000"


# alex <- as.data.frame(t(TrackSig:::alex))
# alexDS <- whichSignatures(tumor.ref = DSformatSigs, contexts.needed = T, signatures.ref = alex)

# cosmic <- deconstructSigs::signatures.cosmic
# cosmicDS <- whichSignatures(tumor.ref = DSformatSigs, contexts.needed = T, signatures.ref = cosmic)

# round(alexDS$weights, 4)
# round(cosmicDS$weights, 4)



