#################
# TrackSig Simulations
#################
library(TrackSig)

# generate simulations
#setwd("~/Desktop/TrackSig_simulations/")
getwd()

#set up
outdir <- "data"
#simulations <- c("sim5000")
#sim_muts <- list(c(2500, 1500, 1000))

simulations <- c("sim940")
sim_muts <- list(c(500, 300, 140))

# create data directory if missing
if (dir.exists(outdir) == FALSE) {
  dir.create(outdir)
}

# remove simulation active signatures (rebuilt upon simulation)
if (file.exists("annotation/sim_active_in_sample.txt")){
  unlink("annotation/sim_active_in_sample.txt")
}

# config (same variables as in in header.R)
TrackSig.options(purity_file = "annotation/sim_purity.txt",
                 signature_file = "annotation/sigProfiler_SBS_signatures.txt",
                 trinucleotide_file = "annotation/trinucleotide.txt",
                 active_signatures_file = "annotation/sim_active_in_sample.txt",
                 tumortype_file = "annotation/sim_tumortypes.txt",
                 sig_amount = "onlyKnownSignatures",
                 compute_bootstrap = FALSE,
                 cancer_type_signatures = FALSE)


#check setup
stopifnot(length(simulations) == length(sim_muts))


for (sim_i in 1:length(simulations)){
  # simulate vcfs
  print(sprintf("generating simulation %s...", simulations[sim_i]))
  create_simulation_set(with_CNA = F, outdir = outdir, simulation_name = simulations[sim_i],
                        n_mut_per_cluster = sim_muts[[sim_i]])
}


for (sim_i in 1:length(simulations)){
    # tracksig - make counts
  system(sprintf("src/run_simulations.sh data/mut_types/ %s", simulations[sim_i]))
}

# tracksig - compute mutational signatures
compute_signatures_for_all_examples(countsDir = "data/counts", bootstrapDir = "data/bootstrap/")

# tracksig - bootstrap not functional yet
#compute_errorbars_for_all_examples()




#################
# deconstructSig Simulations
#################
library(deconstructSigs)


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

fileName <- "data/vaf_DSformat_cosmic.txt"

DSformatSigs <- read.delim(fileName)

# bad habits
DSformatSigs <- as.data.frame(t(as.matrix(table(DSformatSigs))))
DSformatSigs[1,] <- lapply(DSformatSigs[1,], as.numeric)
rownames(DSformatSigs) <- "sim5000"


alex <- as.data.frame(t(TrackSig:::alex))
alexDS <- whichSignatures(tumor.ref = DSformatSigs, contexts.needed = T, signatures.ref = alex)

cosmic <- deconstructSigs::signatures.cosmic
cosmicDS <- whichSignatures(tumor.ref = DSformatSigs, contexts.needed = T, signatures.ref = cosmic)

round(alexDS$weights, 4)
round(cosmicDS$weights, 4)



