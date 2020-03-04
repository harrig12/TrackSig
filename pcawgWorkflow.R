c# pcawgWorkflow.R

# Code to use TrackSig R package with real pcawg counts data.
# Author: Cait Harrigan
# June 2019

# DEFS ##############################################################

# conveinent list stuct
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}


# SIMULATION WORKFLOW ##################################################

reticulate::use_condaenv("tracksig")
library(TrackSig)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(foreach)
library(doParallel)
library(ggplot2)

# Remember: unless registerDoMC is called, foreach will not run in parallel. Simply loading the doParallel package is not enough
registerDoParallel(cores=20)

setwd("~/Desktop/pcawg_sim/")
#TrackSig:::create_simulation_set(outdir = path.expand("~/Desktop/pcawg_sim/simulation_data/"))

# set up
TrackSig.options(purity_file = "~/Desktop/pcawg_sim/annotation/sim_purity.txt",
                 signature_file = "~/Desktop/pcawg_sim/annotation/sigProfiler_SBS_signatures.txt",
                 trinucleotide_file = "~/Desktop/pcawg_sim/annotation/trinucleotide.txt",
                 active_signatures_file = "~/Desktop/pcawg_sim/annotation/sim_active_in_sample.txt",
                 tumortype_file = "~/Desktop/pcawg_sim/annotation/sim_tumortypes.txt",
                 sig_amount = "onlyKnownSignatures",
                 compute_bootstrap = FALSE,
                 cancer_type_signatures = FALSE,
                 pcawg_format = TRUE,
                 DIR_RESULTS = "og_results/",
                 pelt_penalty = expression( (n_sigs + 1) * log(n_bins)),
                 pelt_score_fxn = TrackSig:::sum_gaussian_mixture_multinomials_ll,
                 bin_size = 100)

# load annotation
list[alex, tumortypes, active_signatures, active_signatures.our_samples] <- TrackSig:::load_annotation_pcawg()

# get sim names

simnames <- list.files("~/Desktop/pcawg_sim/simulation_data/")
#sel <- grepl("100$", simnames)
#simnames <- simnames[sel]

foreach (i=1:length(simnames)) %dopar% {
#foreach (i=1:1) %dopar% {

  simname <-simnames[i]

  loadAndScoreIt_pcawg(vcfFile = paste0("~/Desktop/pcawg_sim/simulation_data/", simname, "/", simname, ".vcf"),
                       cnaFile = paste0("~/Desktop/pcawg_sim/simulation_data/", simname, "/", simname, "_cna.txt"),
                       #purityFile = "~/Desktop/pcawg/annotation/sim_purity.txt",
                       tumortypes = tumortypes, acronym = "SIMULATED")

}

# SIMULATION COMPARISON #############################

# want to compare quality of simulaitons.
# need ground truth information, as well as scoring for each simulation.

getSimMeta <- function(dataDir){
  # pickup simulation info

  simUnparsed <- list.files(dataDir)
  simMeta <-data.frame(row.names = simUnparsed)

  simUnparsed <- strsplit(simUnparsed, "_")
  simUnparsed <- t(as.data.frame(simUnparsed, stringsAsFactors = F))
  rownames(simUnparsed) <- NULL
  simUnparsed[,1] <- paste0(simUnparsed[,1], simUnparsed[,2], simUnparsed[,3])
  simUnparsed <- simUnparsed[,c(1, 4:7)]
  colnames(simUnparsed) <- c("type", "depth", "bin", "subCloneAt", "sigChange")
  simUnparsed <- data.frame(simUnparsed, stringsAsFactors = F)

  # populate the meta df
  simMeta$type <- simUnparsed$type
  simMeta$depth <- as.numeric(unlist(strsplit(simUnparsed$depth, "depth"))[c(F,T)])
  simMeta$bin <- as.numeric(unlist(strsplit(simUnparsed$bin, "bin"))[c(F,T)])
  simMeta$subCloneAt <- as.numeric(unlist(strsplit(simUnparsed$subCloneAt, "subCloneAt"))[c(F,T)])
  simMeta$sigChange <- as.numeric(unlist(strsplit(simUnparsed$sigChange, "sigChange"))[c(F,T)])

  return(simMeta)
}


getSimCPs <- function(resultDir, simMeta){
  # pick up the found changepoints and their associated phi

  #simNames <- rownames(simMeta)
  simNames <- list.files(resultDir)
  CPs <- list()

  for (simName in simNames) {

    # try to read cps, handle if file empty
    simCPs <- tryCatch({read.table(paste0(resultDir, "/", simName, "/changepoints.txt"))},
                    error = function(err){ return("*") })

    CPs[[simName]] <- simCPs
  }

  return(CPs)
}


simMeta <- getSimMeta("~/Desktop/pcawg_sim/simulation_data/")

# ground truth cps's placed exactly between 1.0 and dist. Only fair if
# variance is equal.
gtCPs <- (simMeta$dist + 1) / 2

simCPs <- getSimCPs("~/Desktop/sigAddSimulation_results/", simMeta)

correct <- lapply(simCPs, length) == 1
names(correct) <- NULL

simMeta$correctNcp <- correct

# BOLTZ PCAWG WORKFLOW ##################################################

reticulate::use_condaenv("tracksig")
library(TrackSig)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(foreach)
library(doParallel)

# Remember: unless registerDoMC is called, foreach will not run in parallel. Simply loading the doParallel package is not enough
registerDoParallel(cores=20)

setwd("~/Desktop/small25/")
#TrackSig:::create_simulation_set(outdir = "~/Desktop/pcawg/simulation_data")

# set up
TrackSig.options(#purity_file = "~/Desktop/pcawg/annotation/sim_purity.txt",
                 #signature_file = "~/Desktop/pcawg/annotation/sigProfiler_SBS_signatures.txt",
                 #trinucleotide_file = "~/Desktop/pcawg/annotation/trinucleotide.txt",
                 #active_signatures_file = "~/Desktop/pcawg/annotation/sim_active_in_sample.txt",
                 #tumortype_file = "~/Desktop/pcawg/annotation/sim_tumortypes.txt",
                 sig_amount = "onlyKnownSignatures",
                 compute_bootstrap = FALSE,
                 cancer_type_signatures = FALSE,
                 pcawg_format = TRUE,
                 DIR_RESULTS = "results/",
                 pelt_penalty = expression( (n_sigs + 1) * log(n_bins)),
                 pelt_score_fxn = TrackSig:::sum_gaussian_mixture_multinomials_ll,
                 bin_size = 100)

# load annotation
load("~/Desktop/Yulia_current0519/annotation_data.pcawg_sigs.sigProfiler.RData")

simnames <- list.files("./data")


foreach (i=1:length(simnames)) %dopar% {
#foreach (i=1:1) %dopar% {

  simname <-simnames[i]

  loadAndScoreIt_pcawg(vcfFile = paste0("~/Desktop/small25/data/", simname),
                       #cnaFile = "~/Desktop/pcawg/annotation/example_cna.txt",
                       #purityFile = "~/Desktop/pcawg/annotation/sim_purity.txt",
                       tumortypes = tumortypes, acronym = "SIMULATED")

}



# [END] ####
