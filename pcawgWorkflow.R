# pcawgWorkflow.R

# Code to use TrackSig R package with real pcawg counts data.
# Author: Cait Harrigan
# June 2019

#library(TrackSig)
#library(foreach)
#library(doParallel)

# Remember: unless registerDoMC is called, foreach will not run in parallel. Simply loading the doParallel package is not enough
#registerDoParallel(cores=2)

# load the provided annotation
#load("~/Desktop/pcawg/annotation/pcawg_annotation.RData")

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


# first attempt: one sample -> function
# TrackSig::: calles added

#TrackSig.options(purity_file = "~/Desktop/pcawg/data/consensus_purity.txt",
#                 DIR_RESULTS = "~/Desktop/pcawg/results/")

#countsDir <- "~/Desktop/pcawg/data/counts"

loadAndScoreIt_pcawg <- function(vcfFile, tumortypes, acronym) {


  # load_sample() returns (inorder)
  # list(vcfData, vcf, phis, quadratic_phis, phis_sliding_window, assigns_phylo_nodes,
  #      assigns_phylo_nodes_sw, acronym, window, tumor_id, phis_for_plot, bootstrap_vcfs, bootstrap_phis))

  vcfFile <- path.expand(vcfFile)
  tumor_id <- strsplit( unlist(strsplit(vcfFile, "/"))[ length( strsplit(vcfFile, "/")[[1]] ) ] , ".vcf")[[1]]

  list[phis, quadPhis, counts] <- vcfToCounts(vcfFile, refGenome = BSgenome.Hsapiens.UCSC.hg19)

  # following checking is all from within compute_signatures_for_all_examples()
  # will throw a next error if check fails

  if (TrackSig.options()$sig_amount == "onlyKnownSignatures") {
    # Fit only known signatures
    list[alex.t, matched_type, acronym] <- TrackSig:::get_signatures_for_current_sample(tumor_id, active_signatures.our_samples, alex, TrackSig.options()$noise_sig)
  } else {
    alex.t <- alex
  }

  if (is.null(acronym) || acronym == "") {
    print(paste("ERROR: Cancer type not found for ", example))
    next
  }

  if (is.null(alex.t))
  {
    print(paste0("No active signatures for sample", example, " ...."))
    next
  }

  if (is.vector(alex.t)) {
    next
  }

  if (sum(counts[apply(alex.t,1,sum) == 0,] != 0) != 0) {
    print(paste0("Sample ", example, ": some trinucleotides have probability 0 under the model, but their count is non-zero. Sot the count vector is impossible under the model."))
    next
  }

  # this is the work-horse part

  dir_name <- paste0(TrackSig.options()$DIR_RESULTS, acronym, "/", tumor_id, "/")
  suppressWarnings(dir.create(dir_name, recursive = T))

  #method_name <- "iterativeChangePoints"

  if (!file.exists(paste0(dir_name, "mixtures.csv")) || !file.exists(paste0(dir_name, "changepoints.txt")))
  {
    if (TrackSig.options()$changepoint_method == "PELT") {
      list[changepoints, mixtures] <- TrackSig:::find_changepoints_pelt(counts, alex.t, phis, quadPhis)
    } else {
      list[bics, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(counts, alex.t, n_signatures = ncol(alex.t))
    }

    write.csv(mixtures, file=paste0(dir_name, "mixtures.csv"))

    n_col <- ifelse(length(changepoints) > 0, length(changepoints), 1)
    write(changepoints, file=paste0(dir_name, "changepoints.txt"), ncolumns=n_col)
  } else {
    mixtures <- read_mixtures(paste0(dir_name, "mixtures.csv"))
    cp_file = paste0(dir_name, "changepoints.txt")
    if (file.info(cp_file)$size == 1) {
      changepoints <- c()
    } else {
      changepoints <- unlist(read.table(cp_file, header=F))
    }
  }


  plot_name <- paste0(dir_name, "/", acronym, "_", tumor_id, "_", TrackSig.options()$sig_amount, TrackSig.options()$postfix, ".pdf")

  #if (TrackSig.options()$PLOT_FULL_NAME)
  #{
  #  plot_name <- paste0(dir_name, "/", acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_noPrior_", method_name, postfix, ".pdf")
  #}

  mark_cp <- !is.null(changepoints)
  TrackSig:::plot_signatures(mixtures*100, plot_name=plot_name, phis = phis, mark_change_points=mark_cp, change_points=changepoints,
                  transition_points = NULL,
                  scale=1.2)

  mixtures.rescaled = NULL
}

loadAndScoreIt_simulation <- function(vcfFile, countsDir, tumortypes, acronym = "SIMUlATION") {


  # load_sample() returns (inorder)
  # list(vcfData, vcf, phis, quadratic_phis, phis_sliding_window, assigns_phylo_nodes,
  #      assigns_phylo_nodes_sw, acronym, window, tumor_id, phis_for_plot, bootstrap_vcfs, bootstrap_phis))

  vcfFile <- path.expand(vcfFile)
  tumor_id <- strsplit( unlist(strsplit(vcfFile, "/"))[ length( strsplit(vcfFile, "/")[[1]] ) ] , ".vcf")[[1]]

  # list[phis, quadPhis, counts] <- vcfToCounts(vcfFile, refGenome = BSgenome.Hsapiens.UCSC.hg19)
  list[phis, quadPhis, counts] <- load_sample(tumor_id, countsDir, NULL, tumortypes)[c(3, 4, 2)]

  # following checking is all from within compute_signatures_for_all_examples()
  # will throw a next error if check fails

  if (TrackSig.options()$sig_amount == "onlyKnownSignatures") {
    # Fit only known signatures
    list[alex.t, matched_type, acronym] <- TrackSig:::get_signatures_for_current_sample(tumor_id, active_signatures.our_samples, alex, TrackSig.options()$noise_sig)
  } else {
    alex.t <- alex
  }

  if (is.null(acronym) || acronym == "") {
    print(paste("ERROR: Cancer type not found for ", example))
    next
  }

  if (is.null(alex.t))
  {
    print(paste0("No active signatures for sample", example, " ...."))
    next
  }

  if (is.vector(alex.t)) {
    next
  }

  if (sum(counts[apply(alex.t,1,sum) == 0,] != 0) != 0) {
    print(paste0("Sample ", example, ": some trinucleotides have probability 0 under the model, but their count is non-zero. Sot the count vector is impossible under the model."))
    next
  }

  # this is the work-horse part

  dir_name <- paste0(TrackSig.options()$DIR_RESULTS, acronym, "/", tumor_id, "/")
  suppressWarnings(dir.create(dir_name, recursive = T))

  #method_name <- "iterativeChangePoints"

  if (!file.exists(paste0(dir_name, "mixtures.csv")) || !file.exists(paste0(dir_name, "changepoints.txt")))
  {
    if (TrackSig.options()$changepoint_method == "PELT") {
      list[changepoints, mixtures] <- TrackSig:::find_changepoints_pelt(counts, alex.t, phis, quadPhis)
    } else {
      list[bics, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(counts, alex.t, n_signatures = ncol(alex.t))
    }

    write.csv(mixtures, file=paste0(dir_name, "mixtures.csv"))

    n_col <- ifelse(length(changepoints) > 0, length(changepoints), 1)
    write(changepoints, file=paste0(dir_name, "changepoints.txt"), ncolumns=n_col)
  } else {
    mixtures <- read_mixtures(paste0(dir_name, "mixtures.csv"))
    cp_file = paste0(dir_name, "changepoints.txt")
    if (file.info(cp_file)$size == 1) {
      changepoints <- c()
    } else {
      changepoints <- unlist(read.table(cp_file, header=F))
    }
  }


  plot_name <- paste0(dir_name, "/", acronym, "_", tumor_id, "_", TrackSig.options()$sig_amount, TrackSig.options()$postfix, ".pdf")

  #if (TrackSig.options()$PLOT_FULL_NAME)
  #{
  #  plot_name <- paste0(dir_name, "/", acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_noPrior_", method_name, postfix, ".pdf")
  #}

  mark_cp <- !is.null(changepoints)
  TrackSig:::plot_signatures(mixtures*100, plot_name=plot_name, phis = phis, mark_change_points=mark_cp, change_points=changepoints,
                             transition_points = NULL,
                             scale=1.2)

  mixtures.rescaled = NULL
}




## several samples in parallel
#sel <- grep("([^/]*)\\.phi\\.txt", list.files(countsDir))
#sampleNames <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(countsDir)[sel])
#
#foreach (i=1:length(sampleNames)) %dopar% {
#  sample = sampleNames[i]
#  loadAndScoreIt_pcawg(sample)
#}



##################################################
# SIMULATION WORKFLOW
#################################################

library(TrackSig)
library(BSgenome.Hsapiens.UCSC.hg19)
library(foreach)
library(doParallel)

# Remember: unless registerDoMC is called, foreach will not run in parallel. Simply loading the doParallel package is not enough
registerDoParallel(cores=10)

#setwd("~/Desktop/pcawg/")
#TrackSig:::create_simulation_set(outdir = "~/Desktop/pcawg/simulations_data")

# set up
TrackSig.options(purity_file = "~/Desktop/pcawg/annotation/sim_purity.txt",
                 signature_file = "~/Desktop/pcawg/annotation/sigProfiler_SBS_signatures.txt",
                 trinucleotide_file = "~/Desktop/pcawg/annotation/trinucleotide.txt",
                 active_signatures_file = "~/Desktop/pcawg/annotation/sim_active_in_sample.txt",
                 tumortype_file = "~/Desktop/pcawg/annotation/sim_tumortypes.txt",
                 sig_amount = "onlyKnownSignatures",
                 compute_bootstrap = FALSE,
                 cancer_type_signatures = FALSE,
                 pcawg_format = TRUE,
                 DIR_RESULTS = "simulation_results/",
                 pelt_penalty = expression( (n_sigs - 1) * log(n_bins) + log(n_bins) ),
                 pelt_score_fxn = TrackSig:::sum_gaussian_mixture_multinomials_ll,
                 bin_size = 100)

# load annotation
list[alex, tumortypes, active_signatures, active_signatures.our_samples] <- TrackSig:::load_annotation_pcawg()

# get sim names

simnames <- list.files("~/Desktop/pcawg/simulations_data/")
sel <- grepl("100$", simnames)
simnames <- simnames[sel]

#foreach (i=1:length(simnames)) %dopar% {
#  simname = simnames[i]
#  run_simulation(simname, "simulations_data")
#  loadAndScoreIt_pcawg(simname, "~/Desktop/pcawg/simulations_data/counts/", tumortypes)
#}

loadAndScoreIt_simulation(vcfFile = paste0("~/Desktop/pcawg/simulations_data/", simnames[1], "/", simnames[1], ".vcf"),
                          countsDir = "~/Desktop/pcawg/simulations_data/counts/", tumortypes = tumortypes)

#run_simulation(simnames[1], "simulations_data")
#loadAndScoreIt_pcawg(simnames[1], "~/Desktop/pcawg/simulations_data/counts/", tumortypes)


## compare cp's found
resDir <- "~/Desktop/pcawg/simulation_results/"
sim_types <- c()
sim_cps <- c()

for (sim_dir in list.files("~/Desktop/pcawg/simulation_results/")){
  sim_types <- c(sim_types, sim_dir)
  print(sim_types)
  cps <- c()

  for (sim in list.files(paste0("~/Desktop/pcawg/simulation_results/", sim_dir))){

    print(sim)

    path <- paste0(resDir, sim_dir, "/", sim, "/changepoints.txt")
    cp <- NA
    tryCatch(
      cp <- read.delim(path, header = F, sep = " ", blank.lines.skip = T, skipNul = T)$V1
      , error = function(e) NULL
    )

    cps <- c(cps, cp)

  }

  sim_cps[[sim_dir]] <- cps
}

#vcfFile <- "~/Desktop/pcawg/data/Simulation_two_clusters4_depth100_bin100/Simulation_two_clusters4_depth100_bin100.vcf"

##################################################
# EXAMPLE.vcf
#################################################


# set up

library(TrackSig)
library(BSgenome.Hsapiens.UCSC.hg19)

TrackSig.options(purity_file = system.file("extdata", "example_purity.txt", package = "TrackSig"),
                 signature_file = "~/Desktop/pcawg/annotation/sigProfiler_SBS_signatures.txt",
                 #trinucleotide_file = "~/Desktop/pcawg/annotation/trinucleotide.txt",
                 #active_signatures_file = "~/Desktop/pcawg/annotation/sim_active_in_sample.txt",
                 tumortype_file = system.file("extdata", "tumortypes.txt", package = "TrackSig"),
                 sig_amount = "onlyKnownSignatures",
                 compute_bootstrap = FALSE,
                 cancer_type_signatures = FALSE,
                 pcawg_format = TRUE,
                 DIR_RESULTS = "simulation_results/",
                 pelt_penalty = expression( (n_sigs - 1) * log(n_bins) + log(n_bins) ),
                 pelt_score_fxn = TrackSig:::sum_gaussian_mixture_multinomials_ll,
                 bin_size = 100)

load("~/Desktop/pcawg/annotation/pcawg_annotation.Rdata")
active_signatures.our_samples$ID[1] <- "example"
tumortypes <- data.frame(ID = "example", tumor_type = "SIMULATED")

# use pcawg funtion

vcfFile <- system.file("extdata", "example.vcf", package = "TrackSig")
loadAndScoreIt_pcawg(vcfFile, tumortypes = tumortypes, acronym = "SIMULATED")


# [END]
