# pcawgWorkflow.R

# Code to use TrackSig R package with real pcawg counts data.
# Author: Cait Harrigan
# June 2019

library(TrackSig)
library(foreach)
library(doParallel)

# Remember: unless registerDoMC is called, foreach will not run in parallel. Simply loading the doParallel package is not enough
registerDoParallel(cores=2)

# load the provided annotation
load("~/Desktop/pcawg/annotation/pcawg_annotation.RData")

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

TrackSig.options(purity_file = "~/Desktop/pcawg/data/consensus_purity.txt",
                 DIR_RESULTS = "~/Desktop/pcawg/results/")

countsDir <- "~/Desktop/pcawg/data/counts"

loadAndScoreIt_pcawg <- function(samplename, countsDir = "~/Desktop/pcawg/data/counts", tumortypes) {


  # load_sample() returns (inorder)
  # list(vcfData, vcf, phis, quadratic_phis, phis_sliding_window, assigns_phylo_nodes,
  #      assigns_phylo_nodes_sw, acronym, window, tumor_id, phis_for_plot, bootstrap_vcfs, bootstrap_phis))

  list[vcfData, vcf, phis, quadratic_phis, phis_sliding_window, assigns_phylo_nodes,
       assigns_phylo_nodes_sw, acronym, window, tumor_id, phis_for_plot,
       bootstrap_vcfs, bootstrap_phis] <- load_sample(samplename, countsDir, bootstrapDir = NULL, tumortypes)

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

  if (sum(vcf[apply(alex.t,1,sum) == 0,] != 0) != 0) {
    print(paste0("Sample ", example, ": some trinucleotides have probability 0 under the model, but their count is non-zero. Sot the count vector is impossible under the model."))
    next
  }

  # this is the work-horse part

  dir_name <- paste0(TrackSig.options()$DIR_RESULTS, acronym, "/", tumor_id, "/")

  suppressWarnings(dir.create(dir_name, recursive = T))

  if (!is.null(phis_for_plot))
  {
    write(phis_for_plot, file=paste0(dir_name, "phis.txt"), ncolumns=length(phis_for_plot))
  }

  method_name <- "iterativeChangePoints"

  if (!file.exists(paste0(dir_name, "mixtures.csv")) || !file.exists(paste0(dir_name, "changepoints.txt")))
  {
    if (TrackSig.options()$changepoint_method == "PELT") {
      list[changepoints, mixtures] <- TrackSig:::find_changepoints_pelt(vcf, alex.t, phis, quadratic_phis)
    } else {
      list[bics, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(vcf, alex.t, n_signatures = ncol(alex.t))
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

  if (!is.null(assigns_phylo_nodes_sw)) {
    write(assigns_phylo_nodes_sw,  file=paste0(dir_name, "assignments.txt"), ncolumns=length(assigns_phylo_nodes_sw))
  } else  {
    n_clusters = transition_points = assigns_phylo_nodes_sw = NULL
  }

  plot_name <- paste0(dir_name, "/", acronym, "_", tumor_id, "_", TrackSig.options()$sig_amount, TrackSig.options()$postfix, ".pdf")

  if (TrackSig.options()$PLOT_FULL_NAME)
  {
    plot_name <- paste0(dir_name, "/", acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_noPrior_", method_name, postfix, ".pdf")
  }

  mark_cp <- !is.null(changepoints)
  TrackSig:::plot_signatures(mixtures*100, plot_name=plot_name, phis = phis_for_plot, mark_change_points=mark_cp, change_points=changepoints,
                  #assigns_phylo_nodes = assigns_phylo_nodes_sw,
                  transition_points = transition_points,
                  scale=1.2)

  mixtures.rescaled = NULL
}





# several samples in parallel
sel <- grep("([^/]*)\\.phi\\.txt", list.files(countsDir))
sampleNames <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(countsDir)[sel])

foreach (i=1:length(sampleNames)) %dopar% {
  sample = sampleNames[i]
  loadAndScoreIt_pcawg(sample)
}




# SIMULATION WORKFLOW

library(TrackSig)
library(foreach)
library(doParallel)

# Remember: unless registerDoMC is called, foreach will not run in parallel. Simply loading the doParallel package is not enough
registerDoParallel(cores=10)

# setwd("~/Desktop/pcawg/)
# TrackSig:::create_simulation_set()

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
                 pelt_penalty = expression(0),
                 pelt_score_fxn = TrackSig:::gaussian_ll)

# load annotation
list[alex, tumortypes, active_signatures, active_signatures.our_samples] <- TrackSig:::load_annotation_pcawg()

# get sim names

#simnames <- list.files("~/Desktop/pcawg/simulations_data/")
#foreach (i=1:length(simnames)) %dopar% {
#  simname = simnames[i]
#  run_simulation(simname, "simulations_data")
#  loadAndScoreIt_pcawg(simname, "~/Desktop/pcawg/simulations_data/counts/")
#}


#run_simulation(simnames[5], "simulations_data")
loadAndScoreIt_pcawg(simnames[5], "~/Desktop/pcawg/simulations_data/counts/", tumortypes)



#a <- as.numeric(read.delim("~/Desktop/pcawg/simulation_results/SIMULATED/neg_LL/phis.txt", header = F, sep = " "))
#length(a)




# [END]
