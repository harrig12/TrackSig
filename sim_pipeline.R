#################
# TrackSig Simulations
#################
library(TrackSig)
source("~/Documents/Cait-TrackSig/sim_pipeline_helpers.R")

# generate simulations
#setwd("~/Desktop/TrackSig_simulations/")
getwd()

#set up
outdir <- "data"

# create data directory
dir.create(outdir, showWarnings=FALSE)

#generate simulations
####################################################
#annotation files
sim_activities_file <- "annotation/sim_active_in_sample.txt"
sim_purity_file <- "annotation/sim_purity.txt"
sim_tumortype_file <-"annotation/sim_tumortypes.txt"
signature_file <- "annotation/sigProfiler_SBS_signatures.txt"

mut_per_sim <- 5000

tracksig_results_dir = "TS_results_signature_trajectories/"


# Create simulation vcfs
print(sprintf("Generating simulations ..."))
simulations = create_simulation_set(outdir = outdir,
                                    mut_per_sim = mut_per_sim,
                                    sim_activity_file = sim_activities_file,
                                    sim_purity_file = sim_purity_file,
                                    sim_tumortype_file = sim_tumortype_file,
                                    signature_file = signature_file)

####################################################

simulations <- list.files(outdir)
sel <- grep(x = simulations, "^Simulation")
simulations <- simulations[sel]

# simulations by bin size
bin_sizes <- c(25, 50, 100, 125, 250)

for (size_i in seq_along(bin_sizes)){
  # TrackSig - set options (same variables as in in header.R)
  TrackSig.options(purity_file = sim_purity_file,
                   signature_file = "annotation/sigProfiler_SBS_signatures.txt",
                   trinucleotide_file = "annotation/trinucleotide.txt",
                   active_signatures_file = sim_activities_file,
                   tumortype_file = sim_tumortype_file,
                   sig_amount = "onlyKnownSignatures",
                   compute_bootstrap = FALSE,
                   cancer_type_signatures = FALSE,
                   pcawg_format = TRUE,
                   DIR_RESULTS = tracksig_results_dir,
                   bin_size = bin_sizes[size_i])


  # tracksig - make counts
  for (sim_i in 1:length(simulations)){
    print(sprintf("%s", simulations[sim_i]))
    run_simulation(simulations[sim_i])

  }

  # tracksig - compute mutational signatures
  compute_signatures_for_all_examples(countsDir = "data/counts", bootstrapDir = "data/bootstrap/")

  # tracksig - get exposures
  extract_exposures_per_mutation(activities_dir = paste0(tracksig_results_dir, "/SIMULATED/"),
                                 sorted_mutations_dir = "data/mut_types/", bin_size = bin_sizes[size_i])

  # tracksig - bootstrap not functional yet
  #compute_errorbars_for_all_examples()

  res <- compare_simulation_results(simulations,
      ground_truth_dir = outdir,
      method_results_dir = paste0(tracksig_results_dir, "/SIMULATED/"),
      res_file_name = sprintf("TrackSig_simulation_results_post%d.txt", bin_sizes[size_i]))

  pdf(sprintf("TrackSig_KL_post%d.pdf", bin_sizes[size_i]), width = 5, height=5)

  # ploting arguments
  res$pch <- rep_len(c(20, 5), length.out = dim(res)[1])
  res$col <- rep_len(c(2,2,3,3,4,4,5,5,6,6,7,7,8,8), length.out = dim(res)[1])

  plot(res$kl, res$abs_diff_max, main=sprintf("TrackSig KL post_bin_size = %d", bin_sizes[size_i]),
     xlab="KL", ylab="max abs diff", pch = res$pch, col = res$col)

  pre_bin_sizes <- round(seq(30, 300, length.out = 8), 0)

  legend("bottomright", pch = c(20, 5, rep(15, 8)),
         col = c(1, 1, 2:8), legend = c("depth 100", "depth 1000", as.character(pre_bin_sizes)))

  res$labels <- strsplit(as.character(res$sim), "^Simulation_")
  res$labels <- apply(res["labels"], MARGIN = 1, FUN=unlist)[2,]
  res$labels <- strsplit(res$labels, "[[:digit:]]_[[:alnum:]]*$")
  res$labels <- apply(res["labels"], MARGIN = 1, FUN=unlist)

  text(res$kl, res$abs_diff_max, labels = res$labels, cex = 0.8, pos = 4)

  dev.off()
}
