#################
# TrackSig Simulations
#################
library(TrackSig)
source("sim_pipeline_helpers.R")

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
bin_sizes <- c(30, 75, 100, 156, 213, 250, 300)

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
                 bin_size = 30)


# tracksig - make counts
for (sim_i in 1:length(simulations)){
  print(sprintf("%s", simulations[sim_i]))
  run_simulation(simulations[sim_i])

}

# tracksig - compute mutational signatures
compute_signatures_for_all_examples(countsDir = "data/counts", bootstrapDir = "data/bootstrap/")

# tracksig - get exposures
extract_exposures_per_mutation(activities_dir = paste0(tracksig_results_dir, "/SIMULATED/"),
                               sorted_mutations_dir = "data/mut_types/", bin_size = 100)

# tracksig - bootstrap not functional yet
#compute_errorbars_for_all_examples()

res <- compare_simulation_results(simulations,
    ground_truth_dir = outdir,
    method_results_dir = paste0(tracksig_results_dir, "/SIMULATED/"),
    res_file_name = "TrackSig_simulation_results.txt")

pdf("TrackSig_KL.pdf", width = 5, height=5)
plot(res$kl, res$abs_diff_max, main="TrackSig KL",
   xlab="KL", ylab="max abs diff")
dev.off()

