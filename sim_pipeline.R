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

for (sim_i in 1:length(simulations)){
  # tracksig - make counts
  print(sprintf("%s", simulations[sim_i]))
  system(sprintf("src/run_simulations.sh data/mut_types/ data/%s %s", simulations[sim_i], simulations[sim_i]))
}


pcawg_format = T
tumortype_file <- sim_tumortype_file
purity_file <- sim_purity_file
active_signatures_file <- sim_activities_file
tracksig_results_dir = DIR_RESULTS = "TS_results_signature_trajectories/"
source("src/init_tracksig.R")

# source(sprintf("src/init_tracksig.R --tumortypes %s --purity %s  --signatures %s --trinucleotide %s --active %s",
#   sim_tumortype_file, sim_purity_file, signature_file, trinucleotide_file, sim_activities_file))

# tracksig - compute mutational signatures
compute_signatures_for_all_examples(countsDir = "data/counts/", bootstrapDir = "data/bootstrap/")

extract_exposures_per_mutation(activities_dir = paste0(tracksig_results_dir, "/SIMULATED/"),
  sorted_mutations_dir = "data/mut_types/", bin_size = 100)

# tracksig - bootstrap not functional yet
#compute_errorbars_for_all_examples()




