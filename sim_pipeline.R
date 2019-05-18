#################
# TrackSig Simulations
# Author: Cait Harrigan
#################
library(TrackSig)

# set up
tracksig_results_dir = "TS_results_signature_trajectories/"
outdir <- "data"
signature_file <- "annotation/sigProfiler_SBS_signatures.txt"

mut_per_sim <- 5000

# create data directory
dir.create(outdir, showWarnings=FALSE)

# generate simulations
####################################################

# populate annotation files
sim_activities_file <- "annotation/sim_active_in_sample.txt"
sim_purity_file <- "annotation/sim_purity.txt"
sim_tumortype_file <-"annotation/sim_tumortypes.txt"


# Create simulation vcfs
print(sprintf("Generating simulations ..."))
simulations = create_simulation_set(outdir = outdir,
                                   mut_per_sim = mut_per_sim,
                                   sim_activity_file = sim_activities_file,
                                   sim_purity_file = sim_purity_file,
                                   sim_tumortype_file = sim_tumortype_file,
                                   signature_file = signature_file,
                                   rewrite_annotations = T)

# Run TrackSig
####################################################

# grab simulation names
simulations <- list.files(outdir)
sel <- grep(x = simulations, "^Simulation")
simulations <- simulations[sel]

# TrackSig - set options (previously header.R)
TrackSig.options(purity_file = sim_purity_file,
                 signature_file = "annotation/sigProfiler_SBS_signatures.txt",
                 trinucleotide_file = "annotation/trinucleotide.txt",
                 active_signatures_file = sim_activities_file,
                 tumortype_file = sim_tumortype_file,
                 sig_amount = "onlyKnownSignatures",
                 compute_bootstrap = FALSE,
                 cancer_type_signatures = FALSE,
                 pcawg_format = TRUE,
                 DIR_RESULTS = tracksig_results_dir)

# TrackSig - make counts
for (sim_i in 1:length(simulations)){
  print(sprintf("%s", simulations[sim_i]))
  run_simulation(simulations[sim_i], outdir)
}

# TrackSig - compute mutational signatures
compute_signatures_for_all_examples(countsDir = "data/counts", bootstrapDir = "data/bootstrap/",
  samples_to_run = simulations)

# TrackSig - get exposures
extract_exposures_per_mutation(activities_dir = paste0(tracksig_results_dir, "/SIMULATED/"),
                               sorted_mutations_dir = "data/mut_types/")


# Run simulations for different bin sizes
####################################################

# set up
mut_per_sim <- 6000
outdir <- "data_bin_simulations/"
dir.create(outdir, showWarnings=FALSE)

simulations = create_simulation_bin_sizes(outdir = outdir,
                                    mut_per_sim = mut_per_sim,
                                    sim_activity_file = sim_activities_file,
                                    sim_purity_file = sim_purity_file,
                                    sim_tumortype_file = sim_tumortype_file,
                                    signature_file = signature_file)


# simulations by bin size
bin_sizes <- c(25, 50, 75, 100, 150, 200, 300, 500)

simulations <- list.files(outdir)
sel <- grep(x = simulations, "^Simulation")
simulations <- simulations[sel]

for (bin_size in bin_sizes){
  sim_bin_size_idx = as.integer(gsub(".*_bin(.*)_depth.*", "\\1", simulations)) == bin_size
  simulations_w_bin_size = simulations[sim_bin_size_idx]
  simulations_w_bin_size <- sample(simulations_w_bin_size, length(simulations_w_bin_size))

  # TrackSig - set current bin size
  TrackSig.options(bin_size = bin_size)

  # TrackSig - make counts
  for (sim_i in 1:length(simulations_w_bin_size)){
    print(sprintf("%s", simulations_w_bin_size[sim_i]))
    run_simulation(simulations_w_bin_size[sim_i], outdir)
  }

  print("Computing signatures....")
  for (sim in simulations_w_bin_size) {
    tryCatch({
      # TrackSig - compute mutational signatures
      compute_signatures_for_all_examples(countsDir = paste0(outdir,"/counts"),
                                        bootstrapDir = paste0(outdir,"/bootstrap/"),
                                        samples_to_run = c(sim))

       # TrackSig - get exposures
      extract_exposures_per_mutation(activities_dir = paste0(tracksig_results_dir, "/SIMULATED/"),
                                   sorted_mutations_dir = paste0(outdir,"/mut_types/"),
                                   samples_to_run = c(sim), bin_size = bin_size)

      }, warning = function(war) {
      }, error = function(err) {

        print(paste0("Re-running simulation ", sim))

        unlink(paste0(outdir,"/counts/", sim, ".phi.txt"))
        unlink(paste0(outdir,"/counts/", sim, ".quadraticp.txt"))
        unlink(paste0(outdir,"/mut_types/", sim, ".mut_types.txt"))
        unlink(paste0(tracksig_results_dir, "/SIMULATED/", sim), recursive = TRUE)

        run_simulation(sim, outdir)

        tryCatch({
        # TrackSig - compute mutational signatures
        compute_signatures_for_all_examples(countsDir = paste0(outdir,"/counts"),
                                        bootstrapDir = paste0(outdir,"/bootstrap/"),
                                        samples_to_run = c(sim))

           # TrackSig - get exposures
         extract_exposures_per_mutation(activities_dir = paste0(tracksig_results_dir, "/SIMULATED/"),
                                     sorted_mutations_dir = paste0(outdir,"/mut_types/"),
                                     samples_to_run = c(sim), bin_size = bin_size)
          }, warning = function(war) {
          }, error = function(err) {
              # if we failed the second time with this simulation, just skip it
              next
          }, finally = {
          })

      }, finally = {
      })
  }
}
