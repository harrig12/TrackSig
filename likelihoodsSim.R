# scoringSim.R
library(TrackSig)

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


TrackSig.options(purity_file = sim_purity_file,
                 signature_file = "annotation/sigProfiler_SBS_signatures.txt",
                 trinucleotide_file = "annotation/trinucleotide.txt",
                 active_signatures_file = sim_activities_file,
                 tumortype_file = sim_tumortype_file,
                 sig_amount = "onlyKnownSignatures",
                 compute_bootstrap = FALSE,
                 cancer_type_signatures = FALSE,
                 pcawg_format = TRUE)



# tracksig - make counts
for (sim_i in 1:length(simulations)){
  print(sprintf("%s", simulations[sim_i]))
  run_simulation(simulations[sim_i])
}

# show results for 3 likelihood functions; sig only, ccf only, sig+ccf

# options for 3 functions
DIR_RESULTS_list <- c("sig", "ccf", "sigccf")

pelt_penalty_list <- c(expression((n_sigs - 1) * log(n_bins)), expression( 2 * log(n_bins)),
                      expression((n_sigs - 1) * log(n_bins) + 2 * log(n_bins)))

pelt_score_fxn_list <- c(TrackSig:::log_likelihood_mixture_multinomials, TrackSig:::gaussian_ll,
                         TrackSig:::sum_gaussian_mixture_multinomials_ll)


for (likelihood_i in 1:3){

  #set options for this scoring
  TrackSig.options(DIR_RESULTS = paste0(DIR_RESULTS_list[likelihood_i], "/"),
                   pelt_penalty = pelt_penalty_list[likelihood_i],
                   pelt_score_fxn = pelt_score_fxn_list[[likelihood_i]])


  # tracksig - compute mutational signatures with settings
  #compute_signatures_for_all_examples(countsDir = "data/counts", bootstrapDir = "data/bootstrap/")

  # tracksig - get exposures
  extract_exposures_per_mutation(activities_dir = paste0(DIR_RESULTS_list[likelihood_i], "/SIMULATED/"),
                                 sorted_mutations_dir = "data/mut_types/")

  # create plots and compare simulation with ground truth
  capture.output(TrackSig:::compare_simulations(DIR_RESULTS_list[likelihood_i],
                                                dataDir = "data",
                                                outDir = DIR_RESULTS_list[likelihood_i]),
                 file = sprintf("%s/compare_summary_%s.txt",
                                DIR_RESULTS_list[likelihood_i],DIR_RESULTS_list[likelihood_i]))
}

# [END]
