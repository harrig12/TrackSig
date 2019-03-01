# sciClone_deconstructSigs_pipeline.R

#################
# set up
#################

library(sciClone)
library(tidyr)
library(deconstructSigs)
library(TrackSig)

#setwd("~/Desktop/TrackSig_simulations/")
getwd()

# get simulation names
simNames <- list.files("data")
sel <- grep(x = simNames, "^Simulation")
simNames <- simNames[sel]

# annotation files
sim_activities_file <- "annotation/sim_active_in_sample.txt"
signature_file <- "annotation/sigProfiler_SBS_signatures.txt"


# to be parallelized
for (simName in simNames){ #for each simulation

  print(simName)

  # outdir
  resultsDir <- paste0("SCDS_results/SIMULATED", "/", simName)
  dir.create(resultsDir, showWarnings=FALSE, recursive = TRUE)

  #################
  # sciclone
  #################

  # read in vaf data
  vafTable <- read.table(sprintf("data/%s/%s_scicloneDS.txt", simName, simName),
                         stringsAsFactors = F, header = T)

  # format for sciClone input
  vafTable <- separate(vafTable, chr, into = c("V1", "chr"), sep = "chr")
  vafTable$V1 <- NULL

  # scale vaf by 100
  vafTable$vaf <- apply(vafTable[5], MARGIN = 1, prod, 100)

  # read in CNA data
  cnaTable <- read.table(sprintf("data/%s/%s_cna.txt", simName, simName), header = T)

  # sciclone - compute clusters
  sC <- sciClone::sciClone(vafTable[, (1:6)], copyNumberCalls = cnaTable, sampleNames = simName,
                           minimumDepth = 20, useSexChrs = FALSE)

  # get cluster assignments
  vafTable$cluster <- sC@clust$cluster.assignments

  #################
  # deconstructSigs
  #################

  # select by sample name and activity
  activeSigs <- read.delim(sim_activities_file, header = T, stringsAsFactors = F)
  activeSigs <- activeSigs[activeSigs$Sample_Name == simName,]
  sel <- colnames(activeSigs[,3:ncol(activeSigs)])[colSums(activeSigs[,3:ncol(activeSigs)]) > 0]

  # get signatures active for selection
  #allSigs <- read.delim(signature_file, header = T, stringsAsFactors = F)
  #allSigs <- setNames(data.frame(t(allSigs[,-1])), allSigs[,1])
  allSigs <- TrackSig:::load_sim_signatures(signature_file)
  allSigs <- as.data.frame(t(allSigs))
  activeSigs <- allSigs[sel,]

  # collect colnames from active signatures for mutation types
  allMutTypes <- data.frame(rep(0, length(colnames(activeSigs))), row.names = colnames(activeSigs))

  # get mutation type counts for simulation
  vafTable$mutType <- apply(vafTable[, c("ref", "alt", "tri")], MARGIN = 1, paste, collapse = "_")

  exposurePerCluster <- NULL
  exposurePerMut <- vafTable[,c("chr","pos","cluster")]
  exposurePerMut[,row.names(activeSigs)] <- NA

  # subset on custer and deconstructSigs
  for (cluster_i in unique(vafTable$cluster)){
    if (cluster_i == 0) {
      next
    }
    # subset mutations for cluster
    clusterVafTable <- subset(vafTable, cluster == cluster_i)

    # tabulate present mutation types
    activeMutTypes <- as.data.frame(table(clusterVafTable$mutType), row.names = 1)
    typesTable <- allMutTypes
    typesTable[row.names(activeMutTypes), ] <- activeMutTypes

    colnames(typesTable) <- sprintf("%s_cluster%s", simName, cluster_i)

    # transpose for deconstructSigs input
    typesTable <- as.data.frame(t(typesTable))

    # find signatures with deconstructSigs
    foundSigs <- whichSignatures(tumor.ref = typesTable, contexts.needed = T, signatures.ref = activeSigs)
    #write.table(foundSigs$weights, sprintf("%s/signatures_%s", resultsDir, simName), append = T, quote = FALSE, row.names = FALSE)

    exposures_i <- foundSigs$weights
    rownames(exposures_i) <- cluster_i

    exposurePerCluster <- rbind(exposurePerCluster, exposures_i)
    exposurePerMut[exposurePerMut$cluster == cluster_i, colnames(exposures_i)] <- exposures_i

  }

  #################
  # output
  #################

  stopifnot(length(sC@clust$cluster.means) == nrow(exposurePerCluster))
  
  # phis
  phis <- sC@clust$cluster.means
  write.table(t(phis), file = sprintf("%s/%s", resultsDir, "phis.txt"), quote = F, row.names = F, col.names = F)

  # mixtures
  mixtures <- exposurePerCluster[ order(row.names(exposurePerCluster)), ]
  mixtures <- as.data.frame(t(mixtures))
  colnames(mixtures) <- phis

  write.csv(mixtures, file = sprintf("%s/%s", resultsDir, "mixtures.csv"), quote = T, row.names = T, col.names = T)

  # sig_exposures_per_mut
  exposurePerMut$cluster <- NULL
  write.table(exposurePerMut, file = sprintf("%s/%s", resultsDir, "sig_exposures_per_mut.txt"), quote = F, row.names = F, col.names = T)

  sc.plot1d(sC, sprintf("%s/%s", resultsDir, "sciclone.pdf"))

  # plot
  plotName <- sprintf("%s/%s_%s", resultsDir, simName, "trajectory.pdf")
  TrackSig:::plot_signatures(mixtures*100, plot_name = plotName, phis = phis, mark_change_points = F,
                  change_points = NULL, transition_points = NULL, scale=1.2, save = T)

}


# [END]

