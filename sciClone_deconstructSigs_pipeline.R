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

# shuffle
simNames <- sample(simNames, length(simNames))
#simNames <- simNames[(i*60):((i+1)*60)]

# binSize parameter
binSize <- 100

times <- c()
# to be parallelized
for (simName in simNames){ #for each simulation

  print(sprintf("Processing %s: %s/%s", 
    simName, which(simNames == simName), length(simNames)))

  #tic()

  # outdir
  resultsDir <- paste0("SCDS_results/SIMULATED", "/", simName)
  dir.create(resultsDir, showWarnings=FALSE, recursive = TRUE)

  if (file.exists(sprintf("%s/%s", resultsDir, "mixtures.csv"))) {
    print(sprintf("Skipping %s...", simName))
    next
  }
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
  # If using non-zero minimumDepth, sC@clust$cluster.assignments is filtered by depth, 
  # so it has less elements than  vafTable$cluster
  # I don't know what to do with this -- I didn't find the list of remaining mutations in sC object
  skip_sim = FALSE
  tryCatch({
     sC <- sciClone::sciClone(vafTable[, (1:6)], copyNumberCalls = cnaTable, sampleNames = simName,
                           minimumDepth = -1, useSexChrs = FALSE, clusterMethod = 'binomial.bmm')

    }, warning = function(war) {
    }, error = function(err) {
      print("Error in sciClone:")
      print(err)
      skip_sim = TRUE
    }, finally = {
   }) 

  if (skip_sim) {
    next
  }

 
 if (length(sC@clust$cluster.assignments) == nrow(vafTable)) {
    # get cluster assignments
    vafTable$cluster <- sC@clust$cluster.assignments
  } else {
     # Remove mutations with depth 0
    vafTable <- vafTable[vafTable[,3] + vafTable[,4] >= 1,]
    
    # get cluster assignments
    vafTable$cluster <- sC@clust$cluster.assignments
  }


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

  #check signature names line up
  if (all(sel %in% rownames(allSigs)) == FALSE){
    stop("active signatures and reference signature names inconsistant")
  }

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

    mut_cluster_idx <- exposurePerMut$cluster == cluster_i
    exposurePerMut[mut_cluster_idx, colnames(exposures_i)] <- exposures_i[rep(1,sum(mut_cluster_idx)),]

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

  # repeat columns proportional to mutations per cluster
  mutPerClust <- c(table(vafTable$cluster)) %/% binSize

  write.csv(mixtures, file = sprintf("%s/%s", resultsDir, "mixtures.csv"), quote = T, row.names = T)

  # sig_exposures_per_mut
  exposurePerMut$cluster <- NULL
  exposurePerMut[,"chr"] <- paste0("chr", exposurePerMut[,"chr"])
  colnames(exposurePerMut)[c(1,2)] <- c("chromosome" ,"start")
  write.table(exposurePerMut, file = sprintf("%s/%s", resultsDir, "sig_exposures_per_mut.txt"), sep = "\t", row.names=F, quote=F)

  sc.plot1d(sC, sprintf("%s/%s_%s", resultsDir, simName, "sciclone.pdf"))

  # plot
  plotName <- sprintf("%s/%s_%s", resultsDir, simName, "trajectory.pdf")
  
  if (ncol(mixtures) > 1) {
    TrackSig:::plot_signatures(mixtures*100, plot_name = plotName, phis = phis * 2, mark_change_points = F,
                    change_points = NULL, transition_points = NULL, scale=1.2, save = T)
  }

  #toc(log=T)
  #beep(2)
}


# [END]

