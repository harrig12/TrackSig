# sciClone_deconstructSigs_pipeline.R

#################
# set up
#################

library(sciClone)
library(tidyr)
library(deconstructSigs)

#setwd("~/Desktop/TrackSig_simulations/")
getwd()

# get simulation names
simNames <- list.files("data")
sel <- grep(x = simNames, "^Simulation")
simNames <- simNames[sel]

# annotation files
sim_activities_file <- "annotation/sim_active_in_sample.txt"
signature_file <- "annotation/sigProfiler_SBS_signatures.txt"

# outdirs
scicloneDir <- "sciclone_results"
deconstructSigsDir <- "deconstructSigs_results"

dir.create(scicloneDir, showWarnings=FALSE)
dir.create(deconstructSigsDir, showWarnings=FALSE)

# to be parallelized
for (simName in simNames){ #for each simulation

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
  sC <- sciClone::sciClone(vafTable[,(1:6)], copyNumberCalls = cnaTable, sampleNames = simName,
                           minimumDepth = 20, useSexChrs = FALSE)

  writeClusterSummaryTable(sC, sprintf("%s/summary_%s", scicloneDir, simName))
  writeClusterTable(sC, sprintf("%s/clusters_%s", scicloneDir, simName))

  # get cluster assignments
  vafTable$cluster <- sC@clust$cluster.assignments


  #################
  # deconstructSigs
  #################

  # select by sample name and activity
  activeSigs <- read.delim(sim_activities_file, header = T, stringsAsFactors = F)
  activeSigs <- activeSigs[activeSigs$Sample_Name == simName,]
  sel <- colSums(activeSigs[,3:ncol(activeSigs)]) > 0

  # get signatures active for selection
  allSigs <- read.delim(signature_file, header = T, stringsAsFactors = F)
  allSigs <- setNames(data.frame(t(allSigs[,-1])), allSigs[,1])
  activeSigs <- allSigs[sel,]

  # get mutation type counts for simulation
  vafTable$mut_type <- apply(vafTable[, c("tri", "alt")], MARGIN = 1, paste, collapse = ">")

  # subset on custer and deconstructSigs
  for (cluster_i in unique(vafTable$cluster)){
    clusterVafTable <- subset(vafTable, cluster == cluster_i)

    # subset mutations for cluster
    typesTable <- as.data.frame(t(as.matrix(table(clusterVafTable$mut_type))))
    typesTable[1,] <- apply(typesTable[1,], MARGIN = 1, as.numeric)
    rownames(typesTable) <- sprintf("%s_cluster%s", simName, cluster_i)

    # some mutation types could be missing from simulation randomness and subsetting
    simTypes <- colnames(typesTable)
    activeTypes <- colnames(activeSigs)
    if (length(simTypes) < length(activeTypes)){

      # find which missing
      nDiff <- length(activeTypes) - length(simTypes)

      # add them into typesTable with frequency 0
      missingTypes <- activeTypes[!(activeTypes %in% simTypes)]

      for (missingType in missingTypes){ #can apply work here?
        typesTable[[missingType]] <- 0
      }
    }

    foundSigs <- whichSignatures(tumor.ref = typesTable, contexts.needed = T, signatures.ref = activeSigs)
    write.table(foundSigs$weights, sprintf("%s/signatures_%s", deconstructSigsDir, simName),
                append = T, quote = FALSE, row.names = FALSE)
  }


}


# [END]

