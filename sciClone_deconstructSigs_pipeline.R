# sciClone_deconstructSigs_pipeline.R

# sciClone + DeconstructSigs pipeline

library(sciClone)
library(deconstructSigs)

setwd("~/Desktop/TrackSig_simulations/")

# get simulation names
simNames <- list.files("data/mut_types/")
simNames <- strsplit(simNames, ".mut_types.txt")

for (simName in simNames){ #for each simulation

  # read in vaf data
  vafTable <- read.table(sprintf("data/mut_types/%s.mut_types.txt", simName))

  # format for sciClone input
  vafTable <- vafTable[,c(1,2,4,5,3)]
  colnames(vafTable) <- c("chr", "pos", "ref_reads", "var_reads", "vaf")

  # simulate read counts

  # scale vaf by 100
  vafTable$vaf <- apply(vafTable[5], MARGIN = 1, prod, 100)

  # read in CNA data
  cnaTable <- read.table(sprintf("data/%s/%s_cna.txt", simName, simName), header = T)

  sC <- sciClone::sciClone(vafTable, copyNumberCalls = cnaTable, sampleNames = simName)

}



# [END]
