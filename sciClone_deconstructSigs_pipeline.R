# sciClone_deconstructSigs_pipeline.R

# sciClone + DeconstructSigs pipeline

library(sciClone)
library(tidyr)
library(deconstructSigs)

setwd("~/Desktop/TrackSig_simulations/")

# get simulation names
simNames <- list.files("data")
sel <- grep(x = simNames, "^Simulation")
simNames <- simNames[sel]

for (simName in simNames){ #for each simulation

  # read in vaf data
  vafTable <- read.table(sprintf("data/%s/%s_sciclone.txt", simName, simName),
                         stringsAsFactors = F, header = T)

  # format for sciClone input
  vafTable <- separate(vafTable, chr, into = c("V1", "chr"), sep = "chr")
  vafTable <- vafTable[,c(2:6)]

  # read in CNA data
  cnaTable <- read.table(sprintf("data/%s/%s_cna.txt", simName, simName), header = T)

  sC <- sciClone::sciClone(vafTable, copyNumberCalls = cnaTable, sampleNames = simName)

}



# [END]
