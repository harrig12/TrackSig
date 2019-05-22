# Simulations Branch 
For evaluation of TrackSig (Trackature) - A framework to infer mutational signatures over time.

## Background
Cell processes leave a unique signature of mutation types in cancer genome. Using the mutational signatures, it is possible to infer the fraction of mutations contributed by each cell process. Mutational signatures are represented as multinomial distributions over 96 mutation types. Using our framework Trackature, we can infer mutational signatures changing over time.

## See also the TrackSig [R package repo](https://github.com/harrig12/BCB430)

## Usage 
A brief overview page for generating simulated data and testing the performance of TrackSig, as well as comparing it to a pipeline of other packages (SciClone and deconstructSigs) is made available [here](https://harrig12.github.io/TrackSig/). This page assumes some familiarity with the TrackSig packages and design principles. To familiarize yourself in depth, the TrackSig paper (Rubanova et. al) is [available on bioRxiv](https://www.biorxiv.org/content/10.1101/260471v3)

## Important notes
1) **Is not applicable for samples with <600 mutations.** Please note that Trackature does not run on samples with less than 600 mutations. Less than 600 mutations will result in less than 3 time points, and there is no point to analize it as a time series. On tumors with less than 600 mutations, you can compute signature activities without dividing mutations into time points (see "Computing overall signature activities" section).

2) **Results are not re-computed when script is re-started.** Please note that at every step if you stop the script and re-start it, the computations will *continue* instead of re-writing the previous results. It is useful for launching large batches of samples: scripts can be paused when needed; if one sample fails, other samples don't need to be re-computed again. However, if you wish some results to be re-computed, please erase the corresponding directory.

3) **If the tumour names contain a dot, please replace it with another symbol, for example, with underscore.**


