source("sim_pipeline_helpers.R")

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

outdir <- "data"

# ===========================================
# Compare TrackSig exposures to ground truth
simulations <- list.files(outdir)
sel <- grep(x = simulations, "^Simulation")
simulations <- simulations[sel]

tracksig_results_dir = "TS_results_signature_trajectories/"

list[res, gt_exp_l, estim_exp_l] <- compare_simulation_results(
    simulations, 
    ground_truth_dir = outdir, 
    method_results_dir = paste0(tracksig_results_dir, "/SIMULATED/"),
    res_file_name = "TrackSig_simulation_results.txt")

res_TrackSig <- res
plot_kl_results(res, "TrackSig")

print("Percentage of samples with KL larger than 0.05")
mean(res$kl > 0.05)

plot_per_sim_type(res, "TrackSig")




simulations_large_kl <- res[res$kl > 0.05,"sim"]
sigs_in_bad_sims <- get_sig_names_from_list(gt_exp_l[simulations_large_kl])
table(unlist(sigs_in_bad_sims))
# Signature 7 is present in ALL the simulations where TrackSig results differ from ground truth

# has_sig7 <- list()
# for (sim in names(gt_exp_l)) {
#     has_sig7[[sim]] <- "SBS7" %in% colnames(gt_exp_l[[sim]][3:6])
# }
# sim_order <- res_TrackSig[,1]
# has_sig7 = has_sig7[sim_order]
# has_sig7 <- as.numeric(data.frame(has_sig7))
# has_sig7[has_sig7 == 0] <- "darkgreen"
# has_sig7[has_sig7 == 1] <- "magenta"

pdf("TrackSig_simulation_results_KL_vs_mean.pdf", width = 5, height=5)
plot(res_TrackSig$kl, res_TrackSig$abs_diff_mean, main="TrackSig versus Truth", 
   xlab="KL", ylab="mean abs diff",  xlim=c(0, 1), ylim=c(0, 1)) #col=has_sig7)
#legend('topleft', legend = c("No SBS7", "Has SBS7"), col = c("darkgreen", "magenta"), cex = 0.8, pch = 1)
dev.off()

pdf("TrackSig_simulation_results_KL_vs_max.pdf", width = 5, height=5)
plot(res_TrackSig$kl, res_TrackSig$abs_diff_max, main="TrackSig versus Truth", 
   xlab="KL", ylab="max abs diff", xlim=c(0, 1), ylim=c(0, 1)) # col=has_sig7)
#legend('topleft', legend = c("No SBS7", "Has SBS7"), col = c("darkgreen", "magenta"), cex = 0.8, pch = 1)
dev.off()




# ===========================================
# Compare SciClone exposures to ground truth
simulations <- list.files(outdir)
sel <- grep(x = simulations, "^Simulation")
simulations <- simulations[sel]

sciclone_dir <- "SCDS_results/"
list[res, gt_exp_l, estim_exp_l] <- compare_simulation_results(simulations, 
    ground_truth_dir = outdir, 
    method_results_dir = paste0(sciclone_dir, "/SIMULATED/"),
    res_file_name = "sciclone_simulation_results.txt")

res_SciClone <- res
estim_exp_l_sciclone <- estim_exp_l
plot_kl_results(res, "SciClone")

print("sciclone: Percentage of samples with KL larger than 0.05")
mean(res$kl > 0.05)

plot_per_sim_type(res, "SciClone")




# has_sig7 <- list()
# for (sim in names(gt_exp_l)) {
#     has_sig7[[sim]] <- "SBS7" %in% colnames(gt_exp_l[[sim]][3:6])
# }

# sim_order <- res_SciClone[,1]
# has_sig7 = has_sig7[sim_order]

# has_sig7 <- as.numeric(data.frame(has_sig7))
# has_sig7[has_sig7 == 0] <- "darkgreen"
# has_sig7[has_sig7 == 1] <- "magenta"

pdf("SciClone_simulation_results_KL_vs_mean.pdf", width = 5, height=5)
plot(res_SciClone$kl, res_SciClone$abs_diff_mean, main="SciClone versus Truth", 
   xlab="KL", ylab="mean abs diff", xlim=c(0, 1), ylim=c(0, 1)) #col=has_sig7)
#legend('topleft', legend = c("No SBS7", "Has SBS7"), col = c("darkgreen", "magenta"), cex = 0.8, pch = 1)
dev.off()

pdf("SciClone_simulation_results_KL_vs_max.pdf", width = 5, height=5)
plot(res_SciClone$kl, res_SciClone$abs_diff_max, main="SciClone versus Truth", 
   xlab="KL", ylab="max abs diff", xlim=c(0, 1), ylim=c(0, 1)) #col=has_sig7)
#legend('topleft', legend = c("No SBS7", "Has SBS7"), col = c("darkgreen", "magenta"), cex = 0.8, pch = 1)
dev.off()






# ===========================================
# Compare TrackSig and SciClone
sim_order <- intersect(res_SciClone[,1], res_TrackSig[,1])
print("TrackSig has lower KL in this % of simulations:")
mean(res_TrackSig[sim_order,]$kl < res_SciClone[sim_order,]$kl)

print("Mean KL between the method and ground truth")
mean(res_TrackSig$kl)
mean(res_SciClone$kl)

print("mean and median KL diff between two methods")
mean(res_SciClone[sim_order,]$kl - res_TrackSig[sim_order,]$kl)
median(res_SciClone[sim_order,]$kl - res_TrackSig[sim_order,]$kl)

print("TrackSig has lower mean abs diff in this % of simulations:")
mean(res_TrackSig[sim_order,]$abs_diff_mean < res_SciClone[sim_order,]$abs_diff_mean)

print("TrackSig has lower max abs diff in this % of simulations:")
mean(res_TrackSig[sim_order,]$abs_diff_max < res_SciClone[sim_order,]$abs_diff_max)

print("percentage of correct reconstructions")
mean(res_TrackSig$kl < 0.05)
mean(res_SciClone$kl < 0.05)

mean(res_TrackSig$kl < 0.1)
mean(res_SciClone$kl < 0.1)


correct_tracksig = res_TrackSig[sim_order,]$kl < 0.05

correct_tracksig <- as.numeric(correct_tracksig)
correct_tracksig[correct_tracksig == 0] <- "red"
correct_tracksig[correct_tracksig == 1] <- "darkgreen"

pdf("sciclone_results_coloured_TrackSig_mean_diff.pdf", width = 5, height=5)
plot(res_SciClone[sim_order,]$kl, res_SciClone[sim_order,]$abs_diff_mean, 
  main="Sciclone versus Truth", 
   xlab="KL", ylab="mean abs diff", col=correct_tracksig,
   xlim=c(0, 1), ylim=c(0, 1))
legend('topleft', legend = c("TrackSig: KL > 0.05", "TrackSig: KL < 0.05"), col = c("red", "darkgreen"), cex = 0.8, pch = 1)
dev.off()

pdf("sciclone_results_coloured_TrackSig_max_diff.pdf", width = 5, height=5)
plot(res_SciClone[sim_order,]$kl, res_SciClone[sim_order,]$abs_diff_max, 
  main="Sciclone versus Truth", 
   xlab="KL", ylab="max abs diff", col=correct_tracksig,
   xlim=c(0, 1), ylim=c(0, 1))
legend('topleft', legend = c("TrackSig: KL > 0.05", "TrackSig: KL < 0.05"), col = c("red", "darkgreen"), cex = 0.8, pch = 1)
dev.off()




print("Samples where Sciclone makes mistakes but TrackSig doesn't")
sciclone_mistakes <- (res_SciClone[sim_order,]$kl > 0.1) & (res_TrackSig[sim_order,]$kl < 0.1)
sciclone_mistakes <- res_SciClone[sciclone_mistakes,1]

table(unlist(get_sig_names_from_list(estim_exp_l_sciclone[sciclone_mistakes])))
# More 1/3 have signature 2+13, but that does say much


# ===========================================
# Compare number of change-points
# simulations_depth100 <- simulations[grepl("depth100$", simulations)]
# simulations_depth1000 <- simulations[grepl("depth1000$", simulations)]
# Compare change-points
cp_comparison <- compare_changepoints(simulations, 
  ground_truth_dir = outdir,
  tracksig_results_dir = paste0(tracksig_results_dir, "/SIMULATED/"), 
  sciclone_results_dir = paste0(sciclone_dir, "/SIMULATED/"), 
  res_file_name = "cp_comparison.txt") 

{
print("TrackSig agrees with sciclone")
print(mean(cp_comparison$cp_tracksig == cp_comparison$cp_sciclone))

print("TrackSig agrees with GT")
print(mean(cp_comparison$n_gt_created_cp == cp_comparison$cp_tracksig))
print("SciClone agrees with GT")
print(mean(cp_comparison$n_gt_created_cp == cp_comparison$cp_sciclone))

print("TrackSig overestimates # change-points")
print(mean(cp_comparison$n_gt_created_cp < cp_comparison$cp_tracksig))
print("TrackSig underestimates # change-points")
print(mean(cp_comparison$n_gt_created_cp > cp_comparison$cp_tracksig))

print("SciClone overestimates # change-points")
print(mean(cp_comparison$n_gt_created_cp < cp_comparison$cp_sciclone))
print("SciClone underestimates # change-points")
print(mean(cp_comparison$n_gt_created_cp > cp_comparison$cp_sciclone))

print("Mean difference between ground truth and method")
print(mean(cp_comparison$n_gt_created_cp - cp_comparison$cp_tracksig))
print(mean(cp_comparison$n_gt_created_cp - cp_comparison$cp_sciclone))

print("Mean abs difference between ground truth and method")
print(mean(abs(cp_comparison$n_gt_created_cp - cp_comparison$cp_tracksig)))
print(mean(abs(cp_comparison$n_gt_created_cp - cp_comparison$cp_sciclone)))


#print("Examples where TrackSig makes mistakes but SciClone doesn't")
#print(cp_comparison[(cp_comparison$n_gt_created_cp != cp_comparison$cp_tracksig) & (cp_comparison$n_gt_exposure_cp == cp_comparison$cp_sciclone),])

#print(cp_comparison[(cp_comparison$n_gt_created_cp > cp_comparison$cp_tracksig) & (cp_comparison$n_gt_exposure_cp == cp_comparison$cp_sciclone),])
}

print("Comparison of number of CP per sim type")
sim_types <- c("one_cluster", "two_clusters", 
    "branching", "cna_plus", 
    "inf_site_viol_plus")

depth_types <- c("depth30", "depth100", "depth1000")
res_table <- data.frame(matrix(0, ncol = length(sim_types), nrow = length(depth_types)))
colnames(res_table) <- sim_types
rownames(res_table) <- depth_types

TrackSig_cp_summary <- SciClone_cp_summary <- res_table

for (sim_type in sim_types) {
  for (d_type in depth_types) {
    idx <- grepl(sim_type, sapply(cp_comparison[,1],toString)) 
    idx <- idx &  grepl(paste0(d_type,"$"), sapply(cp_comparison[,1],toString))

    # print(sim_type)
    # print(d_type)
    # print(mean(cp_comparison[idx,]$cp_tracksig ))
    # print(mean(cp_comparison[idx,]$n_gt_created_cp ))
    # print("TrackSig agrees with sciclone")
    # print(mean(cp_comparison[idx,]$cp_tracksig == cp_comparison[idx,]$cp_sciclone))

    # print("TrackSig agrees with GT")
    # print(mean(cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,]$cp_tracksig))
    # print("SciClone agrees with GT")
    # print(mean(cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,]$cp_sciclone))
    
    TrackSig_cp_summary[d_type, sim_type] <- mean((cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,]$cp_tracksig) | (cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,]$cp_tracksig +1))
    SciClone_cp_summary[d_type, sim_type] <- mean((cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,]$cp_sciclone) | (cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,]$cp_sciclone -1))
  }
}



