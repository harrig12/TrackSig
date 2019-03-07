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

kldiv_multinomials <- function(multinom1, multinom2) {
  return(apply(multinom1 * log(multinom1/multinom2),1,sum))
}

compare_simulation_results  <- function(simulation_list, 
  ground_truth_dir, method_results_dir, res_file_name) {

  results_df <- c()
  gt_exposures_list <- list()
  estim_exposures_list <- list()
  
  for (sim in simulation_list) {
    print(sprintf("Processing simulation %s ...", sim))

    gt_exposures_file = paste0(ground_truth_dir, "/", sim, "/", sim, "_sig_exp_per_mut.txt")

    if (!file.exists(gt_exposures_file)) {
      print(sprintf("File %s not found", gt_exposures_file))
      next
    }

    gt_exposures <- read.delim(gt_exposures_file, header=T, stringsAsFactors=F)
    gt_pos = paste0(gt_exposures[,"chromosome"], "_", gt_exposures[,"start"])

    estim_exposures_file = paste0(method_results_dir, "/", sim, "/", "sig_exposures_per_mut.txt")

    if (!file.exists(estim_exposures_file)) {
      print(sprintf("File %s not found", estim_exposures_file))
      next
    }

    estim_exposures <- read.delim(estim_exposures_file, header=T, stringsAsFactors=F)
    stopifnot(dim(estim_exposures) == dim(gt_exposures))

    # Making the rows be in the same order
    estim_pos = paste0(estim_exposures[,"chromosome"], "_", estim_exposures[,"start"])
    rownames(estim_exposures) <- estim_pos
    estim_exposures = estim_exposures[gt_pos,]

    # Making the signatures be in the same order
    stopifnot(sum(sort(colnames(estim_exposures)) == sort(colnames(gt_exposures))) == ncol(estim_exposures))
    estim_exposures = estim_exposures[,colnames(gt_exposures)]
    stopifnot(colnames(estim_exposures) == colnames(gt_exposures))

    gt_exposures_d = gt_exposures[,3:ncol(gt_exposures)]
    estim_exposures_d = estim_exposures[,3:ncol(estim_exposures)]

    abs_diff <- abs(gt_exposures_d - estim_exposures_d)

    estim_exposures_d_eps <- estim_exposures_d
    estim_exposures_d_eps[estim_exposures_d_eps == 0] = 0.01
    estim_exposures_d_eps = estim_exposures_d_eps / apply(estim_exposures_d_eps, 1, sum)

    d <- data.frame(sim = sim, 
      abs_diff_mean = mean(unlist(abs_diff), na.rm=TRUE),
      abs_diff_max = max(unlist(abs_diff), na.rm=TRUE),
      kl = mean(kldiv_multinomials(gt_exposures_d, estim_exposures_d_eps), na.rm=TRUE))

    gt_exposures_list[[sim]] <- gt_exposures
    estim_exposures_list[[sim]] <- estim_exposures

    results_df <- rbind(results_df, d)
  }

  results_df <- data.frame(results_df)
  rownames(results_df) <- results_df[,1]
  write.table(results_df, file = res_file_name, sep = "\t", row.names=F, quote=F)
  
  return(list(results_df, gt_exposures_list, estim_exposures_list))
}


plot_kl_results <- function(res, method, sim_type = "") {
  if (sim_type != "") {
    sim_type = paste0(".", sim_type)
  }
  pdf(paste0(method, "_simulation_results_KL_vs_max", sim_type, ".pdf"), width = 5, height=5)
  plot(res$kl, res$abs_diff_max, main=paste0(method, " versus Truth"), 
     xlab="KL", ylab="max abs diff", xlim=c(0, 1), ylim=c(0, 1))
  dev.off()

  pdf(paste0(method,"_simulation_results_KL_vs_mean", sim_type, ".pdf"), width = 5, height=5)
  plot(res$kl, res$abs_diff_mean, main=paste0(method, " versus Truth"), 
     xlab="KL", ylab="mean abs diff", xlim=c(0, 1), ylim=c(0, 1))
  dev.off()

}

get_sig_names_from_list <- function(exp_per_mut_list) {
  return(lapply(exp_per_mut_list, function(x) colnames(x[3:6])))
}


compare_changepoints  <- function(simulation_list, ground_truth_dir,
  tracksig_results_dir, sciclone_results_dir, res_file_name) {

  results_df <- c()
  
  for (sim in simulation_list) {
    print(sprintf("Processing simulation %s ...", sim))
    tracksig_cp_file = paste0(tracksig_results_dir, "/", sim, "/changepoints.txt")

    if (!file.exists(tracksig_cp_file)) {
      print(sprintf("File %s not found", tracksig_cp_file))
      next
    }

   if (file.info(tracksig_cp_file)$size == 1) {
      cp_tracksig <- 0
    } else {
      cp_tracksig <- ncol(read.table(tracksig_cp_file, header=F))
    }

    sciclone_cp_file = paste0(sciclone_results_dir, "/", sim, "/phis.txt")

    if (!file.exists(sciclone_cp_file)) {
      print(sprintf("File %s not found", sciclone_cp_file))
      next
    }
    # counting change-points, not the number of clusters
    cp_sciclone <- nrow(read.table(sciclone_cp_file, header=F, stringsAsFactors=F)) - 1
    
    gt_exposures_file = paste0(ground_truth_dir, "/", sim, "/", sim, "_exposures.txt")

    if (!file.exists(gt_exposures_file)) {
      print(sprintf("File %s not found", gt_exposures_file))
      next
    }

    gt_exposures <- read.csv(gt_exposures_file, header=T, stringsAsFactors=F)
    rownames(gt_exposures) <- gt_exposures[,1]
    gt_exposures <- gt_exposures[,-1]

    gt_delta_exp <- abs(gt_exposures[,2:ncol(gt_exposures)]-gt_exposures[,1:(ncol(gt_exposures)-1)])
    n_gt_exposure_cp <- sum(apply(gt_delta_exp,2,sum) > 0.1)

    if (grepl("two_clusters", sim)) {
      n_gt_created_cp = 1
    } else if (grepl("one_cluster", sim)) {
      n_gt_created_cp = 0
    } else if (grepl("inf_site_viol", sim)) {
      n_gt_created_cp = 1 # !!!!!!!!! change to 2 once I re-run simulations
    } else {
      n_gt_created_cp = 2
    }

    d <- data.frame(sim = sim, 
      cp_tracksig = cp_tracksig,
      cp_sciclone = cp_sciclone,
      n_gt_exposure_cp = n_gt_exposure_cp,
      n_gt_created_cp = n_gt_created_cp,
      stringsAsFactors = F)

    results_df <- rbind(results_df, d)
  }

  rownames(results_df) <- results_df[,1]
  write.table(results_df, file = res_file_name, sep = "\t", row.names=F, quote=F)
  
  return(results_df)
}


plot_per_sim_type <- function(res, method_name) {
  sim_types <- c("one_cluster", "two_clusters", 
    "branching", "cna_plus", 
    "inf_site_viol_plus")

  print("Per sim type")
  for (type in sim_types) {
    idx <- grepl(type, sapply(res[,1],toString))

    if (sum(idx) > 1) {
      plot_kl_results(res[idx,], "SciClone", sim_type = type)
    }
    print(paste0("Percentage of samples with KL larger than 0.05 in type ", type))
    print(mean(res[idx,]$kl > 0.05))
  }

  print("Per depth")
  depth_types <- c("depth30$", "depth100$", "depth1000$")
  for (type in depth_types) {
    idx <- grepl(type, sapply(res[,1],toString))

    if (sum(idx) > 1) {
      plot_kl_results(res[idx,], "SciClone", sim_type = type)
    }
    print(paste0("Percentage of samples with KL larger than 0.05 in type ", type))
    print(mean(res[idx,]$kl > 0.05))
  }
}
