## Phylodynamic simulations with reporting delays

# Load libraries
library(phylodyn2)
library(adapref)
library(gridExtra)
library(INLA)
library(ape)
library(phylotools)
library(here)
library(tidyverse)

#windowsFonts(A = windowsFont("Helvetica")) 

# generate-rd-tree-simulations -------------------------------------------------
sim_trees_rd_comparisons <- function(
    eff_pop_traj, # Function defining true effective population trajectory
    rd_fun,
    max_time, # Numeric specifying max sampling time
    time_grid_length, # Numeric specifying length of time grid
    pref_sample_c, # Numeric c in pref_sample()
    pref_sample_beta, # Numeric beta in pref_sample()
    truncation_time # Vector of times to truncated the data at
) {
  
  
  ## Get simulated data/trees
  # Simulate sampling events with preferential sampling (PS)
  simed_sampling_times <- c(
    0, 
    pref_sample(
      f = eff_pop_traj, 
      lim = c(0, max_time),
      grid.len = time_grid_length - 1,
      c = pref_sample_c, 
      beta = pref_sample_beta
    )
  )
  
  # Generate coalescence
  simed_full_coal <- suppressWarnings(coalsim(
    samp_times = simed_sampling_times, 
    n_sampled = rep(1, length(simed_sampling_times)), 
    traj = eff_pop_traj,
  ))
  
  # Simulate full tree
  full_tree <- suppressWarnings(generate_newick(simed_full_coal))
  
  # Save trees together
  trees <- vector(mode = "list", length = 3)
  names(trees) <- c("full", "obs", paste0("trunc", truncation_time))
  
  # Save full tree
  trees$full <- full_tree$newick
  
  # Match coalescence to sampling times
  tips_matched <- data.frame(
    connect_sample_to_tips(full_tree$newick)$ungrouped_samp_df,
    "true0_samp_times" = simed_sampling_times
  )
  
  # Simulate tips to be dropped due to Reporting Delay (RD)
  sim_reporting_delay <- function(samp_times, prob_reported_func) {
    sim_data <- data.frame("samp_times" = samp_times) 
    sim_data$prob_reported <- prob_reported_func(samp_times) 
    sim_data$reported <- mapply(
        FUN = rbinom, 
        n = 1, 
        p = sim_data$prob_reported,  
        MoreArgs = list(size = 1)
      )
    
    ifelse(sim_data$reported == 0, yes = TRUE, no = FALSE)
  }
  
  tips_matched$drop <- sim_reporting_delay(
    samp_times = tips_matched$true0_samp_times,
    prob_reported_func = rd_fun
  )
  
  all_samp_times_df <- data.frame(
    "samp_time" = tips_matched$true0_samp_times,
    "observed" = !tips_matched$drop
  )
  
  # Drop tips to get observed tree due to RD
  trees$obs <- drop.tip(
    full_tree$newick,
    tip = tips_matched$tip_labels[tips_matched$drop]
  )
  
  # Match observed tree tips to sampling times
  obs_tips_matched <- data.frame(
    connect_sample_to_tips(trees$obs)$ungrouped_samp_df,
    "true0_samp_times" = tips_matched$true0_samp_times[!tips_matched$drop]
  )
  
  # Get truncation sampling time offset
  # Choose tips less than truncation_time to be dropped
  obs_tips_matched$tip_below_trunc <- obs_tips_matched$true0_samp_times < truncation_time
  
  # Drop tips to get truncated tree
  trees[[3]] <- drop.tip(
    trees$obs,
    tip = obs_tips_matched$tip_labels[obs_tips_matched$tip_below_trunc]
  )
  
  # Collect untruncated tips matched to sampling times
  trunc_tips_matched <- data.frame(
    connect_sample_to_tips(trees[[3]])$ungrouped_samp_df,
    "true0_samp_times" = obs_tips_matched$true0_samp_times[!obs_tips_matched$tip_below_trunc]
  )
  
  list(
    "trees" = trees, 
    time0_offsets = list(
      "obs" = min(obs_tips_matched$true0_samp_times), 
      "trunc" = min(trunc_tips_matched$true0_samp_times)
    ),
    all_samp_times_df = all_samp_times_df
  )
}


# Define function to correct BNPR output time 0 ---------------------------
correct_BNPR_time_zero_offset <- function(BNPR_output, time0_offset) {
  adj_BNPR_output <- BNPR_output
  adj_BNPR_output$grid <- BNPR_output$grid + time0_offset
  adj_BNPR_output$x <- BNPR_output$x + time0_offset
  adj_BNPR_output$samp_times <- BNPR_output$samp_times + time0_offset
  adj_BNPR_output$coal_times <- BNPR_output$coal_times + time0_offset
  adj_BNPR_output
}

# Ne(t) Estimation --------------------------------------------------------
est_eff_pops <- function(
    simed_trees_list, 
    historic_reporting_delay_data, 
    lengthout, 
    beta1_mean,
    beta1_prec
  ) {
  
  est_list <- vector("list", length = (length(simed_trees_list$trees) * 2 + 2))
  
  simed_trees_names <- names(simed_trees_list$trees)
  trunc_names <- simed_trees_names[str_detect(
    simed_trees_names, 
    pattern = "trunc"
  )]
  
  names(est_list) <- c(
    "full_nops", "full_ps", 
    "obs_nops", "obs_ps",
    "obs_ps_rd_as_offset", "obs_ps_rd_as_fn",
    paste0(trunc_names, c("_nops", "_ps"))
  )
  
  # Full tree
  est_list$full_nops <- suppressWarnings(phylodyn2::BNPR(
    simed_trees_list$trees$full, 
    lengthout = lengthout
  ))
  # Warning suppressed because we know we have coincident sampling and coalescent times
  
  est_list$full_ps <- suppressWarnings(phylodyn2::BNPR_PS(
    simed_trees_list$trees$full, 
    lengthout = lengthout,
    beta1_mean = beta1_mean,
    beta1_prec = beta1_prec
  ))
  
  # Observed tree due to RD
  est_list$obs_nops <- correct_BNPR_time_zero_offset(
    suppressWarnings(phylodyn2::BNPR(
      simed_trees_list$trees$obs, 
      lengthout = lengthout
    )), 
    time0_offset = simed_trees_list$time0_offsets$obs
  )
  
  est_list$obs_ps <- correct_BNPR_time_zero_offset(
    suppressWarnings(phylodyn2::BNPR_PS(
      simed_trees_list$trees$obs, 
      lengthout = lengthout,
      beta1_mean = beta1_mean, 
      beta1_prec = beta1_prec
    )), 
    time0_offset = simed_trees_list$time0_offsets$obs
  )
  
  # Observed tree due to RD with RD prob offset in estimation
  est_list$obs_ps_rd_as_offset <- correct_BNPR_time_zero_offset(
    suppressWarnings(phylodyn2::BNPR_PS_with_RD(
      data = simed_trees_list$trees$obs,
      rd_as_offset = TRUE, 
      beta1_mean = beta1_mean,
      beta1_prec = beta1_prec,
      historic_reporting_delays = historic_reporting_delay_data,
      lengthout = lengthout,
      time_offset = simed_trees_list$time0_offsets$obs
    )),
    time0_offset = simed_trees_list$time0_offsets$obs
  )
  
  # Observed tree due to RD with RD prob covariate fn in estimation
  est_list$obs_ps_rd_as_fn <- correct_BNPR_time_zero_offset(
    suppressWarnings(phylodyn2::BNPR_PS_with_RD(
      data = simed_trees_list$trees$obs,
      rd_as_offset = FALSE,
      beta1_mean = beta1_mean,
      beta1_prec = beta1_prec,
      historic_reporting_delays = historic_reporting_delay_data,
      lengthout = lengthout,
      time_offset = simed_trees_list$time0_offsets$obs
    )),
    time0_offset = simed_trees_list$time0_offsets$obs
  )
  
  # Obs truncated trees
  est_list[[paste0(trunc_names, "_nops")]] <- correct_BNPR_time_zero_offset(
    suppressWarnings(phylodyn2::BNPR(
      simed_trees_list$trees[[trunc_names]],
      lengthout = lengthout
    )), 
    time0_offset = as.numeric(simed_trees_list$time0_offsets$trunc)
  )
  
  est_list[[paste0(trunc_names, "_ps")]] <- correct_BNPR_time_zero_offset(
    suppressWarnings(phylodyn2::BNPR_PS(
      simed_trees_list$trees[[trunc_names]], 
      beta1_mean = beta1_mean,
      beta1_prec = beta1_prec,
      lengthout = lengthout
    )), 
    time0_offset = as.numeric(simed_trees_list$time0_offsets$trunc)
  )
    
  
  est_list
}

# Summarize and save inference results df --------------------------------
save_eff_pop_est_summary <- function(est_list, eff_pop_estimation_path, sim_iter) {
  make_res_df <- function(est, scenario_type, time_of_sim, sim_iter) {
    data.frame(
      x = est$x,
      est_effpop = est$effpop, 
      est_effpop025 = est$effpop025,
      est_effpop975 = est$effpop975
    ) |> 
      mutate(scenario = scenario_type) |> 
      mutate(time_of_sim = time_of_sim) |> 
      mutate(sim_iter = sim_iter)
  }
  
  est_df <- NULL
  est_names <- names(est_list)
  time_of_sim <- Sys.time()
  
  for (i in 1:length(est_names)) {
    est_df <- bind_rows(
      est_df, 
      make_res_df(
        est_list[[i]], 
        scenario_type = est_names[i], 
        time_of_sim = time_of_sim,
        sim_iter = sim_iter
      ),
    )
  }
  
  if (file.exists(eff_pop_estimation_path)) {
    prev_est_df <- readr::read_csv(eff_pop_estimation_path, show_col_types = FALSE)
    
    est_df <- prev_est_df |> 
      data.frame() |> 
      bind_rows(est_df)
  }
  
  readr::write_csv(est_df, file = eff_pop_estimation_path)
  
}

save_sampling_times <- function(simed_trees_list, scenario_name, samp_times_path) {
  samp_times_df <- simed_trees_list$all_samp_times_df |> 
    mutate(scenario = scenario_name)
  
  write_csv(samp_times_df, samp_times_path)
  
}

# Summarize and save beta1 inference results df --------------------------
save_beta1_est_summary <- function(
    est_list, 
    pref_samp_coeff_estimation_path,
    sim_iter
  ){
  
  beta1_est_df <- NULL
  beta1_summ <- est_list$obs_ps_rd_as_offset$beta1summ
  
  names_ps_est <- names(est_list)[which(!str_detect(
    names(est_list), 
    pattern = "_nops"
  ))]
  
  for (ps_res in names_ps_est) {
    beta1_est_df <- bind_rows(
      beta1_est_df, 
      est_list[[ps_res]]$beta1summ |> 
        rename(
          quant0.025 = '0.025quant', 
          quant0.5 = '0.5quant', 
          quant0.975 = '0.975quant'
        ) |> 
        mutate(scenario = ps_res) |> 
        mutate(sim = sim_iter) |> 
        remove_rownames(), 
    )
  }
  
  if (file.exists(pref_samp_coeff_estimation_path)) {
    prev_beta1_est_df <- readr::read_csv(
      pref_samp_coeff_estimation_path, 
      show_col_types = FALSE
    )
    
    beta1_est_df <- prev_beta1_est_df |> 
      data.frame() |> 
      bind_rows(beta1_est_df)
  }
  
  readr::write_csv(beta1_est_df, file = pref_samp_coeff_estimation_path)
  
}


# Evaluate estimation -----------------------------------------------------
eval_sim_inference <- function(
    eff_pop_estimation_path, 
    eff_pop_traj, 
    ave_window = 7,
    simulation_scenario
  ) {
  
  eff_pop_est_summary <- readr::read_csv(
    eff_pop_estimation_path, 
    show_col_types = FALSE
  ) |> 
    select(-sim_iter, -time_of_sim) |> 
    mutate(true_effpop = eff_pop_traj(x))
  
  max_time <- ceiling(max(eff_pop_est_summary$x))
  
  int_breaks <- c(
    seq(0, max_time, by = ave_window), 
    max_time + ave_window
  )
  
  est_by_interval <- eff_pop_est_summary |> 
    mutate(interval_midpoint = cut(
      x, 
      breaks = int_breaks, 
      labels = as.character(roll::roll_mean(int_breaks, width = 2)[-1]),
      include.lowest = TRUE, 
      right = FALSE
    )) |> 
    mutate(interval_label = cut(
      x, 
      breaks = int_breaks, 
      include.lowest = TRUE, 
      right = FALSE
    )) |> 
    summarize(
      mse = mean((est_effpop - true_effpop)^2),
      mean_relative_deviation = mean((est_effpop - true_effpop) / true_effpop),
      per_coverage = 100 * mean(
        true_effpop >= est_effpop025 & true_effpop <= est_effpop975
      ),
      rev_per_coverage = 100 - per_coverage,
      mean_band_width = mean(est_effpop975 - est_effpop025), 
      .by = c(scenario, interval_midpoint, interval_label)
    ) |> 
    mutate(interval_midpoint = as.character(interval_midpoint)) |> 
    ungroup()
  
  overall_est <- eff_pop_est_summary |> 
    summarize(
      mse = mean((est_effpop - true_effpop)^2),
      mean_relative_deviation = mean((est_effpop - true_effpop) / true_effpop),
      per_coverage = 100 * mean(
        true_effpop >= est_effpop025 & true_effpop <= est_effpop975
      ),
      rev_per_coverage = 100 - per_coverage,
      mean_band_width = mean(est_effpop975 - est_effpop025), 
      .by = c(scenario)
    ) |> 
    mutate(interval_midpoint = "all") |> 
    mutate(interval_label = paste0("[0,", max_time, "]"))
  
  bind_rows(overall_est, est_by_interval) |> 
    mutate(sim_scenario = simulation_scenario)
}
