library(tidyverse)
library(knitr)
library(kableExtra)


# Plot inference on simulated/real data --------------------------------------

prep_est_fits_comparisons_data <- function(
    est_fits,
    time_zero_date = NULL # Null if simulations
) {
  
  prep_single_BNPR_data <- function(
    BNPR_out, traj = NULL, scenario_name
  ) {
    
    BNPR_data <- data.frame(
      t = BNPR_out$x,
      eff_pop_median = BNPR_out$effpop,
      eff_pop_95ub = BNPR_out$effpop975,
      eff_pop_95lb = BNPR_out$effpop025
    )
    
    if (!is.null(traj)) {
      BNPR_data <- BNPR_data |> 
        mutate(eff_pop_true = traj(t))
    }
    
    BNPR_data
    
  }
  
  formatted_data <- map_dfr(
    est_fits, 
    .f = prep_single_BNPR_data, 
    .id = "scenario"
  ) |>
    mutate(data_with_rds = factor(
      case_when(
        scenario %in% c("full_no_ps", "full_ps") ~ "full",
        scenario %in% c("trunc_no_ps", "trunc_ps") ~ "trunc",
        .default = "obs"
      ),
      levels = c("obs", "trunc", "full"),
      labels = c(
        "Real-time Inference with Partially-reported Data",
        "Real-time Inference with Partially-reported, and Truncated Data",
        "Retrospective Inference with Fully-reported Data"
      )
    )) |>
    mutate(model = factor(
      case_when(
        scenario %in% c("full_no_ps", "trunc_no_ps", "obs_no_ps") ~ "BNPR",
        scenario %in% c("full_ps", "trunc_ps", "obs_ps") ~ "BNPR PS",
        scenario == "obs_ps_rd_offset" ~ "Delay-Aware BNPR PS (RP offset implementation)",
        scenario == "obs_ps_rd_covariate" ~ "Delay-Aware BNPR PS (RP covariate implementation)",
        TRUE ~ "error"
      ),
      levels = c(
        "Model without PS", 
        "PS Model", 
        "Delay-Aware BNPR PS (RP offset implementation)", 
        "Delay-Aware BNPR PS (RP covariate implementation)"
      )
    )) |> 
    mutate(scenario = factor(
      scenario,
      levels = c(
        "obs_no_ps",
        "trunc_no_ps", 
        "obs_ps",
        "trunc_ps",
        "obs_ps_rd_offset",
        "obs_ps_rd_covariate",
        "full_no_ps", 
        "full_ps"
      )
    ))
  
  if (!is.null(time_zero_date)) {
    formatted_data <- formatted_data |> 
      mutate(x_date = as_date(-t * 365, origin = time_zero_date))
  }
  
  formatted_data
  
}

plot_single_comparisons_main <- function(
    est_fits,
    truncation_time,
    min_x, # time for simulation or date for real data 
    y_lims,
    font_size,
    eff_pop_traj = NULL,
    time_zero_date_for_real_data = NULL,
    legend_label = "Analysis",
    show_facet_labels = TRUE
) {
  
  if (is.null(time_zero_date_for_real_data)) { # Simulations
    truncation_x <- truncation_time
    
    inference_results <- prep_est_fits_comparisons_data(
      est_fits = est_fits
    ) |> 
      mutate(plot_x_value = t) |> 
      filter(plot_x_value <= min_x)
    
    data_group_details <- data.frame(
      scenario = c("obs_no_ps", "obs_ps", "obs_ps_rd_offset", "full_ps", "truth"),
      color = c("blue", "red", "darkgreen", "black", "black"),
      linetype = c("dashed", "dashed", "dashed", "dashed", "solid"),
      labels = c(
        "Real-Time BNPR", 
        "Real-Time BNPR PS", 
        "Real-Time Delay-Aware BNPR PS", 
        "Retrospective BNPR PS",
        "Truth"
      )
    )
    
    x_tick_label_slant <- 0
    
  } else { # Real data
    truncation_x <- time_zero_date_for_real_data - days(truncation_time) 
    
    inference_results <- prep_est_fits_comparisons_data(
      est_fits = est_fits,
      time_zero_date = time_zero_date_for_real_data
    ) |> 
      mutate(plot_x_value = x_date) |> 
      filter(plot_x_value >= min_x)
    
    data_group_details <- data.frame(
      scenario = c("obs_no_ps", "obs_ps", "obs_ps_rd_offset", "full_ps"),
      color = c("blue", "red", "darkgreen", "black"),
      linetype = c("dashed", "dashed", "dashed", "dashed"),
      labels = c(
        "Real-Time BNPR", 
        "Real-Time BNPR PS", 
        "Real-Time Delay-Aware BNPR PS", 
        "Retrospective BNPR PS"
      )
    )
    
    x_tick_label_slant <- 45
    
  }
  
  subsetted_scenarios <- c(
    "obs_no_ps",
    "obs_ps", 
    "obs_ps_rd_offset" 
  )
  
  # named vectors for scale_._manual()
  group_col_vec <- data_group_details$color
  names(group_col_vec) <- data_group_details$scenario
  
  group_linetype_vec <- data_group_details$linetype
  names(group_linetype_vec) <- data_group_details$scenario
  
  group_labels_vec <- data_group_details$labels
  names(group_labels_vec) <- data_group_details$scenario
  
  # We want this one plotted in every facet of ps
  inference_results_full_ps_only <- inference_results |>
    filter(scenario == "full_ps") |> 
    select(-scenario)
  
  final_plot <- inference_results |>
    filter(scenario %in% subsetted_scenarios) |>
    ggplot(aes(x = plot_x_value, y = eff_pop_median, color = scenario)) +
    geom_line(aes(linetype = scenario)) +
    geom_ribbon(
      aes(
        ymin = eff_pop_95lb, 
        ymax = eff_pop_95ub, 
        color = scenario, 
        fill = scenario
      ), 
      alpha = 0.4, 
      colour = NA
    ) +
    geom_line(
      data = inference_results_full_ps_only,
      aes(
        x = plot_x_value, 
        y = eff_pop_median, 
        color = "full_ps", 
        linetype = "full_ps"
      )
    ) +
    geom_ribbon(
      data = inference_results_full_ps_only,
      aes(
        ymin = eff_pop_95lb, 
        ymax = eff_pop_95ub, 
        color = "full_ps", 
        fill = "full_ps"
      ), 
      alpha = 0.3, 
      colour = NA
    ) +
    annotate(
      'rect',
      xmin = min_x, xmax = truncation_x,
      ymin = -Inf, ymax = Inf,
      fill = '#00000010'
    ) +
    # geom_text(
    #   data = filter(data_group_details, !(scenario %in% c("full_ps", "truth"))),
    #   aes(
    #     x = 80, 
    #     y = 275, 
    #     label = scenario, 
    #     color = scenario
    #   ),
    #   size = font_size,
    #   hjust = 0
    # ) +
    labs(
      y = "Effective Population Size",
      color = legend_label,
      fill = legend_label,
      linetype = legend_label
    ) +
    scale_color_manual(
      values = group_col_vec, 
      labels = group_labels_vec
    ) +
    scale_fill_manual(
      values = group_col_vec,
      labels = group_labels_vec
    ) +
    scale_linetype_manual(
      values = group_linetype_vec, 
      labels = group_labels_vec
    ) +
    coord_cartesian(ylim = y_lims) +
    guides(
      color = guide_legend(
        title.position = "top", 
        title.hjust = 0.5, 
        nrow = 2, byrow = TRUE
      ),
      fill = guide_legend(
        title.position = "top", 
        title.hjust = 0.5, 
        nrow = 2, byrow = TRUE
      ),
      linetype = guide_legend(
        title.position = "top", 
        title.hjust = 0.5, 
        nrow = 2, byrow = TRUE
      )
    ) +
    theme_bw(base_size = font_size) +
    theme(
      axis.text.x = element_text(
        angle = x_tick_label_slant, 
        size = font_size, 
        hjust = 1
      ),
      axis.title.x = element_text(size = font_size),
      axis.text.y = element_text(size = font_size),
      axis.title.y = element_text(size = font_size),
      text = element_text(size = font_size),
      legend.text = element_text(size = font_size),
      legend.position = "bottom",
      legend.direction = "horizontal",
      strip.background = element_blank(),
      legend.key.size = unit(1, "lines")
    ) + 
    facet_wrap(
      ~  factor(
        scenario, 
        levels = c(
          "obs_no_ps",
          #"trunc_no_ps",
          "obs_ps",
          #"trunc_ps",
          "obs_ps_rd_offset"#,
          #"obs_ps_rd_covariate"
        ), 
        labels = c(
          "Real-Time BNPR\nVersus Retrospective BNPR PS",
          #"Real-time BNPR\nversus retrospective BNPR PS ",
          "Real-Time BNPR PS\nVersus Retrospective BNPR PS",
          #"Real-time BNPR PS\nversus retrospective BNPR PS ",
          "Real-Time Delay-Aware BNPR PS\nVersus Retrospective BNPR PS"#,
          #"Real-time BNPR PS with RP covariate\nversus retrospective BNPR PS"
        )
      ), 
      nrow = 1
    )
  
  if (is.null(time_zero_date_for_real_data)) { # Real data
    true_eff_pop_df <- data.frame(
      plot_x_value = seq(0, min_x, by = 0.01)
    ) |> 
      mutate(eff_pop_median = eff_pop_traj(plot_x_value))
    
    final_plot <- final_plot +
      xlab("Time Prior to Present (days)") +
      scale_x_reverse(expand = c(0, 0)) +
      geom_line(
        data = true_eff_pop_df,
        aes(
          x = plot_x_value,
          y = eff_pop_median,
          color = "truth",
          fill = "truth",
          linetype = "truth"
        )
      )
    
  } else {
    final_plot <- final_plot +
      scale_x_date(
        date_breaks = "1 months", 
        date_labels = "%b %d, %Y", 
        expand = c(0, 0)
      ) +
      xlab("Date")
    
  }
  
  if(show_facet_labels) {
    final_plot <- final_plot +
      theme(strip.text.x = element_text(size = font_size-1))
      
  } else {
    final_plot <- final_plot +
      theme(strip.text.x = element_blank())
  }
  
  final_plot
  
}

plot_single_comparisons_supplemental <- function(
    est_fits,
    truncation_time,
    min_x, # oldest time for simulation or date for real data 
    y_lims,
    eff_pop_traj = NULL,
    time_zero_date_for_real_data = NULL,
    legend_label = "Inference"
) {
  
  if (is.null(time_zero_date_for_real_data)) { # Simulations
    truncation_x <- truncation_time
    
    inference_results <- prep_est_fits_comparisons_data(
      est_fits = est_fits
    ) |> 
      mutate(plot_x_value = t) |> 
      filter(plot_x_value <= min_x)
    
    data_group_details <- data.frame(
      scenario = c("obs", "trunc", "full", "truth"),
      color = c("darkorange2", "purple", "black", "black"),
      linetype = c("dashed", "dashed", "dashed", "solid"),
      labels = c(levels(inference_results$data_with_rds), "truth")
    )
    
  } else { # Real data
    truncation_x <- time_zero_date_for_real_data - days(truncation_time) 
    
    inference_results <- prep_est_fits_comparisons_data(
      est_fits = est_fits,
      time_zero_date = time_zero_date_for_real_data
    ) |> 
      mutate(plot_x_value = x_date) |> 
      filter(plot_x_value >= min_x)
    
    data_group_details <- data.frame(
      scenario = c("obs", "trunc", "full"),
      color = c("darkorange2", "purple", "black"),
      linetype = c("dashed", "dashed", "dashed"),
      labels = levels(inference_results$data_with_rds)
    )
    
  }
  
  subsetted_scenarios <- c(
    "obs_no_ps", "trunc_no_ps",
    "obs_ps", 
    "trunc_ps", 
    "obs_ps_rd_offset", 
    "obs_ps_rd_covariate"
  )
  
  # named vectors for scale_._manual()
  group_col_vec <- data_group_details$color
  #group_col_vec <- my_colors
  names(group_col_vec) <- data_group_details$labels
  
  group_linetype_vec <- data_group_details$linetype
  names(group_linetype_vec) <- data_group_details$labels
  
  group_labels_vec <- data_group_details$labels
  names(group_labels_vec) <- data_group_details$labels
  
  # We want this one plotted in every facet of no ps
  results_full_no_ps_only <- inference_results |>
    filter(scenario == "full_no_ps")
    
  results_full_no_ps_for_obs_no_ps <- results_full_no_ps_only |> 
    mutate(scenario = "obs_no_ps")
  
  results_full_no_ps_for_trunc_no_ps <- results_full_no_ps_only |> 
    mutate(scenario = "trunc_no_ps")
  
  # We want this one plotted in every facet of ps
  results_full_ps_only <- inference_results |>
    filter(scenario == "full_ps")
  
  results_full_ps_for_obs_ps_only <- results_full_ps_only |>
    mutate(scenario = "obs_ps")
  
  results_full_ps_for_trunc_ps_only <- results_full_ps_only |>
    mutate(scenario = "trunc_ps")
  
  results_full_ps_for_obs_ps_rd_offset_only <- results_full_ps_only |>
    mutate(scenario = "obs_ps_rd_offset")
  
  results_full_ps_for_obs_ps_rd_covariate_only <- results_full_ps_only |>
    mutate(scenario = "obs_ps_rd_covariate")
  
  final_plot <- inference_results |>
    filter(scenario %in% subsetted_scenarios) |>
    ggplot(
      aes(x = plot_x_value, y = eff_pop_median, color = data_with_rds)
    ) +
    geom_line(aes(linetype = data_with_rds)) +
    geom_ribbon(
      aes(
        ymin = eff_pop_95lb, 
        ymax = eff_pop_95ub, 
        color = data_with_rds, 
        fill = data_with_rds
      ), 
      alpha = 0.4, 
      colour = NA
    ) +
    geom_line(
      data = results_full_no_ps_for_obs_no_ps,
      aes(
        x = plot_x_value, 
        y = eff_pop_median, 
        color = data_group_details$labels[data_group_details$scenario == "full"], 
        linetype = data_group_details$labels[data_group_details$scenario == "full"]
      )
    ) +
    geom_ribbon(
      data = results_full_no_ps_for_obs_no_ps,
      aes(
        ymin = eff_pop_95lb, 
        ymax = eff_pop_95ub, 
        color = data_group_details$labels[data_group_details$scenario == "full"], 
        fill = data_group_details$labels[data_group_details$scenario == "full"]
      ), 
      alpha = 0.3, 
      colour = NA
    ) +
    geom_line(
      data = results_full_no_ps_for_trunc_no_ps,
      aes(
        x = plot_x_value,
        y = eff_pop_median,
        color = data_group_details$labels[data_group_details$scenario == "full"],
        linetype = data_group_details$labels[data_group_details$scenario == "full"]
      )
    ) +
    geom_ribbon(
      data = results_full_no_ps_for_trunc_no_ps,
      aes(
        ymin = eff_pop_95lb,
        ymax = eff_pop_95ub,
        color = data_group_details$labels[data_group_details$scenario == "full"],
        fill = data_group_details$labels[data_group_details$scenario == "full"]
      ),
      alpha = 0.3, 
      colour = NA
    ) +
    geom_line(
      data = results_full_ps_for_obs_ps_only,
      aes(
        x = plot_x_value,
        y = eff_pop_median,
        color = data_group_details$labels[data_group_details$scenario == "full"],
        linetype = data_group_details$labels[data_group_details$scenario == "full"]
      )
    ) +
    geom_ribbon(
      data = results_full_ps_for_obs_ps_only,
      aes(
        ymin = eff_pop_95lb,
        ymax = eff_pop_95ub,
        color = data_group_details$labels[data_group_details$scenario == "full"],
        fill = data_group_details$labels[data_group_details$scenario == "full"]
      ),
      alpha = 0.3, 
      colour = NA
    ) +
    geom_line(
      data = results_full_ps_for_trunc_ps_only,
      aes(
        x = plot_x_value,
        y = eff_pop_median,
        color = data_group_details$labels[data_group_details$scenario == "full"],
        linetype = data_group_details$labels[data_group_details$scenario == "full"]
      )
    ) +
    geom_ribbon(
      data = results_full_ps_for_trunc_ps_only,
      aes(
        ymin = eff_pop_95lb,
        ymax = eff_pop_95ub,
        color = data_group_details$labels[data_group_details$scenario == "full"],
        fill = data_group_details$labels[data_group_details$scenario == "full"]
      ),
      alpha = 0.3, 
      colour = NA
    ) +
    geom_line(
      data = results_full_ps_for_obs_ps_rd_offset_only,
      aes(
        x = plot_x_value,
        y = eff_pop_median,
        color = data_group_details$labels[data_group_details$scenario == "full"],
        linetype = data_group_details$labels[data_group_details$scenario == "full"]
      )
    ) +
    geom_ribbon(
      data = results_full_ps_for_obs_ps_rd_offset_only,
      aes(
        ymin = eff_pop_95lb,
        ymax = eff_pop_95ub,
        color = data_group_details$labels[data_group_details$scenario == "full"],
        fill = data_group_details$labels[data_group_details$scenario == "full"]
      ),
      alpha = 0.3, 
      colour = NA
    ) +
    geom_line(
      data = results_full_ps_for_obs_ps_rd_covariate_only,
      aes(
        x = plot_x_value,
        y = eff_pop_median,
        color = data_group_details$labels[data_group_details$scenario == "full"],
        linetype = data_group_details$labels[data_group_details$scenario == "full"]
      )
    ) +
    geom_ribbon(
      data = results_full_ps_for_obs_ps_rd_covariate_only,
      aes(
        ymin = eff_pop_95lb,
        ymax = eff_pop_95ub,
        color = data_group_details$labels[data_group_details$scenario == "full"],
        fill = data_group_details$labels[data_group_details$scenario == "full"]
      ),
      alpha = 0.3, 
      colour = NA
    ) +
    annotate(
      'rect',
      xmin = min_x, xmax = truncation_x,
      ymin = -Inf, ymax = Inf,
      fill = '#00000010'
    ) +
    labs(
      y = "Effective Population Size",
      color = legend_label,
      fill = legend_label,
      linetype = legend_label
    ) +
    #ylim(0, 50) +
    scale_color_manual(
      values = group_col_vec, 
      labels = group_labels_vec
    ) +
    scale_fill_manual(
      values = group_col_vec,
      labels = group_labels_vec
    ) +
    scale_linetype_manual(
      values = group_linetype_vec, 
      labels = group_labels_vec
    ) +
    coord_cartesian(ylim = y_lims) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.title=element_blank()
    ) + 
    facet_wrap(
      ~  factor(
        scenario, 
        levels = c(
          "obs_no_ps",
          "trunc_no_ps",
          "obs_ps",
          "trunc_ps",
          "obs_ps_rd_offset",
          "obs_ps_rd_covariate"
        ), 
        labels = c(
          "BNPR Model",
          "BNPR Model ",
          "BNPR PS Model",
          "BNPR PS Model ",
          "Delay-Aware BNPR PS Model\nwith (RP offset implementation)",
          "Delay-Aware BNPR PS Model\nwith (RP covariate implementation)"
        )
      ),
      nrow = 3
    ) +
    guides(
      # all layers stack, so it gets darker with each layer
      fill = guide_legend(override.aes = c(alpha = 0.075)) 
    )
  
  if (is.null(time_zero_date_for_real_data)) { # Simulations
    true_eff_pop_df <- data.frame(
      plot_x_value = seq(0, min_x, by = 0.01)
    ) |> 
      mutate(eff_pop_median = eff_pop_traj(plot_x_value))
    
    final_plot <- final_plot +
      xlab("Time Since Most Recently Collected Sample (days)") +
      scale_x_reverse(expand = c(0, 0)) +
      geom_line(
        data = true_eff_pop_df,
        aes(
          x = plot_x_value,
          y = eff_pop_median,
          color = data_group_details$labels[data_group_details$scenario == "truth"],
          fill = data_group_details$labels[data_group_details$scenario == "truth"],
          linetype = data_group_details$labels[data_group_details$scenario == "truth"]
        )
      )
      
  } else { # Real data analysis
    final_plot <- final_plot +
      scale_x_date(
        date_breaks = "1 months", 
        date_labels = "%b %d, %Y", 
        expand = c(0, 0)
      ) +
      xlab("Date")
    
  }
  
  final_plot
  
}

# Plot inference evaluation ----------------------------------------------
plot_main_sim_inference_eval <- function(
    all_est_evaluation, 
    truncation_time, 
    x_lim,
    font_size
  ) {
  
  evaluation_names_df <- data.frame(
    level = c(
      #"mse", 
      "mean_relative_deviation",
      "per_coverage", 
      "mean_band_width"
    ),
    label = c(
      #"Mean Square Error", 
      "Mean\nRelative Devation",
      "Mean % \nCoverage by BCI", 
      "Mean\n95% BCI Width"
    )
  )
  
  estimation_names_df <- data.frame(
    level = c(
      "obs_nops", 
      "obs_ps", 
      "obs_ps_rd_as_offset"
    ),
    label = c(
      "BNPR", 
      "BNPR PS", 
      "Delay-Aware BNPR PS"
    ),
    linetype = c("dotted", "dotted", "dotted"),
    shape = c(24, 25, 23),
    color = c("blue", "red", "darkgreen"),
    fill = c("blue", "red", "darkgreen")
  )
    
  
  plot_data <- all_est_evaluation |> 
    filter(scenario %in% c("obs_nops", "obs_ps", "obs_ps_rd_as_offset")) |> 
    filter(interval_midpoint != "all") |> 
    mutate(interval_midpoint = as.numeric(interval_midpoint)) |> 
    select(-rev_per_coverage) |> 
    select(-mse) |> 
    pivot_longer(
      cols = c(mean_relative_deviation, per_coverage, mean_band_width),
      names_to = "evaluation",
      values_to = "value"
    ) |> 
    mutate(evaluation = factor(
      evaluation,
      levels = evaluation_names_df$level,
      labels = evaluation_names_df$label
    )) |> 
    mutate(estimation = factor(
      scenario,
      levels = estimation_names_df$level,
      labels = estimation_names_df$label
    ))
  
  data_hline <- data.frame(
    evaluation = factor(evaluation_names_df$label), 
    hline = c(0, 95, NA)
  )
  
  plot_data |> 
    filter(interval_midpoint < 100) |> 
    ggplot(aes(
      x = interval_midpoint, 
      y = value, 
      linetype = estimation, 
      shape = estimation, 
      fill = estimation, 
      color = estimation
    )) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = estimation_names_df$linetype) +
    scale_shape_manual(values = estimation_names_df$shape) +
    scale_fill_manual(values = estimation_names_df$fill) +
    scale_color_manual(values = estimation_names_df$color) +
    facet_grid(evaluation ~ sim_scenario, scales = "free_y", switch = "y") +
    geom_hline(data = data_hline, aes(yintercept = hline), linetype = "dashed") +
    annotate(
      'rect',
      xmin = truncation_time, xmax = Inf,
      ymin = -Inf, ymax = Inf,
      fill = '#00000010'
    ) +
    theme_bw(base_size = font_size) +
    theme(
      strip.background = element_blank(), 
      strip.placement = "outside",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      #text = element_text(family = "A"),
      plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
      legend.position = "bottom",
      panel.spacing = unit(0, "cm"),
      axis.title.x.top = element_text(size = plot_font_size),
      axis.text.x = element_text(size = plot_font_size),
      axis.title.x = element_text(size = plot_font_size),
      axis.text.y = element_text(size = plot_font_size),
      axis.title.y = element_text(size = plot_font_size, margin = margin(l = 0, r = 0)),
      text = element_text(size = plot_font_size),
      legend.text = element_text(size = plot_font_size),
      strip.text = element_text(size = plot_font_size)
    ) +
    guides(
      linetype = guide_legend(title.position = "top", title.hjust = 0.5),
      shape = guide_legend(title.position = "top", title.hjust = 0.5),
      fill = guide_legend(title.position = "top", title.hjust = 0.5),
      color = guide_legend(title.position = "top", title.hjust = 0.5)
    ) + 
    scale_x_reverse() +
    labs(
      x = "Time Since Most Recently Collected Sample (days)",
      y = "",
      linetype = "Real-Time Inferential Method",
      shape = "Real-Time Inferential Method",
      fill = "Real-Time Inferential Method",
      color = "Real-Time Inferential Method"
    )
  
}

plot_all_real_time_sim_inference_eval <- function(
    all_est_evaluation, 
    truncation_time, 
    x_lim
) {
  
  evaluation_names_df <- data.frame(
    level = c(
      #"mse", 
      "mean_relative_deviation",
      "per_coverage", 
      "mean_band_width"
    ),
    label = c(
      #"Mean Square Error", 
      "Mean\nRelative Devation",
      "Mean % \nCoverage by BCI", 
      "Mean\n95% BCI Width"
    )
  )
  
  estimation_names_df <- data.frame(
    level = c(
      str_subset(unique(all_est_evaluation$scenario), pattern = "trunc"),
      "obs_nops",
      "obs_ps",
      "obs_ps_rd_as_offset",
      "obs_ps_rd_as_fn"
    ),
    label = c(
      "BNPR with Truncated Data",
      "BNPR PS with Truncated Data",
      "BNPR",
      "BNPR PS",
      "Delay-Aware BNPR PS (RP offset implementation)",
      "Delay-Aware BNPR PS (RP covariate implementation)"
    ),
    linetype = c("dotted", "dotted", "dotted", "dotted", "dotted", "dotted"),
    shape = c(2, 6, 24, 25, 23, 22),
    color = c("brown", "limegreen", "blue", "#DB3A07FF", "darkgreen", "purple"),
    fill = c("brown", "limegreen", "blue", "#DB3A07FF", "darkgreen", "purple")
  )
  
  plot_data <- all_est_evaluation |> 
    filter(interval_midpoint != "all") |> 
    mutate(interval_midpoint = as.numeric(interval_midpoint)) |> 
    select(-rev_per_coverage) |> 
    select(-mse) |> 
    filter(!(scenario %in% c("full_nops", "full_ps"))) |> 
    pivot_longer(
      cols = c(mean_relative_deviation, per_coverage, mean_band_width),
      names_to = "evaluation",
      values_to = "value"
    ) |> 
    mutate(evaluation = factor(
      evaluation,
      levels = evaluation_names_df$level,
      labels = evaluation_names_df$label
    )) |> 
    mutate(estimation = factor(
      scenario,
      levels = estimation_names_df$level,
      labels = estimation_names_df$label
    ))
  
  data_hline <- data.frame(
    evaluation = factor(evaluation_names_df$label), 
    hline = c(0, 95, NA)
  )
  
  plot_data |> 
    filter(interval_midpoint <= 100) |> 
    ggplot(aes(
      x = interval_midpoint, 
      y = value, 
      linetype = estimation, 
      shape = estimation, 
      fill = estimation, 
      color = estimation
    )) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = estimation_names_df$linetype) +
    scale_shape_manual(values = estimation_names_df$shape) +
    scale_fill_manual(values = estimation_names_df$fill) +
    scale_color_manual(values = estimation_names_df$color) +
    facet_grid(evaluation ~ sim_scenario, scales = "free_y", switch = "y") +
    geom_hline(data = data_hline, aes(yintercept = hline), linetype = "dashed") +
    annotate(
      'rect',
      xmin = truncation_time, xmax = Inf,
      ymin = -Inf, ymax = Inf,
      fill = '#00000010'
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(), 
      strip.placement = "outside",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      #text = element_text(family = "A"),
      plot.margin = unit(c(0.5, 0.5, 0, 0), "pt"),
      axis.title.y = element_text(margin = margin(l = 0, r = 0)),
      legend.position = "bottom",
      panel.spacing = unit(0, "cm"),
    ) +
    guides(
      linetype = guide_legend(title.position = "top", title.hjust = 0.5),
      shape = guide_legend(title.position = "top", title.hjust = 0.5),
      fill = guide_legend(title.position = "top", title.hjust = 0.5),
      color = guide_legend(title.position = "top", title.hjust = 0.5)
    ) + 
    scale_x_reverse() +
    labs(
      x = "Time Since Most Recently Collected Sample (days)",
      y = "",
      linetype = "Real-Time Inferential Method",
      shape = "Real-Time Inferential Method",
      fill = "Real-Time Inferential Method",
      color = "Real-Time Inferential Method"
    )
  
}

plot_all_retrospective_sim_inference_eval <- function(
    all_est_evaluation, 
    truncation_time, 
    x_lim
) {
  
  evaluation_names_df <- data.frame(
    level = c(
      #"mse", 
      "mean_relative_deviation",
      "per_coverage", 
      "mean_band_width"
    ),
    label = c(
      #"Mean Square Error", 
      "Mean\nRelative Devation",
      "Mean % \nCoverage by BCI", 
      "Mean\n95% BCI Width"
    )
  )
  
  estimation_names_df <- data.frame(
    level = c(
      "full_nops",
      "full_ps"
    ),
    label = c(
      "BNPR",
      "BNPR PS"
    ),
    linetype = c("dotted", "dotted"),
    shape = c(24, 25),
    color = c("peru", "cornflowerblue"),
    fill = c("peru", "cornflowerblue")
  )
  
  plot_data <- all_est_evaluation |> 
    filter(interval_midpoint != "all") |> 
    mutate(interval_midpoint = as.numeric(interval_midpoint)) |> 
    select(-rev_per_coverage) |> 
    select(-mse) |> 
    filter(scenario %in% c("full_nops", "full_ps")) |> 
    pivot_longer(
      cols = c(mean_relative_deviation, per_coverage, mean_band_width),
      names_to = "evaluation",
      values_to = "value"
    ) |> 
    mutate(evaluation = factor(
      evaluation,
      levels = evaluation_names_df$level,
      labels = evaluation_names_df$label
    )) |> 
    mutate(estimation = factor(
      scenario,
      levels = estimation_names_df$level,
      labels = estimation_names_df$label
    ))
  
  data_hline <- data.frame(
    evaluation = factor(evaluation_names_df$label), 
    hline = c(0, 95, NA)
  )
  
  plot_data |> 
    filter(interval_midpoint <= 100) |> 
    ggplot(aes(
      x = interval_midpoint, 
      y = value, 
      linetype = estimation, 
      shape = estimation, 
      fill = estimation, 
      color = estimation
    )) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = estimation_names_df$linetype) +
    scale_shape_manual(values = estimation_names_df$shape) +
    scale_fill_manual(values = estimation_names_df$fill) +
    scale_color_manual(values = estimation_names_df$color) +
    facet_grid(evaluation ~ sim_scenario, scales = "free_y", switch = "y") +
    geom_hline(data = data_hline, aes(yintercept = hline), linetype = "dashed") +
    annotate(
      'rect',
      xmin = truncation_time, xmax = Inf,
      ymin = -Inf, ymax = Inf,
      fill = '#00000010'
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(), 
      strip.placement = "outside",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      #text = element_text(family = "A"),
      plot.margin = unit(c(0.5, 0.5, 0, 0), "pt"),
      axis.title.y = element_text(margin = margin(l = 0, r = 0)),
      legend.position = "bottom",
      panel.spacing = unit(0, "cm"),
    ) +
    guides(
      linetype = guide_legend(title.position = "top", title.hjust = 0.5),
      shape = guide_legend(title.position = "top", title.hjust = 0.5),
      fill = guide_legend(title.position = "top", title.hjust = 0.5),
      color = guide_legend(title.position = "top", title.hjust = 0.5)
    ) + 
    scale_x_reverse() +
    labs(
      x = "Time Since Most Recently Collected Sample (days)",
      y = "",
      linetype = "Retrospective Inferential Method",
      shape = "Retrospective Inferential Method",
      fill = "Retrospective Inferential Method",
      color = "Retrospective Inferential Method"
    )
  
}

# Plot preferential sampling coefficient ------------------------------------
plot_pref_samp_coeff_est <- function(pref_samp_paths) {
  beta1_est_df_a <- readr::read_csv(
    file = pref_samp_paths[1],
    show_col_types = FALSE
  ) |> 
    mutate(sim_scenario = "Scenario A")
  
  beta1_est_df_b <- readr::read_csv(
    file = pref_samp_paths[2],
    show_col_types = FALSE
  ) |> 
    mutate(sim_scenario = "Scenario B")
  
  beta1_est_df_c <- readr::read_csv(
    file = pref_samp_paths[3],
    show_col_types = FALSE
  ) |> 
    mutate(sim_scenario = "Scenario C")
  
  beta1_est_df <- beta1_est_df_a |> 
    bind_rows(
      beta1_est_df_b
    ) |> 
    bind_rows(
      beta1_est_df_c
    ) |> 
    filter(stringr::str_detect(scenario, pattern = "_ps")) |> 
    mutate(scenario = ifelse(
      stringr::str_detect(scenario, pattern = "trunc"),
      yes = "trunc_ps",
      no = scenario
    )) |> 
    mutate(scenario = factor(
      scenario,
      levels = c(
        "full_ps", 
        "obs_ps", 
        "obs_ps_rd_as_offset", 
        "obs_ps_rd_as_fn",
        "trunc_ps"
      ),
      labels = c(
        "Retrospective BNPR PS", 
        "Real-Time BNPR PS", 
        "Real-Time Delay Aware BNPR PS\n(RP offset implementation)", 
        "Real-Time Delay Aware BNPR PS\n(RP covariate implementation)",
        "Real-Time BNPR PS\nwith Truncated Data"
      )
    ))
  
  ggplot(beta1_est_df, aes(x = scenario, y = quant0.5)) +
    geom_boxplot() +
    facet_wrap(~ sim_scenario) +
    theme_bw() +
    geom_hline(yintercept = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
      x = "Inferential Strategy",
      y = "Preferential Sampling\nCoefficient Median Estimate",
      title = "Boxplot of Posterior Median of Preferential Sampling Coefficient,\nby Inferential Strategy and Simulation Scenario"
    )
}

# Sim evaluation tables ---------------------------------------------------


write_eval_table <- function(
    all_est_evaluation,
    est_eval_table_path,
    eval_metric, 
    ave_window, 
    truncation_time
) {
  
  est_evaluation_intervals <- all_est_evaluation |>
    filter(interval_midpoint != "all") |>
    mutate(interval_midpoint = as.numeric(interval_midpoint))
  
  interval_mids <- est_evaluation_intervals |>
    filter(interval_midpoint <= truncation_time) |>
    pull(interval_midpoint) |>
    unique() |>
    sort()
  
  intervals <- paste0(
    "[", interval_mids - ave_window / 2, ", ",
    interval_mids + ave_window / 2, ")"
  )
  
  table_eval_general <- all_est_evaluation |> 
    filter(interval_midpoint %in% c(as.character(interval_mids), "all")) |> 
    mutate(interval_midpoint = as.numeric(replace(
      interval_midpoint, 
      interval_midpoint == "all", 
      "99999"
    ))) |>
    select(interval_midpoint, interval_label, sim_scenario, scenario, all_of(eval_metric)) |> 
    pivot_wider(
      names_from = scenario, 
      names_glue = "{scenario}",
      values_from = eval_metric
    ) |> 
    arrange(sim_scenario, interval_midpoint)
  
  table_eval_scenario_A <- table_eval_general |> 
    filter(sim_scenario == "Scenario A") |> 
    #select(-sim_scenario) |> 
    filter(!is.na(obs_nops))
  
  table_eval_scenario_B <- table_eval_general |> 
    filter(sim_scenario == "Scenario B") |> 
    # select(-sim_scenario) |> 
    filter(!is.na(obs_nops))
  
  table_eval_scenario_C <- table_eval_general |> 
    filter(sim_scenario == "Scenario C") |> 
    #select(-sim_scenario) |> 
    filter(!is.na(obs_nops))
  
  full_table <- table_eval_scenario_A |> 
    bind_rows(table_eval_scenario_B) |> 
    bind_rows(table_eval_scenario_C) |> 
    select(
      sim_scenario,
      interval_label,
      full_nops, full_ps,
      str_subset(unique(all_est_evaluation$scenario), pattern = "trunc"), 
      obs_nops, obs_ps, 
      obs_ps_rd_as_offset, obs_ps_rd_as_fn
    )
  
  write_csv(full_table, est_eval_table_path)
  
}

print_eval_table <- function(
    est_eval_table_path,
    num_digits,
    row_nums = c(1, 9, 10, 15, 16, 22)
) {
  
  full_table <- read_csv(est_eval_table_path, show_col_types = FALSE)
  
  options(knitr.kable.NA = '')
  
  full_table |> 
    select(-sim_scenario) |> 
    kable(
      booktabs = TRUE,
      col.names = c(
        "Time Period (days)", 
        "BNPR", "BNPR PS", 
        "Trunc. BNPR", "Tunc. BNPR PS", 
        "BNPR", "BNPR PS",
        "Delay-Aware BNPR PS (RP offset implementation)", "BNPR PS (RP covariate implementation)"
      ), 
      digits = num_digits,
      escape = FALSE
    ) |> 
    kable_styling(
      font_size = 7,
      latex_options = c("striped", "scale_down")
    ) |> 
    add_header_above(c(" ", "Retrospective" = 2, "Real-Time Inference" = 6)) |> 
    pack_rows("Scenario A", row_nums[1], row_nums[2], underline = FALSE) |> 
    pack_rows("Scenario B", row_nums[3], row_nums[4], underline = FALSE) |> 
    pack_rows("Scenario C", row_nums[5], row_nums[6], underline = FALSE) |> 
    column_spec(1, width = "0.5in") |> 
    column_spec(2:9, width = "0.325in")
}

write_eval_reduced_table <- function(full_table_paths, reduced_table_path){
  
  tab1 <- read_csv(full_table_paths[1], show_col_types = FALSE) |> 
    select(-full_nops, -full_ps, -trunc55_nops, -trunc55_ps)
  
  tab2 <- read_csv(full_table_paths[2], show_col_types = FALSE) |> 
    select(-full_nops, -full_ps, -trunc55_nops, -trunc55_ps)
  
  tab3 <- read_csv(full_table_paths[3], show_col_types = FALSE) |> 
    select(-full_nops, -full_ps, -trunc55_nops, -trunc55_ps)
  
  reduced_table <- data.frame(
    tab1, 
    select(tab2, -interval_label, -sim_scenario),
    select(tab3, -interval_label, -sim_scenario)
  ) |> 
    filter(!(interval_label %in% c("[0,156]", "[0,229]", "[0,307]")))
  
  write_csv(reduced_table, reduced_table_path)
  
}