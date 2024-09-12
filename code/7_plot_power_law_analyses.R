# Create plots from the power-law analysis done in 6_model_cq_power_law_no3_bdc.R

### NOTE: THIS SCRIPT IS ONCOMPLETE ###

# Load libraries ----
  library("tidyverse")    # general workhorse
  library("here")         # takes the guesswork out of dealing with file paths
  library("patchwork")    # makes laying out plots much easier
  library("grid")         # more plotting assistance
  library("scales")       # used for plotting on log scales
  library("scico")        # color palette for ggplot2


# Load results from 6_model_cq_power_law_XXX_XXX.R
  # CV results
  dat_path <- here("data", "out")
  cv_results <- 
    list.files(path = dat_path, pattern = "cv", full.names = TRUE) %>% 
    # Without set_names(), .id= will use integer indicators, instead of actual file names.
    set_names() %>% 
    map_dfr(read_csv, .id = "filename") %>% 
    # Keep just the short filename w/o full path
    mutate(filename = basename(filename)) %>% 
    # Extract site abbrev from filename
    mutate(site = substr(filename, nchar(filename) - 7 + 1,  nchar(filename) - 4)) %>% 
    # Extract solute abbrev from filename
    mutate(solute = substr(filename, nchar(filename) - 11 + 1,  nchar(filename) - 8)) %>% 
    # Create n source column
    mutate(n_source = case_when(
      grepl("nFit", eq) ~ "fit",
      .default = "input"
    ))
  
  # Model residuals
    # Using C0 from baserlow
  dat_path <- here("data", "out")
  resids <- 
    list.files(path = dat_path, pattern = "resids_baseflow", full.names = TRUE) %>% 
    # Without set_names(), .id= will use integer indicators, instead of actual file names.
    set_names() %>% 
    map_dfr(read_csv, .id = "filename") %>% 
    # Keep just the short filename w/o full path
    mutate(filename = basename(filename)) %>% 
    # Extract site abbrev from filename
    mutate(site = substr(filename, nchar(filename) - 7 + 1,  nchar(filename) - 4)) %>% 
    # Extract solute abbrev from filename
    mutate(solute = substr(filename, nchar(filename) - 11 + 1,  nchar(filename) - 8)) %>% 
    # Create n source column
    mutate(n_source = case_when(
      grepl("nFit", model) ~ "fit",
      .default = "input"
    ))

# Plot cross-validation results ----
  # cat(paste(shQuote(unique(cv_results$eq), type="cmd"), collapse=", "))
  
  ## baseflow vs. wetdep ----
  ## Compare C-V for C0 from baseflow vs. wetdep for models w/ C0 parameter using 15-min data
  cv_results %>% 
    filter(!is.na(eq)) %>% 
    filter(eq %in% c("e2", "e2_nFit", "e5", "e5_nFit", "e8", "e8_nFit")) %>% 
    filter(df == "15min") %>% 
    mutate(model = factor(eq,
                          levels = c("e2", "e2_nFit", "e5", "e5_nFit", "e8", "e8_nFit"),
                          labels = c("2", "2-2", "5", "5-2", "8", "8-2"))) %>% 
    mutate(solute = factor(solute,
                           levels = c("doc", "no3"),
                           labels = c("DOC", "NO3"))) %>% 
    mutate(c0_source = factor(c0_source,
                              levels = c("baseflow", "wet_deposition"),
                              labels = c("Baseflow", "Wet dep."))) %>% 
    ggplot(aes(x = model, y = nrmse_mean, color = c0_source)) +
    facet_grid(rows = vars(site), cols = vars(solute), scales = "free") +
    geom_point(size = 2, position = position_dodge(0.8)) +
    geom_errorbar(aes(ymin = nrmse_ci_lower, ymax = nrmse_ci_upper), width = 0, position = position_dodge(0.8)) + 
    scale_color_manual(name = "C0 source",
                       values = c("#1f78b4", "#a6cee3")) +
    ylab("NRMSE") +
    xlab("Model") +
    theme_bw() +
    theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
          strip.background = element_blank())
  
  ggsave(here("plots", "fig_c0_comparison_cv.tiff"), height = 6, width = 5.5, units = "in", dpi = 150)    
  
  ## n source ----
  ## Compare C-V for models w/ n provided vs. estimated for models w/ param 'n' using C0 = baseflow and 15-min data
  cv_results %>% 
    filter(!is.na(eq)) %>% 
    # Remove outliers
    filter(nrmse_mean < 10) %>% 
    # Choose dataset
    filter(df == "15min") %>% 
    # Choose C0 source
    filter(c0_source == "baseflow") %>% 
    filter(eq %in% c("e1", "e1_nFit", "e2", "e2_nFit", "e4", "e4_nFit",
                     "e5", "e5_nFit", "e6", "e6_nFit", "e7", "e7_nFit", "e8", "e8_nFit")) %>%
    mutate(solute = factor(solute, levels = c("doc", "no3"), labels = c("DOC", "NO3"))) %>% 
    mutate(n_source = factor(n_source, levels = c("input", "fit"), labels = c("Input", "Fit"))) %>% 
    mutate(model = parse_number(eq)) %>% 
    ggplot(aes(x = factor(model), y = nrmse_mean, color = n_source)) +
    facet_grid(rows = vars(site), cols = vars(solute), scales = "free") +
    geom_point(size = 2, position = position_dodge(0.8)) +
    geom_errorbar(aes(ymin = nrmse_ci_lower, ymax = nrmse_ci_upper), width = 0, position = position_dodge(0.8)) + 
    scale_color_manual(name = "'n' source",
                       values = c("#264653", "#2a9d8f")) +
    ylab("NRMSE") +
    xlab("Model") +
    theme_bw() +
    theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
          strip.background = element_blank())
  
  ggsave(here("plots", "fig_n_source_comparison_cv.tiff"), height = 6, width = 5.5, units = "in", dpi = 150)
    
  ## all models ----
  ## Compare C-V for all models across data resolution
    # n_source = Input
    cv_results %>%
      filter(!is.na(eq)) %>% 
      # Remove outliers
      filter(nrmse_mean < 10) %>% 
      # Choose C0 source
      filter(c0_source == "baseflow") %>% 
      # Choose n source
      filter(n_source == "input") %>% 
      mutate(df = factor(df,
                         levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                         labels = c("G-\nWeekly", "S-\n15-min", "S-\nDaily",
                                    "S-\nWeekly", "S-\nMonthly"))) %>%
      mutate(solute = factor(solute, levels = c("doc", "no3"), labels = c("DOC", "NO3"))) %>% 
      mutate(n_source = factor(n_source, levels = c("input", "fit"), labels = c("Input", "Fit"))) %>%
      mutate(model = factor(parse_number(eq))) %>%       
      mutate(model = factor(model,
                            levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                            labels = c("Eq. 0", "Eq. 1", "Eq. 2", "Eq. 3", "Eq. 4", "Eq. 5", "Eq. 6", "Eq. 7", "Eq. 8"))) %>%
      ggplot(aes(x = df, y = nrmse_mean, color = model)) +
      facet_grid(rows = vars(site), cols = vars(solute), scales = "free") +
      geom_point(size = 2, position = position_dodge(0.8)) +
      geom_errorbar(aes(ymin = nrmse_ci_lower, ymax = nrmse_ci_upper), width = 0, position = position_dodge(0.8)) +
      scico::scale_colour_scico_d(palette = "batlow", name = "Model", drop = FALSE) +
      # scico::scale_fill_scico_d(palette = "batlow", name = "Model", drop = FALSE) +
      ylab("NRMSE") +
      xlab("Data set") +
      coord_cartesian(ylim=c(0, 1.25)) +
      theme_bw() +
      theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
            strip.background = element_blank(),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 16),
            axis.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))
    
    ggsave(here("plots", "fig_all_model_comparison_cv.tiff"), height = 7, width = 8, units = "in", dpi = 150)
    

    
# Plot residuals for 15-min model (n source = input; C0 = baseflow) ----
  # Need to figure out how to get backtransformed residuals for log-log model (probably just need to add C and predicted to output)
    
  # Bin residuals by flow percentile
  # Reminder about flow-duration percentiles: lower % values = higher flows; higher % values = lower flows
  # A 5-percent exceedance probability represents a high flow that has been exceeded only 5-percent of all 
  # days of the flow record. Conversely, a 95-percent exceedance probability would characterize low-flow 
  # conditions in a stream, because 95 percent of all daily mean flows in the record are greater than that amount.    
    
  ## DOC ----
  resids %>%
    filter(solute == "doc") %>% 
    # For now filter out log-log model
    filter(model != "eq0") %>% 
    filter(df == "15min") %>% 
    # Choose n source
    filter(n_source == "input") %>%
    # Tidy data
    mutate(solute = factor(solute, levels = c("doc", "no3"), labels = c("DOC", "NO3"))) %>% 
    mutate(n_source = factor(n_source, levels = c("input", "fit"), labels = c("Input", "Fit"))) %>%
    mutate(model = factor(parse_number(model))) %>%       
    mutate(model = factor(model,
                          levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                          labels = c("Eq. 1", "Eq. 2", "Eq. 3", "Eq. 4", "Eq. 5", "Eq. 6", "Eq. 7", "Eq. 8"))) %>%     
    group_by(site, model) %>% 
    # Rank daily discharges from largest to smallest value fr each site
    # Largest discharge value has a rank value of 1
    arrange(model, desc(q)) %>% 
    mutate(rank = row_number()) %>% 
    # Total the number of discharge records for each site
    mutate(n = n()) %>% 
    # Calculate exceedence probability (P)/flow-duration percentile:
    mutate(ex_prob = 100 * (rank / (n + 1))) %>% 
    # Divide into EPA flow percentile zones
    mutate(fp_zone_epa = ifelse(ex_prob <= 10, "High", NA),
           fp_zone_epa = ifelse(is.na(fp_zone_epa) & (ex_prob > 10 & ex_prob <= 40), "Moist", fp_zone_epa),
           fp_zone_epa = ifelse(is.na(fp_zone_epa) & (ex_prob > 40 & ex_prob <= 60), "Mid-range", fp_zone_epa),
           fp_zone_epa = ifelse(is.na(fp_zone_epa) & (ex_prob > 60 & ex_prob <= 90), "Dry", fp_zone_epa),
           fp_zone_epa = ifelse(is.na(fp_zone_epa) & ex_prob > 90, "Low", fp_zone_epa)) %>% 
    mutate(fp_zone_epa = factor(fp_zone_epa, levels = c("High", "Moist", "Mid-range", "Dry", "Low"))) %>% 
    ggplot(aes(x = fp_zone_epa, y = resids, color = model)) +
    facet_grid(rows = vars(site), cols = vars(model), scales = "free") +
    geom_boxplot(outliers = FALSE) +
    scale_colour_scico_d(palette = "batlow", name = "Model") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Flow conditions",
         y = "Residuals") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16, margin = margin(7, 0, 0, 0)),
          axis.title.y = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14))
    
  ggsave(here("plots", "fig_resids_doc.tiff"), height = 7.5, width = 10, units = "in", dpi = 150)

  ## NO3 ----
  resids %>%
    filter(solute == "no3") %>% 
    # For now filter out log-log model
    filter(model != "eq0") %>% 
    filter(df == "15min") %>% 
    # Choose n source
    filter(n_source == "input") %>%
    # Tidy data
    mutate(solute = factor(solute, levels = c("doc", "no3"), labels = c("DOC", "NO3"))) %>% 
    mutate(n_source = factor(n_source, levels = c("input", "fit"), labels = c("Input", "Fit"))) %>%
    mutate(model = factor(parse_number(model))) %>%       
    mutate(model = factor(model,
                          levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                          labels = c("Eq. 1", "Eq. 2", "Eq. 3", "Eq. 4", "Eq. 5", "Eq. 6", "Eq. 7", "Eq. 8"))) %>%     
    group_by(site, model) %>% 
    # Rank daily discharges from largest to smallest value fr each site
    # Largest discharge value has a rank value of 1
    arrange(model, desc(q)) %>% 
    mutate(rank = row_number()) %>% 
    # Total the number of discharge records for each site
    mutate(n = n()) %>% 
    # Calculate exceedence probability (P)/flow-duration percentile:
    mutate(ex_prob = 100 * (rank / (n + 1))) %>% 
    # Divide into EPA flow percentile zones
    mutate(fp_zone_epa = ifelse(ex_prob <= 10, "High", NA),
           fp_zone_epa = ifelse(is.na(fp_zone_epa) & (ex_prob > 10 & ex_prob <= 40), "Moist", fp_zone_epa),
           fp_zone_epa = ifelse(is.na(fp_zone_epa) & (ex_prob > 40 & ex_prob <= 60), "Mid-range", fp_zone_epa),
           fp_zone_epa = ifelse(is.na(fp_zone_epa) & (ex_prob > 60 & ex_prob <= 90), "Dry", fp_zone_epa),
           fp_zone_epa = ifelse(is.na(fp_zone_epa) & ex_prob > 90, "Low", fp_zone_epa)) %>% 
    mutate(fp_zone_epa = factor(fp_zone_epa, levels = c("High", "Moist", "Mid-range", "Dry", "Low"))) %>% 
    ggplot(aes(x = fp_zone_epa, y = resids, color = model)) +
    facet_grid(rows = vars(site), cols = vars(model), scales = "free") +
    geom_boxplot(outliers = FALSE) +
    scale_colour_scico_d(palette = "batlow", name = "Model") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Flow conditions",
         y = "Residuals") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16, margin = margin(7, 0, 0, 0)),
          axis.title.y = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14))
  
  ggsave(here("plots", "fig_resids_no3.tiff"), height = 7.5, width = 10, units = "in", dpi = 150)  
    
    
    
    
    
    
    
    

# OLD PLOTTING CODE ----

# Aggregate and plot results from above ----      
## Plot reduced sensor data & model fits ----
# Create a single df for cq data
df_all <-
  df_grab %>% 
  mutate(df = "lt_grab") %>% 
  bind_rows(df_15min %>% 
              mutate(df = "15min"),
            df_15min_daily_list[[1]] %>% 
              mutate(df = "daily"),
            df_15min_weekly_list[[1]] %>% 
              mutate(df = "weekly"),
            df_15min_monthly_list[[1]] %>% 
              mutate(df = "monthly"))
# Write to CSV for future use
df_all %>% 
  write_csv(here("data", "cached_model_data", paste0("cache_cq_all_datasets", "_", sol_choice, "_", site_choice, ".csv")))

# Read in past results:
# df_all <- read_csv(here("data", "cached_model_data", paste0("cache_cq_all_datasets", "_", sol_choice, "_", site_choice, ".csv")))

# Faceted plot of CQ data
p_cq_all <-
  df_all %>%
  mutate(df = factor(df, 
                     levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                     labels = c("Grab samples", "Sensor: 15-min", "Sensor: Daily", "Sensor: Weekly", "Sensor: Monthly"))) %>% 
  ggplot(aes(x = q, y = c)) +
  facet_wrap(~df, ncol = 1) +
  geom_point(shape = 1, alpha = 0.5) +
  # scale_y_log10(limits = c(-0.1, 1.5), breaks = c(0.0001, 0.0100, 1.0000)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", math_format(10^.x)),
                limits = c(10^-5, 10^2)) +    
  scale_x_log10() +
  labs(y = expression("NO"[3]~"-N (mg L"^-1~")"),
       x = expression(Discharge~(L~s^-1)),
       title = "Raw CQ data") +    
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(size = 13),
        plot.margin = margin(t = 5.5,  # Top margin
                             r = 9,  # Right margin
                             b = 5.5,  # Bottom margin
                             l = 5.5)) # Left margin)

# Create a single df for all bdf data
df_bdf_all <-
  bdf_grab %>% 
  mutate(df = "lt_grab") %>% 
  bind_rows(bdf_15min %>% 
              mutate(df = "15min"),
            bdf_daily %>% 
              mutate(df = "daily"),
            bdf_weekly %>% 
              mutate(df = "weekly"),
            bdf_monthly %>% 
              mutate(df = "monthly"))

# Write to CSV for future use
df_bdf_all %>% 
  write_csv(here("data", "cached_model_data",  paste0("cache_cq_all_bdf", "_", sol_choice, "_", site_choice, ".csv")))


# Create a single df for plot_df_X data
df_fit_all <-
  plot_df_grab %>% 
  mutate(df = "lt_grab") %>% 
  bind_rows(plot_df_15min %>% 
              mutate(df = "15min"),
            plot_df_daily %>% 
              mutate(df = "daily"),
            plot_df_weekly %>% 
              mutate(df = "weekly"),
            plot_df_monthly %>% 
              mutate(df = "monthly"))

# Write to CSV for future use
df_fit_all %>% 
  write_csv(here("data", "cached_model_data",  paste0("cache_cq_all_fits", "_", sol_choice, "_", site_choice, ".csv")))

# Read in past results:
# df_bdf_all <- read_csv(here("data", "cached_model_data",  paste0("cache_cq_all_bdf", "_", sol_choice, "_", site_choice, ".csv")))
# df_fit_all <- read_csv(here("data", "cached_model_data",  paste0("cache_cq_all_fits", "_", sol_choice, "_", site_choice, ".csv")))

# Faceted plot of fits
p_fit_all <-
  df_fit_all %>% 
  mutate(df = factor(df, 
                     levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                     labels = c("Grab samples: Weekly", "Sensor: 15-min", "Sensor: Daily",
                                "Sensor: Weekly", "Sensor: Monthly"))) %>% 
  mutate(model = factor(model,
                        levels = c("eq1", "eq2", "eq3", "eq4", "eq5", "eq6", "eq7", "eq8"),
                        labels = c("Eq. 1", "Eq. 2", "Eq. 3", "Eq. 4", "Eq. 5", "Eq. 6", "Eq. 7", "Eq. 8"))) %>% 
  filter(value > 0) %>%
  ggplot(aes(x = q, y = value)) +
  geom_point(data = df_bdf_all %>% 
               mutate(df = factor(df, 
                                  levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                                  labels = c("Grab samples: Weekly", "Sensor: 15-min", "Sensor: Daily",
                                             "Sensor: Weekly", "Sensor: Monthly"))), 
             mapping = aes(bqs, bcs), shape = 1, size = 1.7, color = "gray30") +
  geom_line(aes(colour = model, linetype =  model), linewidth = 0.8) +
  facet_wrap(~df, ncol = 1) +
  scale_colour_scico_d(palette = "batlow", name = "Model") +
  scale_linetype_discrete(name = "Model") +
  geom_hline(yintercept = c0_bf, linetype = "dashed", color = "gray25") +
  geom_hline(yintercept = c0_wd, linetype = "dashed", color = "lightskyblue1") +      
  scale_y_log10() +
  # scale_y_log10(limits = c(-0.1, 1.5), breaks = c(0.0001, 0.0100, 1.0000)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", math_format(10^.x)),
                limits = c(10^-5, 10^2)) +
  scale_x_log10() +
  # scale_x_log10(limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.2)) +
  labs(x = expression(Discharge~(L~s^-1)),
       y = expression("NO"[3]~"-N (mg L"^-1~")"),
       title = "Model fits",
       subtitle = "Complete Q range") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        # axis.text.y = element_blank(),
        plot.margin = margin(t = 5.5,  # Top margin
                             r = 7,  # Right margin
                             b = 5.5,  # Bottom margin
                             l = 5.5)) # Left margin

# Try making zoomed in plots of "whiskered ends" for first two datasets
# Low end
p_low <-
  df_fit_all %>% 
  mutate(df = factor(df, 
                     levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                     labels = c("Grab samples: Weekly", "Sensor: 15-min", "Sensor: Daily",
                                "Sensor: Weekly", "Sensor: Monthly"))) %>% 
  mutate(model = factor(model,
                        levels = c("eq1", "eq2", "eq3", "eq4", "eq5", "eq6", "eq7", "eq8"),
                        labels = c("Eq. 1", "Eq. 2", "Eq. 3", "Eq. 4", "Eq. 5", "Eq. 6", "Eq. 7", "Eq. 8"))) %>% 
  ggplot(aes(x = q, y = value)) +
  geom_line(aes(colour = model, linetype =  model), size = 0.8) +
  geom_point(data = df_bdf_all %>% 
               mutate(df = factor(df, 
                                  levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                                  labels = c("Grab samples: Weekly", "Sensor: 15-min", "Sensor: Daily",
                                             "Sensor: Weekly", "Sensor: Monthly"))), 
             mapping = aes(bqs, bcs), shape = 1, size = 1.7, color = "gray10") +
  facet_wrap(~df, ncol = 1) +
  scale_colour_scico_d(palette = "batlow", name = "Model") +
  scale_linetype_discrete(name = "Model") +
  # scale_y_log10() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", math_format(10^.x)),
                limits = c(10^-2, 10^2)) +
  # scale_y_log10(limits = c(175, 250), breaks = c(200, 250)) +
  # scale_x_log10() +
  scale_x_log10(limits = c(0.15, 6), breaks = c(0.3, 3)) +
  # scale_x_log10(limits = c(0.1, 6)) +
  labs(x = expression(Discharge~(L~s^-1)),
       y = expression("NO"[3]~"-N (mg L"^-1~")"),
       subtitle = "Low Q") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank())

# High end
p_high <-
  df_fit_all %>% 
  mutate(df = factor(df, 
                     levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                     labels = c("Grab samples", "Sensor: 15-min", "Sensor: Daily",
                                "Sensor: Weekly", "Sensor: Monthly"))) %>% 
  mutate(model = factor(model,
                        levels = c("eq1", "eq2", "eq3", "eq4", "eq5", "eq6", "eq7", "eq8"),
                        labels = c("Eq. 1", "Eq. 2", "Eq. 3", "Eq. 4", "Eq. 5", "Eq. 6", "Eq. 7", "Eq. 8"))) %>% 
  ggplot(aes(x = q, y = value)) +
  geom_line(aes(colour = model, linetype =  model), size = 0.8) +
  geom_point(data = df_bdf_all %>% 
               mutate(df = factor(df, 
                                  levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                                  labels = c("Grab samples", "Sensor: 15-min", "Sensor: Daily",
                                             "Sensor: Weekly", "Sensor: Monthly"))), 
             mapping = aes(bqs, bcs), shape = 1, size = 1.7, color = "gray10") +
  facet_wrap(~df, ncol = 1, strip.position = "right") +
  scale_colour_scico_d(palette = "batlow", name = "Model") +
  scale_linetype_discrete(name = "Model") +
  # scale_y_log10() +
  # scale_x_log10() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", math_format(10^.x)),
                limits = c(10^-6, 10^1)) +
  # scale_x_log10(limits = c(25, 150)) +
  scale_x_log10(limits = c(25, 150), breaks = c(30, 100)) +
  labs(x = expression(Discharge~(L~s^-1)),
       y = expression("NO"[3]~"-N (mg L"^-1~")"),
       subtitle = "High Q") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_blank())       


# Combine plots with patchwork
xlab <- p_cq_all$labels$x
p_cq_all$labels$x <- p_fit_all$labels$x <- " "

tiff(here("plots", paste0("fig_model_fits", "_", sol_choice, "_", site_choice, ".tiff")), width = 10, height = 8.5, units = "in", res = 600)
p_cq_all + p_fit_all + p_low + p_high + 
  plot_layout(ncol = 4, widths = c(0.7, 1, 0.4, 0.4))
grid::grid.draw(grid::textGrob(xlab, y = 0.02, gp = gpar(fontsize = 13)))
dev.off()


## Plot cross-validation results ----
# Create a single df for out2_X (C-V) data
df_cv_all <-
  out2_grab %>% 
  mutate(df = "lt_grab") %>% 
  bind_rows(out2_15min %>% 
              mutate(df = "15min"),
            out2_daily %>% 
              mutate(df = "daily"),
            out2_weekly %>% 
              mutate(df = "weekly"),
            out2_monthly %>% 
              mutate(df = "monthly"))

# Write to CSV for future use
df_cv_all %>% 
  write_csv(here("data", "cached_model_data",  paste0("cache_cq_all_cv", "_", sol_choice, "_", site_choice, ".csv")))

# Read in past results
df_cv_all <-
  read_csv(here("data", "cached_model_data",  paste0("cache_cq_all_cv", "_", sol_choice, "_", site_choice, ".csv")))

# Plot all at once
# df_cv_all %>%
#   ggplot(aes(x = eq, y = rmse_mean, fill = df, color = df)) +
#   geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +    
#   geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, position = position_dodge(0.9)) +
#   scale_colour_scico_d(palette = "batlow", name = "Model") +
#   scale_fill_scico_d(palette = "batlow", name = "Model") +
#   ylab("RMSE") +
#   xlab("Equation") +
#   theme_classic() +
#   theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))

med_c <-
  df_15min %>% 
  summarize(med_c = median(c, na.rm = TRUE)) %>% 
  pull(med_c)    

df_cv_all %>%
  mutate(df = factor(df,
                     levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                     labels = c("G-\nWeekly", "S-\n15-min", "S-\nDaily",
                                "S-\nWeekly", "S-\nMonthly"))) %>%
  mutate(eq = factor(eq,
                     levels = c("e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8"),
                     labels = c("Eq. 1", "Eq. 2", "Eq. 3", "Eq. 4", "Eq. 5", "Eq. 6", "Eq. 7", "Eq. 8"))) %>%
  ggplot(aes(x = df, y = rmse_mean, color = eq)) +
  geom_point(size = 2, position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, position = position_dodge(0.8)) +
  geom_hline(yintercept = med_c, linetype = "dashed", color = "gray50") +
  scico::scale_colour_scico_d(palette = "batlow", name = "Model") +
  scale_fill_scico_d(palette = "batlow", name = "Model") +
  ylab("RMSE") +
  xlab("Data set") +
  theme_classic() +
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))

ggsave(here("plots", paste0("fig_model_cv", "_", sol_choice, "_", site_choice, ".tiff")), height = 3, width = 5, units = "in", dpi = 600)


# On average, how much greater was the RMSE for eq 6 vs the lowest RMSE for that dataset
df_cv_all %>% 
  select(eq, rmse_mean, df) %>% 
  pivot_wider(names_from = eq, values_from = rmse_mean) %>% 
  rowwise() %>% 
  mutate(min_RMSE = min(e1, e2, e3, e4, e5, e6, e7, e8)) %>% 
  mutate(per_inc_e6 = (e6 - min_RMSE)/min_RMSE*100)

# On average, how much greater were the RMSEs for the the weekly grab samples vs the full sensor dataset?
df_cv_all %>% 
  group_by(df) %>% 
  summarize(RMSE_mean = mean(rmse_mean)) %>% 
  pivot_wider(names_from = df, values_from = RMSE_mean) %>% 
  mutate(per_inc_gr_weekly = (lt_grab - `15min`)/`15min`*100)

# Table of CV results ----
df_cv_all %>% 
  mutate(Source = ifelse(df == "lt_grab", "Grab samples", "Sensor")) %>% 
  mutate(Frequency = ifelse(df == "lt_grab", "Weekly", NA),
         Frequency = ifelse(is.na(Frequency) & df == "15min", "15 min", Frequency),
         Frequency = ifelse(is.na(Frequency) & df == "daily", "Daily", Frequency),
         Frequency = ifelse(is.na(Frequency) & df == "weekly", "Weekly", Frequency),
         Frequency = ifelse(is.na(Frequency) & df == "monthly", "Monthly", Frequency)) %>% 
  mutate(Eq. = parse_number(eq)) %>% 
  select(Source, Frequency, `Eq.`, rmse_mean, ci_lower, ci_upper) %>% 
  mutate(rmse_mean = format(round(rmse_mean, 2), nsmall = 2),
         ci_lower = format(round(ci_lower, 2), nsmall = 2),
         ci_upper = format(round(ci_upper, 2), nsmall = 2)) %>% 
  rename(RMSE = rmse_mean, `95% CI, lower` = ci_lower, `95% CI, upper` = ci_upper) %>% 
  write_csv(here("tables", paste0("table_cv_results_summary", "_", sol_choice, "_", site_choice, ".csv")))



# Plot residuals ----
# Low Flows (0–10th percentile),
# Dry Conditions (10–40th percentile), 
# Mid-Range Flows (40–60th percentile), 
# Moist Conditions (60–90th percentile)
# High Flows (90–100th percentile; USEPA 2007).
# US EPA (2007) An approach using load duration curves in the development of TMDLs. EPA 841-B-07–006.
# Arial suggests the package hydrotsm will calc flow percentiles
# Shannon also has code she will share

# Create a single df for resids_X data
df_resids_all <-
  resids_grab %>% 
  mutate(df = "lt_grab") %>% 
  bind_rows(resids_15min %>% 
              mutate(df = "15min"),
            resids_daily %>% 
              mutate(df = "daily"),
            resids_weekly %>% 
              mutate(df = "weekly"),
            resids_monthly %>% 
              mutate(df = "monthly"))
# Write to CSV for future use
df_resids_all %>% 
  write_csv(here("data", "cached_model_data",  paste0("cache_cq_all_resids", "_", sol_choice, "_", site_choice, ".csv")))

# Read past results
# df_resids_all <-
#   read_csv(here("data", "cached_model_data",  paste0("cache_cq_all_resids", "_", sol_choice, "_", site_choice, ".csv")))


# Plot just the 15-min sensor data: bin by flow percentile
# Reminder about flow-duration percentiles: lower % values = higher flows; higher % values = lower flows
# A 5-percent exceedance probability represents a high flow that has been exceeded only 5-percent of all 
# days of the flow record. Conversely, a 95-percent exceedance probability would characterize low-flow 
# conditions in a stream, because 95 percent of all daily mean flows in the record are greater than that amount. 

df_resids_all %>% 
  filter(df == "15min") %>% 
  group_by(model) %>% 
  # Rank daily discharges from largest to smallest value fr each site
  # Largest discharge value has a rank value of 1
  arrange(model, desc(q)) %>% 
  mutate(rank = row_number()) %>% 
  # Total the number of discharge records for each site
  mutate(n = n()) %>% 
  # Calculate exceedence probability (P)/flow-duration percentile:
  mutate(ex_prob = 100 * (rank / (n + 1))) %>% 
  # Divide into EPA flow percentile zones
  mutate(fp_zone_epa = ifelse(ex_prob <= 10, "High", NA),
         fp_zone_epa = ifelse(is.na(fp_zone_epa) & (ex_prob > 10 & ex_prob <= 40), "Moist", fp_zone_epa),
         fp_zone_epa = ifelse(is.na(fp_zone_epa) & (ex_prob > 40 & ex_prob <= 60), "Mid-range", fp_zone_epa),
         fp_zone_epa = ifelse(is.na(fp_zone_epa) & (ex_prob > 60 & ex_prob <= 90), "Dry", fp_zone_epa),
         fp_zone_epa = ifelse(is.na(fp_zone_epa) & ex_prob > 90, "Low", fp_zone_epa)) %>% 
  mutate(fp_zone_epa = factor(fp_zone_epa, levels = c("High", "Moist", "Mid-range", "Dry", "Low"))) %>% 
  # Create bins by flow percentile
  # mutate(ex_prob_bin = cut(ex_prob, breaks = seq(0, 100, by = 10))) %>% 
  mutate(model = factor(model,
                        levels = c("eq1", "eq2", "eq3", "eq4", "eq5", "eq6", "eq7", "eq8"),
                        labels = c("Eq. 1", "Eq. 2", "Eq. 3", "Eq. 4", "Eq. 5", "Eq. 6", "Eq. 7", "Eq. 8"))) %>%
  ggplot(aes(x = fp_zone_epa, y = resids, color = model)) +
  facet_wrap(~model, ncol = 2) +
  geom_boxplot() +
  scale_colour_scico_d(palette = "batlow", name = "Model") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Flow conditions",
       y = "Residuals") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(margin = margin(7, 0, 0, 0)),
        axis.text.x = element_text(size = 7))

ggsave(here("plots", paste0("fig_model_resids_epa_flow_zones", "_", sol_choice, "_", site_choice, ".tiff")), height = 7, width = 4.5, units = "in", dpi = 300)



# CHECK DISTRIBUTIONS OF C AND Q FOR ALL DATASETS USED ----

# Read in past results  
# df_all <- 
#   read_csv(here("data", "cached_model_data", paste0("cache_cq_all_datasets", "_", sol_choice, "_", site_choice, ".csv")))  

# Plot c and q separately so that I can keep x-axis ranges same
# Calculate medians to plot on histograms
meds_cq <-
  df_all %>% 
  mutate(df = factor(df,
                     levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                     labels = c("G-Weekly", "S-15min", "S-Daily", "S-Weekly", "S-Monthly"))) %>%
  group_by(df) %>% 
  summarize(med_c = median(c),
            q1_c = quantile(c, 1/4),
            q3_c = quantile(c, 3/4),
            med_q = median(q),
            q1_q = quantile(q, 1/4),
            q3_q = quantile(q, 3/4))

# Conc
p_c_dist <-
  df_all %>% 
  mutate(df = factor(df,
                     levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                     labels = c("G-Weekly", "S-15min", "S-Daily", "S-Weekly", "S-Monthly"))) %>% 
  ggplot(aes(x = c)) +
  facet_wrap(~df, ncol = 1, scales = "free_y") +
  geom_histogram() +
  geom_vline(data = meds_cq, aes(xintercept = med_c)) +
  geom_vline(data = meds_cq, aes(xintercept = q1_c), linetype = "dashed") +
  geom_vline(data = meds_cq, aes(xintercept = q3_c), linetype = "dashed") +
  labs(x = expression("NO"[3]~"-N (mg L"^-1~")"),
       y = "Count") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

# Discharge
p_q_dist <-
  df_all %>% 
  mutate(df = factor(df,
                     levels = c("lt_grab", "15min", "daily", "weekly", "monthly"),
                     labels = c("G-Weekly", "S-15min", "S-Daily", "S-Weekly", "S-Monthly"))) %>% 
  ggplot(aes(x = q)) +
  facet_wrap(~df, ncol = 1, scales = "free_y") +
  geom_vline(data = meds_cq, aes(xintercept = med_q)) +
  geom_vline(data = meds_cq, aes(xintercept = q1_q), linetype = "dashed") +
  geom_vline(data = meds_cq, aes(xintercept = q3_q), linetype = "dashed") +      
  geom_histogram() +
  labs(x = expression(Discharge~(L~s^-1)),
       y = "Count") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank())

# Combine plots
p_c_dist + p_q_dist + plot_annotation(tag_levels = "a", tag_suffix = ")")

ggsave(here("plots", paste0("fig_compare_cq_histograms", "_", sol_choice, "_", site_choice, ".tiff")), width = 7.5, height = 6, units = "in", dpi = 600)    
