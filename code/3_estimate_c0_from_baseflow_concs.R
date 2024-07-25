# Calculate C0 from baseflow concentrations


# Load libraries ----
  library("tidyverse")
  library("data.table")
  library("plotly")

# Options
  options(scipen=999)

# Load filtered data (from 1_prep_cq_data.R)
  dat_all <- 
    data.table::fread("data/nh_cq_data_filtered.csv")


# Estimate C0 as concentration at 99% flow (lowest flows) ----
  # Note: below, I summarized C0 based on flows less than the inflection point 
  # at the high end of the exceedence probability curve (~99%); results weren't
  # different enough than just filtering at 99%, so I'm just going to use that
  
  ## Calculate flow percentiles
  q_perc <- 
    # Use 15-min sensor data
    dat_all %>% 
    pivot_longer(cols = c(spc, no3, fdom, doc), names_to = "var", values_to = "conc") %>%
    filter(!is.na(conc)) %>%
    filter(!is.na(q)) %>%
    group_by(site, var) %>%
    # Rank daily discharges from largest to smallest value fr each site
    # Largest discharge value has a rank value of 1
    arrange(desc(q)) %>% 
    mutate(rank = row_number()) %>% 
    # Total the number of discharge records for each site
    mutate(n = n()) %>% 
    # Calculate exceedence probability (P)/flow-duration percentile:
    mutate(ex_prob = 100 * (rank / (n + 1))) %>% 
    # Create bins by flow percentile
    mutate(ex_prob_bin = cut(ex_prob, breaks = seq(0, 100, by = 10))) %>%
    select(site, datetime, var, conc, q, ex_prob, ex_prob_bin) %>%
    ungroup()
  
  ## Summarize concentration at lowest flows (ex_prob >= 99%)
    # Visualize the concentration range at baseflow
      # All sites and vars
      q_perc %>% 
        group_by(site, var) %>% 
        # Filter for lowest flows (99% flow probability)
        filter(ex_prob >= 99) %>% 
        ggplot(aes(x = site, y = conc)) +
        facet_wrap(~var, scales = "free_y") +
        geom_boxplot() +
        theme_bw()
      
      # Just NO3
      q_perc %>% 
        group_by(site, var) %>% 
        # Filter for lowest flows (99% flow probability)
        filter(ex_prob >= 99) %>% 
        filter(var == "no3") %>% 
        ggplot(aes(x = site, y = conc)) +
        facet_wrap(~site, scales = "free") +
        geom_boxplot() +
        theme_bw()      
  
  C0_baseflow <-
    q_perc %>% 
    group_by(site, var) %>% 
    # Filter for lowest flows (99% flow probability)
    filter(ex_prob >= 99) %>% 
    # Summarize concentration
    summarize(med_c = median(conc))
  # mutate(across(where(is.numeric), ~ signif(.x, digits = 3)))
  
  # Write to CSV for discussion
  C0_baseflow %>% 
    write_csv(file = "data/c0_baseflow_estimates.csv")  
  
  
  
  
# # ALTERNATIVE METHOD: Estimate C0 as concentration at inflection point of the high end
# # of the exceedence probability cuve (~99%) ----
#   ## Calculate flow percentiles
#   q_perc <- 
#     # Use 15-min sensor data
#     dat_all %>% 
#     # pivot_longer(cols = c(spc, no3, fdom), names_to = "var", values_to = "conc") %>% 
#     # filter(!is.na(conc)) %>% 
#     filter(!is.na(q)) %>%
#     # group_by(site, var) %>% 
#     group_by(site) %>% 
#     # Rank daily discharges from largest to smallest value fr each site
#     # Largest discharge value has a rank value of 1
#     arrange(desc(q)) %>% 
#     mutate(rank = row_number()) %>% 
#     # Total the number of discharge records for each site
#     mutate(n = n()) %>% 
#     # Calculate exceedence probability (P)/flow-duration percentile:
#     mutate(ex_prob = 100 * (rank / (n + 1))) %>% 
#     # Create bins by flow percentile
#     mutate(ex_prob_bin = cut(ex_prob, breaks = seq(0, 100, by = 10))) %>%
#     # select(site, datetime, var, conc, q, ex_prob, ex_prob_bin) %>% 
#     select(site, datetime, q, ex_prob, ex_prob_bin) %>% 
#     ungroup()
#   
#   ## Plot the flow exceedence probability curves to find breakpoint for "low flows"
#   unique(q_perc$site)
#   # "BEF" "LMP" "DCF" "WHB" "HBF" "BDC"
# 
#   pl_ex_prob <-
#     q_perc %>% 
#     filter(site == "BEF") %>% 
#     filter(ex_prob > 80) %>% 
#     ggplot(aes(x = ex_prob, y = q)) +
#     facet_wrap(~site, scales = "free_y", ncol = 1) +
#     geom_point(shape = 1) +
#     scale_x_reverse() +
#     # scale_y_log10() +
#     theme_bw()
#   
#   ggplotly(pl_ex_prob)
#   
#   # Site specific inflection points for low flows
#   low_flows <- tibble(site = c("BEF", "LMP", "DCF", "WHB", "HBF", "BDC"),
#                       ex_prob_lf = c(99.8, 99.0, 99.6, 99.9, 99, 99.2))
#   
#   ## Find flow associated with exceedence probability inflection point
#   # test <-
#   q_perc %>% 
#     mutate(low_flows = ifelse(site == "BEF" & ex_prob > 99.8, "low",
#                               ifelse(site == "LMP" & ex_prob > 99.0, "low",
#                                      ifelse(site == "DCF" & ex_prob > 99.6, "low",
#                                             ifelse(site == "WHB" & ex_prob > 99.9, "low",
#                                                    ifelse(site == "HBF" & ex_prob > 99.0, "low",
#                                                           ifelse(site == "BDC" & ex_prob > 99.0, "low", "not_low"))))))) %>% 
#     filter(low_flows == "low") %>% 
#     arrange(site, desc(q)) %>% 
#     group_by(site) %>% 
#     slice(1) %>% 
#     select(site, q)
#   
#   ## Summarize concentration at lowest flows (ex_prob >= 99%)  
#   C0_baseflow <-
#     dat_all %>% 
#     filter(case_when(site == "BEF" ~ q <= 11.6,
#                      site == "LMP" ~ q <= 119.0,
#                      site == "DCF" ~ q <= 5.51,
#                      site == "WHB" ~ q <= 1.52,
#                      site == "HBF" ~ q <= 0.08,
#                      site == "BDC" ~ q <= 0.777)) %>% 
#     pivot_longer(cols = c(spc, no3, fdom), names_to = "var", values_to = "conc") %>% 
#     group_by(site, var) %>% 
#     summarize(med_c = median(conc, na.rm = TRUE))