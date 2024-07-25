# Create figures of log10(C) ~ log10(Q) comparing grabs vs sensor data

# Load packages ----
  library("tidyverse")
  library("data.table")
  library("patchwork")

# Options ----
  options(scipen = 999) #disable printing of numbers in scientific notation  

# Read in data ----
  ## CQ data from 1_prep_cq_data.R ----
  dat_all <- 
    data.table::fread("data/nh_cq_data_filtered.csv")

# Build plot functions
  plot_logC_vs_logQ <- function(df, site, var, var_grab, file_nm_var, min_c, max_c, axis_c_divisor, min_q, max_q, axis_q_divisor, y_lab) {
    # Filter for site
    df_sub <-
      df %>% 
      filter(site == {{ site }})
    
    # Filter for 15-min sensor data
    sensor <-
      df_sub %>% 
      filter(!is.na({{ var }})) %>%
      filter(!is.na(q)) %>% 
      rename(c = {{ var }})
    
    # Filter for long-term grab sample data
    grab <-
      df_sub %>% 
      filter(!is.na({{ var_grab }})) %>%
      filter(!is.na(q)) %>% 
      rename(c = {{ var_grab }})
    
    # Long-term grab data
      # Regress log10(C)~log10(Q)
      if(nrow(grab) >= 5){
        lm_grab <- lm(log10(grab$c) ~ log10(grab$q))
        b_grab <- round(coef(lm_grab)[2], 2)
        b_grab_ci05 <- round(confint(lm_grab, 'log10(grab$q)', level = 0.95)[1], 2)
        b_grab_ci95 <- round(confint(lm_grab, 'log10(grab$q)', level = 0.95)[2], 2)
      }

      # Plot
      if(nrow(grab) < 5){
        p_cq_grab <-
          grab %>% 
          ggplot(aes(x = log10(q), y = log10(c))) +
          geom_point(shape = 1, alpha = 0.2, color = "gray30") +
          # geom_smooth(method = "lm", se = FALSE) +
          scale_y_continuous(limits = c({{ min_c }}, {{ max_c }}), breaks = seq({{ min_c }}, {{ max_c }}, by = ({{ max_c }}-{{ min_c }}) / axis_c_divisor)) +
          scale_x_continuous(limits = c({{ min_q }}, {{ max_q }}), breaks = seq({{ min_q }}, {{ max_q }}, by = ({{ max_q }}-{{ min_q }}) / axis_q_divisor)) +
          theme_classic() +
          labs(subtitle = "Grabs",
               y = y_lab,
               x = expression(log[10](discharge))) +
          theme(axis.title.x = element_blank())        
      } else {
        p_cq_grab <-
          grab %>% 
          ggplot(aes(x = log10(q), y = log10(c))) +
          geom_point(shape = 1, alpha = 0.2, color = "gray30") +
          geom_smooth(method = "lm", se = FALSE) +
          scale_y_continuous(limits = c({{ min_c }}, {{ max_c }}), breaks = seq({{ min_c }}, {{ max_c }}, by = ({{ max_c }}-{{ min_c }}) / axis_c_divisor)) +
          scale_x_continuous(limits = c({{ min_q }}, {{ max_q }}), breaks = seq({{ min_q }}, {{ max_q }}, by = ({{ max_q }}-{{ min_q }}) / axis_q_divisor)) +
          theme_classic() +
          labs(subtitle = bquote(Grabs~(italic(b)==.(b_grab)~'['*.(b_grab_ci05)*','~.(b_grab_ci95)*']')),
               y = y_lab,
               x = expression(log[10](discharge))) +
          theme(axis.title.x = element_blank())        
      }
      

  
    # 15-min sensor data
      # Regress log10(C)~log10(Q)
      lm_sensor <- lm(log10(sensor$c) ~ log10(sensor$q))
      b_sensor <- round(coef(lm_sensor)[2], 2)
      b_sensor_ci05 <- round(confint(lm_sensor, 'log10(sensor$q)', level = 0.95)[1], 2)
      b_sensor_ci95 <- round(confint(lm_sensor, 'log10(sensor$q)', level = 0.95)[2], 2)      
    
      # Plot
      p_cq_15min <-
        sensor %>% 
        ggplot(aes(x = log10(q), y = log10(c))) +
        geom_point(shape = 1, alpha = 0.2, color = "gray30") +
        geom_smooth(method = "lm", se = FALSE) +
        scale_y_continuous(limits = c({{ min_c }}, {{ max_c }}), breaks = seq({{ min_c }}, {{ max_c }}, by = ({{ max_c }}-{{ min_c }}) / axis_c_divisor)) +
        scale_x_continuous(limits = c({{ min_q }}, {{ max_q }}), breaks = seq({{ min_q }}, {{ max_q }}, by = ({{ max_q }}-{{ min_q }}) / axis_q_divisor)) +
        theme_classic() +
        labs(subtitle = bquote(Sensor~(italic(b)==.(b_sensor)~'['*.(b_sensor_ci05)*','~.(b_sensor_ci95)*']')),
             y = y_lab,
             x = expression(log[10](discharge))) +
        theme(axis.title.y = element_blank())
      
    # Join plots together using patchwork
    ylab <- p_cq_grab$label$y
    p_cq_grab$labels$y <- p_cq_15min$labels$y <- " "

    # TIFF version
    tiff(paste0("plots/fig_CQ_logC_vs_logQ", "_", file_nm_var, "_", site, ".tiff"), height = 5, width = 3.5, units = "in", res = 600)
    pl_comb <- p_cq_grab / p_cq_15min + plot_annotation(tag_levels = "a", tag_suffix = ")") &
      # theme(plot.tag = element_text(hjust = -2))
      theme(plot.tag.position = c(0.05, 1),
            plot.margin = margin(t = 5,  # Top margin
                                 r = 5,  # Right margin
                                 b = 2,  # Bottom margin
                                 l = 10)) # Left margin
    plot(pl_comb)
    grid::grid.draw(grid::textGrob(ylab, x = 0.03, rot = 90))
    dev.off()

  }    
  
# For deciding on axis scales below  
  # Summarize max and mins of c and q
  dat_summ <-
    dat_all %>% 
    pivot_longer(cols = c(q, spc, no3, fdom, doc, spc_grab, no3_grab, fdom_grab, doc_grab), names_to = "var", values_to = "val") %>% 
    group_by(site, var) %>% 
    summarize(min = min(log10(val), na.rm = TRUE),
              max = max(log10(val), na.rm = TRUE))
  
  dat_summ %>% 
    filter(var %in% c("q", "spc", "spc_grab"))  

  
# Create plots ----
  unique(dat_all$site)
  # "BDC" "BEF" "DCF" "HBF" "LMP" "WHB"
  
  ## NO3 ----
    # BDC
    plot_logC_vs_logQ(df = dat_all, site = "BDC", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                      min_c = -3.5, max_c = 1, axis_c_divisor = 6, min_q = -1, max_q = 2.5, axis_q_divisor = 5, 
                      y_lab = expression(log[10]*"(NO"[3]~"-N)"))
    # BEF
    plot_logC_vs_logQ(df = dat_all, site = "BEF", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                      min_c = -6, max_c = 0, axis_c_divisor = 4, min_q = 0.5, max_q = 4.5, axis_q_divisor = 4,
                      y_lab = expression(log[10]*"(NO"[3]~"-N)"))
    # DCF
    plot_logC_vs_logQ(df = dat_all, site = "DCF", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                      min_c = -5.5, max_c = 0, axis_c_divisor = 5, min_q = 0.5, max_q = 4.5, axis_q_divisor = 4,
                      y_lab = expression(log[10]*"(NO"[3]~"-N)"))
    # HBF
    plot_logC_vs_logQ(df = dat_all, site = "HBF", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                      min_c = -3.5, max_c = 0, axis_c_divisor = 5, min_q = -1.5, max_q = 3.5, axis_q_divisor = 5,
                      y_lab = expression(log[10]*"(NO"[3]~"-N)"))
    # LMP
    plot_logC_vs_logQ(df = dat_all, site = "LMP", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                      min_c = -5.5, max_c = 0, axis_c_divisor = 5, min_q = 1.5, max_q = 5, axis_q_divisor = 5,
                      y_lab = expression(log[10]*"(NO"[3]~"-N)"))
    # WHB
    plot_logC_vs_logQ(df = dat_all, site = "WHB", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                      min_c = -3.5, max_c = 0, axis_c_divisor = 5, min_q = 0, max_q = 3.5, axis_q_divisor = 5,
                      y_lab = expression(log[10]*"(NO"[3]~"-N)"))
    
  ## DOC ----
    # BDC
    plot_logC_vs_logQ(df = dat_all, site = "BDC", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                      min_c = -2, max_c = 2, axis_c_divisor = 4, min_q = -1, max_q = 2.5, axis_q_divisor = 5,
                      y_lab = expression(log[10]*"(DOC)"))
    # BEF
    plot_logC_vs_logQ(df = dat_all, site = "BEF", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                      min_c = -0.5, max_c = 1, axis_c_divisor = 3, min_q = 0.5, max_q = 4.5, axis_q_divisor = 4, 
                      y_lab = expression(log[10]*"(DOC)"))
    # DCF
    plot_logC_vs_logQ(df = dat_all, site = "DCF", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                      min_c = 0, max_c = 1.5, axis_c_divisor = 3, min_q = 0.5, max_q = 4.5, axis_q_divisor = 4,
                      y_lab = expression(log[10]*"(DOC)"))
    # HBF
    plot_logC_vs_logQ(df = dat_all, site = "HBF", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                      min_c = 0, max_c = 1, axis_c_divisor = 4, min_q = -1.5, max_q = 3.5, axis_q_divisor = 5, 
                      y_lab = expression(log[10]*"(DOC)"))
    # LMP
    plot_logC_vs_logQ(df = dat_all, site = "LMP", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                      min_c = 0, max_c = 1.5, axis_c_divisor = 3, min_q = 1.5, max_q = 5, axis_q_divisor = 5, 
                      y_lab = expression(log[10]*"(DOC)"))
    # WHB
    plot_logC_vs_logQ(df = dat_all, site = "WHB", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                      min_c = -3.5, max_c = 1, axis_c_divisor = 6, min_q = 0, max_q = 3.5, axis_q_divisor = 5,
                      y_lab = expression(log[10]*"(DOC)"))
    
  ## DOC ----
    # BDC
    plot_logC_vs_logQ(df = dat_all, site = "BDC", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                      min_c = 1.5, max_c = 3, axis_c_divisor = 3, min_q = -1, max_q = 2.5, axis_q_divisor = 5,
                      y_lab = expression(log[10]*"(specific conductance)"))
    # BEF
    plot_logC_vs_logQ(df = dat_all, site = "BEF", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                      min_c = 0.5, max_c = 2, axis_c_divisor = 3, min_q = 0.5, max_q = 4.5, axis_q_divisor = 4, 
                      y_lab = expression(log[10]*"(specific conductance)"))
    # DCF
    plot_logC_vs_logQ(df = dat_all, site = "DCF", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                      min_c = 1, max_c = 3, axis_c_divisor = 4, min_q = 0.5, max_q = 4.5, axis_q_divisor = 4, 
                      y_lab = expression(log[10]*"(specific conductance)"))
    # HBF
    plot_logC_vs_logQ(df = dat_all, site = "HBF", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                      min_c = 0.5, max_c = 2, axis_c_divisor = 3, min_q = -1.5, max_q = 3.5, axis_q_divisor = 5,
                      y_lab = expression(log[10]*"(specific conductance)"))
    # LMP
    plot_logC_vs_logQ(df = dat_all, site = "LMP", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                      min_c = 1.5, max_c = 2.5, axis_c_divisor = 4, min_q = 1.5, max_q = 5, axis_q_divisor = 5, 
                      y_lab = expression(log[10]*"(specific conductance)"))
    # WHB
    plot_logC_vs_logQ(df = dat_all, site = "WHB", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                      min_c = 1.5, max_c = 3, axis_c_divisor = 3, min_q = 0, max_q = 3.5, axis_q_divisor = 5, 
                      y_lab = expression(log[10]*"(specific conductance)"))     
  