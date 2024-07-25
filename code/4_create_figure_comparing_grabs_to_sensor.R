# Create figures comparing grab sample to sensor-derived data

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
  plot_grab_vs_sensor <- function(df, site, var, var_grab,  file_nm_var = "no3", max_c, axis_c_divisor, max_q, axis_q_divisor, y_lab) {
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
      # Get n
      n_grab <-
        grab %>% 
        tally() %>% 
        pull() %>% 
        as.character()

      # Plot
      p_cq_grab <-
        grab %>% 
        ggplot(aes(x = q, y = c)) +
        geom_point(shape = 1, alpha = 0.5) +
        scale_y_continuous(limits = c(0, {{ max_c }}), breaks = seq(0, {{ max_c }}, by = {{ max_c }} / axis_c_divisor)) +
        scale_x_continuous(limits = c(0, {{ max_q }}), breaks = seq(0, {{ max_q }}, by = {{ max_q }} / axis_q_divisor)) +
        theme_classic() +
        labs(subtitle = bquote(Grab~samples:~Weekly~(italic(n)==.(n_grab))),
             y = y_lab,
             x = expression(Discharge~(L~s^-1))) +
        theme(axis.title.x = element_blank())
  
    # 15-min sensor data
      # Get n
      n_15min <-
        sensor %>% 
        tally() %>% 
        pull() %>% 
        as.numeric() %>% 
        formatC(format="d", big.mark=",")
    
      # Plot
      p_cq_15min <-
        sensor %>% 
        ggplot(aes(x = q, y = c)) +
        geom_point(shape = 1, alpha = 0.2) +
        scale_y_continuous(limits = c(0, {{ max_c }}), breaks = seq(0, {{ max_c }}, by = {{ max_c }} / axis_c_divisor)) +
        scale_x_continuous(limits = c(0, {{ max_q }}), breaks = seq(0, {{ max_q }}, by = {{ max_q }} / axis_q_divisor)) +
        theme_classic() +
        labs(subtitle = bquote(Sensor:~"15-min"~(italic(n)==.(n_15min))),
             y = y_lab,
             x = expression(Discharge~(L~s^-1))) +
        theme(axis.title.y = element_blank())
      
    # Join plots together using patchwork
    ylab <- p_cq_grab$label$y
    p_cq_grab$labels$y <- p_cq_15min$labels$y <- " "

    # TIFF version
    tiff(paste0("plots/fig_CQ_grab_vs_sensor", "_", file_nm_var, "_", site, ".tiff"), height = 5, width = 3.5, units = "in", res = 600)
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
    summarize(min = min(val, na.rm = TRUE),
              max = max(val, na.rm = TRUE))
  
  dat_summ %>% 
    filter(var %in% c("q", "spc", "spc_grab"))  
  
  dat_summ %>% 
    filter(var %in% c("q", "doc", "doc_grab"))   
  

  dat_all %>% 
    filter(site == "WHB") %>% 
    filter(!is.na(spc_grab)) %>% 
    filter(!is.na(q)) %>% 
    rename(c = spc_grab) %>% 
    ggplot(aes(x = q, y = c)) +
    geom_point(shape = 1, alpha = 0.5) +
    theme_classic()
  
  # ggsave("plots/troubleshoot_bef_no3_grab.tiff", height = 3, width = 3, units = "in", dpi = 150)

  test <-
    dat_all %>% 
    filter(site == "BEF") %>% 
    filter(date > "2014-06-24" & date < "2014-06-28") %>% 
    ggplot(aes(x = datetime, y = q)) +
    geom_point(shape = 1) +
    theme_classic()
  
  dat_all %>% 
    filter(site == "BEF") %>% 
    filter(!is.na(no3)) %>% 
    filter(!is.na(q)) %>% 
    rename(c = no3) %>% 
    ggplot(aes(x = q, y = c)) +
    geom_point(shape = 1, alpha = 0.5) +
    theme_classic()
  
  # ggsave("plots/troubleshoot_bef_no3_sensor.tiff", height = 3, width = 3, units = "in", dpi = 150)
  
  # bef_test <-
  #   dat_all %>% 
  #   filter(site == "BEF") %>% 
  #   filter(q > 5000)   
  
  
# Create plots ----
  unique(dat_all$site)
  # "BDC" "BEF" "DCF" "HBF" "LMP" "WHB"
  
  ## NO3 ----
    # BDC
    plot_grab_vs_sensor(df = dat_all, site = "BDC", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                        max_c = 8, axis_c_divisor = 4, max_q = 160, axis_q_divisor = 5, 
                        y_lab = expression("NO"[3]~"-N (mg L"^-1~")"))
    # BEF
    plot_grab_vs_sensor(df = dat_all, site = "BEF", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                        max_c = 0.3, axis_c_divisor = 3, max_q = 30000, axis_q_divisor = 5, 
                        y_lab = expression("NO"[3]~"-N (mg L"^-1~")"))
    # DCF
    plot_grab_vs_sensor(df = dat_all, site = "DCF", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                        max_c = 0.8, axis_c_divisor = 4, max_q = 15000, axis_q_divisor = 5, 
                        y_lab = expression("NO"[3]~"-N (mg L"^-1~")"))
    # HBF
    plot_grab_vs_sensor(df = dat_all, site = "HBF", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                        max_c = 1, axis_c_divisor = 4, max_q = 1500, axis_q_divisor = 5, 
                        y_lab = expression("NO"[3]~"-N (mg L"^-1~")"))
    # LMP
    plot_grab_vs_sensor(df = dat_all, site = "LMP", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                        max_c = 1, axis_c_divisor = 4, max_q = 100000, axis_q_divisor = 5, 
                        y_lab = expression("NO"[3]~"-N (mg L"^-1~")"))
    # WHB
    plot_grab_vs_sensor(df = dat_all, site = "WHB", var = no3, var_grab = no3_grab, file_nm_var = "no3",
                        max_c = 2, axis_c_divisor = 4, max_q = 1600, axis_q_divisor = 4, 
                        y_lab = expression("NO"[3]~"-N (mg L"^-1~")"))
    
  ## DOC ----
    # BDC
    plot_grab_vs_sensor(df = dat_all, site = "BDC", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                        max_c = 35, axis_c_divisor = 5, max_q = 180, axis_q_divisor = 6, 
                        y_lab = expression("DOC (mg L"^-1~")"))
    # BEF
    plot_grab_vs_sensor(df = dat_all, site = "BEF", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                        max_c = 10, axis_c_divisor = 5, max_q = 30000, axis_q_divisor = 5, 
                        y_lab = expression("DOC (mg L"^-1~")"))
    # DCF
    plot_grab_vs_sensor(df = dat_all, site = "DCF", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                        max_c = 15, axis_c_divisor = 3, max_q = 15000, axis_q_divisor = 5, 
                        y_lab = expression("DOC (mg L"^-1~")"))
    # HBF
    plot_grab_vs_sensor(df = dat_all, site = "HBF", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                        max_c = 10, axis_c_divisor = 5, max_q = 1500, axis_q_divisor = 5, 
                        y_lab = expression("DOC (mg L"^-1~")"))
    # LMP
    plot_grab_vs_sensor(df = dat_all, site = "LMP", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                        max_c = 15, axis_c_divisor = 3, max_q = 100000, axis_q_divisor = 5, 
                        y_lab = expression("DOC (mg L"^-1~")"))
    # WHB
    plot_grab_vs_sensor(df = dat_all, site = "WHB", var = doc, var_grab = doc_grab, file_nm_var = "doc",
                        max_c = 18, axis_c_divisor = 6, max_q = 1600, axis_q_divisor = 4, 
                        y_lab = expression("DOC (mg L"^-1~")"))
    
  ## DOC ----
    # BDC
    plot_grab_vs_sensor(df = dat_all, site = "BDC", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                        max_c = 600, axis_c_divisor = 4, max_q = 180, axis_q_divisor = 6, 
                        y_lab = expression(Specific~conductance~(mu*S~cm^-1)))
    # BEF
    plot_grab_vs_sensor(df = dat_all, site = "BEF", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                        max_c = 40, axis_c_divisor = 4, max_q = 30000, axis_q_divisor = 5, 
                        y_lab = expression(Specific~conductance~(mu*S~cm^-1)))
    # DCF
    plot_grab_vs_sensor(df = dat_all, site = "DCF", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                        max_c = 700, axis_c_divisor = 4, max_q = 15000, axis_q_divisor = 5, 
                        y_lab = expression(Specific~conductance~(mu*S~cm^-1)))
    # HBF
    plot_grab_vs_sensor(df = dat_all, site = "HBF", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                        max_c = 60, axis_c_divisor = 4, max_q = 1500, axis_q_divisor = 5, 
                        y_lab = expression(Specific~conductance~(mu*S~cm^-1)))
    # LMP
    plot_grab_vs_sensor(df = dat_all, site = "LMP", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                        max_c = 280, axis_c_divisor = 4, max_q = 100000, axis_q_divisor = 5, 
                        y_lab = expression(Specific~conductance~(mu*S~cm^-1)))
    # WHB
    plot_grab_vs_sensor(df = dat_all, site = "WHB", var = spc, var_grab = spc_grab, file_nm_var = "spc",
                        max_c = 1600, axis_c_divisor = 4, max_q = 1600, axis_q_divisor = 4, 
                        y_lab = expression(Specific~conductance~(mu*S~cm^-1)))      
  