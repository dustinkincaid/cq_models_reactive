# Summarize (find the mean of wet deposition concentrations)
  
  # Notes:
    # The following sites should use data from LRHO (Lamprey watershed sites):
      # BDC, DCF, LMP
    # The following sites should use data from HBEF
      # BEF, HBF

# Load libraries ----
  library("tidyverse")
  library("data.table")

# Which dates will we filter on?
  start_date = "2012-09-07"
  end_date = "2019-01-02"

# Load data
  wetd_lrho <-
    data.table::fread("data/Wet Deposition LRHO TF NO3 DOC.csv") %>% 
    mutate(datetime = lubridate::mdy_hm(DateTime),
           date = lubridate::as_date(datetime))
  
  wetd_hbef <-
    data.table::fread("data/HBEF_RG22 NO3 DOC wet dep 2009 2020.csv") %>% 
    mutate(date = lubridate::mdy(date)) %>% 
    rename(no3_mgL = NO3,
           doc_mgCL = DOC) %>% 
    # Convert NO3 from mg/L to mg N/L
    mutate(no3_mgNL = no3_mgL * (14.007/62.005))
 
# Get the means after filtering the correct period
  mean_lrho <-
    wetd_lrho %>% 
    filter(date >= start_date & date <= end_date) %>% 
    summarize(no3_mgNL = mean(NO3mgNL, na.rm = TRUE),
              doc_mgCL = mean(DOCmgCL, na.rm = TRUE)) %>% 
    mutate(source = "lrho")
  
  mean_hbef <-
    wetd_hbef %>% 
    filter(date >= start_date & date <= end_date) %>% 
    summarize(no3_mgNL = mean(no3_mgNL, na.rm = TRUE),
              doc_mgCL = mean(doc_mgCL, na.rm = TRUE)) %>% 
    mutate(source = "hbef")  

# Combine and write to CSV
  mean_lrho %>% 
    bind_rows(mean_hbef) %>% 
    select(source, everything()) %>% 
    write_csv(file = "data/wet_dep_mean_concs.csv")
  