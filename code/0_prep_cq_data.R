# Prepare New Hampshire watershed data for C-Q analysis

  # Initial units: Q = L/s, NO3_corrected = mg N/L; DOC = mg C/L
  # This script prepares NO3 and DOC data

# Load libraries ----
  library("tidyverse")
  library("data.table")   # to read in certain columns from CSV


# Dates to filter on
  start_date = "2012-09-07"
  end_date = "2019-01-02"

# Functions ----
  # Function for hack to align facets with different year data by month-day
  same_year <- function(x) {
    year(x) <- 2000
    x
  }

# Read in 15-min and grab data from NH ----
  ## Lamprey main stem ----
  dat_lmp <-
    data.table::fread("data/LMP_WQual_Level4.csv",
                      select = c("DATETIME", "Site", "Q", "SpConductivity", "NO3_corrected", "FDOM_corrected_QSU", "Spec_Cond_grab", "Nitrate_grab", "FDOM_lab_grab", "NPOC_grab")) %>% 
    mutate(DATETIME = lubridate::mdy_hm(DATETIME)) %>%
    mutate(date = lubridate::date(DATETIME)) %>% 
    rename(datetime = DATETIME,
           site = Site,
           q = Q,
           spc = SpConductivity,
           no3 = NO3_corrected,
           fdom = FDOM_corrected_QSU,
           spc_grab = Spec_Cond_grab,
           no3_grab = Nitrate_grab,
           fdom_grab = FDOM_lab_grab,
           doc_grab = NPOC_grab) %>% 
    select(datetime, date, everything()) %>% 
    # There are two LMP sites (LMP and LMP72), but Adam Wymore says these are the same site
    mutate(site = ifelse(site == "LMP72", "LMP", site)) %>% 
    # Filter out dam release periods (same as we did for Wymore et al. 2023)
    filter(!(date >= ymd("2012-10-15") & date <= ymd("2012-11-15"))) %>%
    filter(!(date >= ymd("2013-10-15") & date <= ymd("2013-11-15"))) %>%
    filter(!(date >= ymd("2014-10-15") & date <= ymd("2014-11-15"))) %>%
    filter(!(date >= ymd("2015-10-15") & date <= ymd("2015-11-15"))) %>%
    filter(!(date >= ymd("2015-9-10") & date <= ymd("2015-10-15")))
  
    # unique(dat_lmp$site)

  ## DCF ----
  dat_dcf <-
    data.table::fread("data/DCF_WQual_Level4.csv",
                      select = c("DATETIME", "Site", "Q", "SpConductivity", "NO3_corrected", "FDOM_corrected_QSU", "Spec_Cond_grab", "Nitrate_grab", "FDOM_lab_grab", "NPOC_grab")) %>% 
    # mutate(DATETIME = lubridate::mdy_hm(DATETIME)) %>% 
    mutate(date = lubridate::date(DATETIME)) %>% 
    rename(datetime = DATETIME,
           site = Site,
           q = Q,
           spc = SpConductivity,
           no3 = NO3_corrected,
           fdom = FDOM_corrected_QSU,
           spc_grab = Spec_Cond_grab,
           no3_grab = Nitrate_grab,
           fdom_grab = FDOM_lab_grab,
           doc_grab = NPOC_grab) %>% 
    select(datetime, date, everything())
  
  ## HBF ----
  dat_hbf <-
    data.table::fread("data/HBF_WQual_Level4.csv",
                      select = c("DATETIME", "Site", "Q", "SpConductivity", "NO3_corrected", "FDOM_corrected_QSU", "Spec_Cond_grab", "Nitrate_grab", "FDOM_lab_grab", "NPOC_grab")) %>% 
    # mutate(DATETIME = lubridate::mdy_hm(DATETIME)) %>% 
    mutate(date = lubridate::date(DATETIME)) %>% 
    rename(datetime = DATETIME,
           site = Site,
           q = Q,
           spc = SpConductivity,
           no3 = NO3_corrected,
           fdom = FDOM_corrected_QSU,
           spc_grab = Spec_Cond_grab,
           no3_grab = Nitrate_grab,
           fdom_grab = FDOM_lab_grab,
           doc_grab = NPOC_grab) %>% 
    select(datetime, date, everything())
  
  ## BEF ----
    # Level 4 continuous data
    dat_bef_level4 <-
      data.table::fread("data/BEF_WQual_Level4.csv",
                        select = c("DATETIME", "Site", "Q", "SpConductivity", "NO3_corrected", "FDOM_corrected_QSU")) %>% 
      mutate(date = lubridate::date(DATETIME),
             DATETIME = round_date(DATETIME, unit = "15 mins")) %>% 
      rename(datetime = DATETIME,
             site = Site,
             q = Q,
             spc = SpConductivity,
             no3 = NO3_corrected,
             fdom = FDOM_corrected_QSU) %>% 
      select(datetime, date, everything())                          
  
    # Grab samples
    dat_bef_grab <-
      data.table::fread("data/CR1000_BEF_WQual_Level3_Grab.csv",
                        select = c("DATETIME", "Site", "Spec_Cond_grab", "Nitrate_grab", "FDOM_lab_grab", "NPOC_grab")) %>% 
      # mutate(DATETIME = lubridate::mdy_hm(DATETIME)) %>% 
      mutate(date = lubridate::date(DATETIME)) %>% 
      rename(datetime = DATETIME,
             site = Site,
             spc_grab = Spec_Cond_grab,
             no3_grab = Nitrate_grab,
             fdom_grab = FDOM_lab_grab,
             doc_grab = NPOC_grab) %>% 
      select(datetime, date, everything()) %>% 
      filter_at(vars(spc_grab, no3_grab, fdom_grab, doc_grab), any_vars(!is.na(.)))
    
    # Combine
    dat_bef <-
      dat_bef_level4 %>% 
      left_join(dat_bef_grab) %>% 
      # Remove two extreme precipitation events
      filter(!(datetime >= ymd_hms("2014-04-15 11:00:00") & date <= ymd_hms("2014-04-16 12:00:00"))) %>%
      filter(!(datetime >= ymd_hms("2014-06-25 19:00:00") & date <= ymd_hms("2014-06-27 23:00:00"))) %>%
      filter(!(datetime >= ymd_hms("2017-10-29 19:00:00") & date <= ymd_hms("2017-10-31 12:00:00"))) %>% 
      filter(!(datetime >= ymd_hms("2018-01-12 13:00:00") & date <= ymd_hms("2018-01-14 23:00:00")))
    
    # Check for duplicated timestamps
    nrow(any_dupes <-
      dat_bef %>%        
      group_by(datetime) %>% 
      filter(row_number() > 1))
    
    rm(dat_bef_level4, dat_bef_grab)
  
  ## BDC ----
    # Level 4 continuous data
    dat_bdc_level4 <-
      data.table::fread("data/BDC_WQual_Level4.csv",
                        select = c("DATETIME", "Site", "Q", "SpConductivity", "NO3_corrected", "FDOM_corrected_QSU")) %>% 
      mutate(date = lubridate::date(DATETIME),
             DATETIME = round_date(DATETIME, unit = "15 mins")) %>% 
      rename(datetime = DATETIME,
             site = Site,
             q = Q,
             spc = SpConductivity,
             no3 = NO3_corrected,
             fdom = FDOM_corrected_QSU) %>% 
      select(datetime, date, everything())        
    
    # Grab samples
    dat_bdc_grab <-
      data.table::fread("data/CR1000_BDC_WQual_Level3_Grab.csv",
                        select = c("DATETIME", "Site", "Spec_Cond_grab", "Nitrate_grab", "FDOM_lab_grab", "NPOC_grab")) %>% 
      mutate(date = lubridate::date(DATETIME)) %>% 
      rename(datetime = DATETIME,
             site = Site,
             spc_grab = Spec_Cond_grab,
             no3_grab = Nitrate_grab,
             fdom_grab = FDOM_lab_grab,
             doc_grab = NPOC_grab) %>% 
      select(datetime, date, everything()) %>% 
      filter_at(vars(spc_grab, no3_grab, fdom_grab, doc_grab), any_vars(!is.na(.)))
    
    # Combine
    dat_bdc <-
      dat_bdc_level4 %>% 
      left_join(dat_bdc_grab)
    
    # Check for duplicated timestamps
    nrow(any_dupes <-
           dat_bdc %>%        
           group_by(datetime) %>% 
           filter(row_number() > 1))
    
    rm(dat_bdc_level4, dat_bdc_grab)    
  
  ## WHB ----
  dat_whb <-
    data.table::fread("data/WHB_WQual_Level4.csv",
                      select = c("DATETIME", "Site", "Q", "SpConductivity", "NO3_corrected", "FDOM_corrected_QSU", "Spec_Cond_grab", "Nitrate_grab", "FDOM_lab_grab", "NPOC_grab")) %>% 
    # mutate(DATETIME = lubridate::ymd_hms(DATETIME)) %>% 
    mutate(date = lubridate::date(DATETIME)) %>% 
    rename(datetime = DATETIME,
           site = Site,
           q = Q,
           spc = SpConductivity,
           no3 = NO3_corrected,
           fdom = FDOM_corrected_QSU,
           spc_grab = Spec_Cond_grab,
           no3_grab = Nitrate_grab,
           fdom_grab = FDOM_lab_grab,
           doc_grab = NPOC_grab) %>% 
    select(datetime, date, everything())
  
# Read in derived (from fDOM) 15-min DOC time series ----
  # BDC
  doc_bdc <- 
    data.table::fread("data/BDC_DOC.csv") %>% 
    mutate(DATETIME = lubridate::mdy_hm(DATETIME),
           DATETIME = lubridate::round_date(DATETIME, unit = "15 mins"),
           date = lubridate::date(DATETIME)) %>% 
    select(datetime = DATETIME,
           site = Site,
           date,
           # fdom = FDOM_corrected_QSU,
           doc = DOC_mgC_L) %>% 
    filter(!is.na(doc))

  # BEF
  doc_bef <- 
    data.table::fread("data/BEF_DOC.csv") %>% 
    mutate(DATETIME = lubridate::mdy_hm(DATETIME),
           DATETIME = lubridate::round_date(DATETIME, unit = "15 mins"),
           date = lubridate::date(DATETIME)) %>% 
    select(datetime = DATETIME,
           site = Site,
           date,
           # fdom = FDOM_corrected_QSU,
           doc = DOC_mgC_L) %>% 
    filter(!is.na(doc))  
  
  # DCF
  doc_dcf <- 
    data.table::fread("data/DCF_DOC.csv") %>% 
    mutate(DATETIME = lubridate::mdy_hm(DATETIME),
           DATETIME = lubridate::round_date(DATETIME, unit = "15 mins"),
           date = lubridate::date(DATETIME)) %>% 
    select(datetime = DATETIME,
           site = Site,
           date,
           # fdom = FDOM_corrected_QSU,
           doc = DOC_mgC_L) %>% 
    filter(!is.na(doc))   
  
  # HBF
  doc_hbf <- 
    data.table::fread("data/HBF_DOC.csv") %>% 
    mutate(DATETIME = lubridate::mdy_hm(DATETIME),
           DATETIME = lubridate::round_date(DATETIME, unit = "15 mins"),
           date = lubridate::date(DATETIME)) %>% 
    select(datetime = DATETIME,
           site = Site,
           date,
           # fdom = FDOM_corrected_QSU,
           doc = DOC_mgC_L) %>% 
    filter(!is.na(doc)) 
  
  # LMP
  doc_lmp <- 
    data.table::fread("data/LMP_DOC.csv") %>% 
    mutate(DATETIME = lubridate::mdy_hm(DATETIME),
           DATETIME = lubridate::round_date(DATETIME, unit = "15 mins"),
           date = lubridate::date(DATETIME)) %>% 
    select(datetime = DATETIME,
           site = Site,
           date,
           # fdom = FDOM_corrected_QSU,
           doc = DOC_mgC_L) %>% 
    filter(!is.na(doc))
  
  # WHB
  doc_whb <- 
    data.table::fread("data/WHB_DOC.csv") %>% 
    mutate(DATETIME = lubridate::mdy_hm(DATETIME),
           DATETIME = lubridate::round_date(DATETIME, unit = "15 mins"),
           date = lubridate::date(DATETIME)) %>% 
    select(datetime = DATETIME,
           site = Site,
           date,
           # fdom = FDOM_corrected_QSU,
           doc = DOC_mgC_L) %>% 
    filter(!is.na(doc))
  
# Join DOC time series to other CQ data ----
  dat_bdc <-
    dat_bdc %>% 
    left_join(doc_bdc)

  dat_bef <-
    dat_bef %>% 
    left_join(doc_bef)
  
  dat_dcf <-
    dat_dcf %>% 
    left_join(doc_dcf)
  
  dat_hbf <-
    dat_hbf %>% 
    left_join(doc_hbf)
  
  dat_lmp <-
    dat_lmp %>% 
    left_join(doc_lmp)
  
  dat_whb <-
    dat_whb %>% 
    left_join(doc_whb)
  
  rm(doc_bdc, doc_bef, doc_dcf, doc_hbf, doc_lmp, doc_whb)
  
  
# Combine data frames into one & filter dates to desired range ----
  dat_all <-
    bind_rows(dat_bdc, dat_bef, dat_dcf, dat_hbf, dat_lmp, dat_whb) %>% 
    filter(date >= start_date & date <= end_date)
  
  rm(dat_bdc, dat_bef, dat_dcf, dat_hbf, dat_lmp, dat_whb)
  
  
# Examine data ----
  # Look for values < 0
  inspect <-
    dat_all %>% 
    pivot_longer(cols = c(q, spc, no3, fdom, doc, spc_grab, no3_grab, fdom_grab, doc_grab), names_to = "var", values_to = "value") %>%
    group_by(site, var) %>% 
    summarize(min = min(value, na.rm = TRUE),
              med = median(value, na.rm = TRUE),
              max = max(value, na.rm = TRUE))
  
  # BEF has INF and -INF for spc_grab
  bef <-
    dat_all %>% 
    filter(site == "BEF")
  # There are no grab sample data from SpC at this site
  
  # A couple of sites have values that look like 0 for min no3
  no3_min <-
    dat_all %>% 
    group_by(site) %>% 
    slice(which.min(no3))
  
  # View(no3_min)
  
  # How many values are 0?
  no3_zero <-
    dat_all %>% 
    filter(no3 == 0) %>% 
    group_by(site) %>% 
    tally()
  
  # View(no3_zero)
  
  # Let's set no3 0 values to 1/2 the minimum (non-zero) value
  no3_min_vals <-
    dat_all %>% 
    filter(no3 > 0,
           !is.na(no3)) %>% 
    group_by(site) %>% 
    summarize(no3_min = min(no3)) %>% 
    ungroup() %>% 
    mutate(no3_zero_correction = 0.5 * no3_min) %>% 
    select(-no3_min)
  
  # Replace no3 0 values with no3_zero_correction
  dat_all <- 
    dat_all %>% 
    left_join(no3_min_vals, by = join_by(site)) %>% 
    mutate(no3 = ifelse(no3 <= 0, no3_zero_correction, no3))
  
  # How many values are 0?
  no3_zero <-
    dat_all %>% 
    filter(no3 == 0) %>% 
    group_by(site) %>% 
    tally()
  
  # View(no3_zero)  
  

# Write dat_all to CSV
  dat_all %>% 
    write_csv(file = "data/nh_cq_data_filtered.csv")
  