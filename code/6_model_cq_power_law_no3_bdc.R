# Fit the CQ data to Hall's models

# Load libraries ----
  library("tidyverse")    # general workhorse
  library("here")         # takes the guesswork out of dealing with file paths
  library("data.table")   # read in large datasets faster
  library("reshape2")     # manipulate dataframes
  library("patchwork")    # makes laying out plots much easier
  library("grid")         # more plotting assistance
  library("rsample")      # for model cross-validation
  library("doParallel")   # for faster computing
  library("scico")        # color palette for ggplot2
  library("minpack.lm")   # alternative to 'nls'; much more flexible!
  # library("scales")       # used for plotting on log scales

# Which site and solute would you like to model?
  site_choice <- "BDC"
  sol_choice <- "no3"
  grab_choice <- "no3_grab"

# Options ----
  options(scipen = 999) #disable printing of numbers in scientific notation
  
# Functions for model fitting and plotting ----
  RMSE <- function(error) { sqrt(mean(error^2)) }
  
  KirchnerBinning = function(df, min_per_bin = 20){
    
    logQ = log(df$q)
    
    logRange = max(logQ) - min(logQ)
    minBinSize = logRange*.01
    
    binBoundaries = c(1)
    for (i in 2:dim(df)[1]){ 
      if (abs(logQ[i] - logQ[tail(binBoundaries,n=1)]) < minBinSize){
        next
      }
      
      if (abs(i-tail(binBoundaries,n=1)) < min_per_bin){
        next
      }
      
      curr = as.numeric(unlist(df[tail(binBoundaries,n=1):i,'c']))
      
      if (sd(curr,na.rm=TRUE)/sqrt(abs(i-tail(binBoundaries,n=1))) > mean(curr)/2){
        next
      }
      
      binBoundaries=c(binBoundaries,i)
    }
    return(binBoundaries)
  }  
  
# Read in data ----
  ## CQ data from 1_prep_cq_data.R ----
  dat_all <- 
    data.table::fread(here("data", "nh_cq_data_filtered.csv"))
  
  ## C0 choices
    # Baseflow
    C0_baseflow <-
      read_csv(here("data", "c0_baseflow_estimates.csv"))
    
    # Wet deposition
    C0_wetdep <-
      read_csv(here("data", "wet_dep_mean_concs.csv")) %>% 
      full_join(tibble(site = c("BDC", "DCF", "LMP", "WHB", "BEF", "HBF"),
                       source = c("lrho", "lrho", "lrho", "lrho", "hbef", "hbef")))
    
  ## Recession 'n' values
    rec_n_all <-
      read_csv(here("data", "recessionA.csv")) %>% 
      rename(site = Site)
    
# Pull C0 and n values ----
  c0_bf <- 
    C0_baseflow %>% 
    filter(site == site_choice,
           var == sol_choice) %>% 
    pull()
  
  c0_wd <-
    C0_wetdep %>% 
    rename(no3 = no3_mgNL,
           doc = doc_mgCL) %>% 
    filter(site == site_choice) %>% 
    pull(sol_choice)
  
  rec_n <-
    rec_n_all %>% 
    filter(site == site_choice) %>% 
    pull(n)

# Create grab sample and sensor datasets ----
  # Subset for the site  
  df_sub <-
    dat_all %>% 
    filter(site == site_choice)    
    
  # Filter for 15-min sensor data
  df_15min <-
    df_sub %>% 
    filter(!is.na(!!as.symbol(sol_choice))) %>% 
    filter(!is.na(q)) %>% 
    rename(c = !!as.symbol(sol_choice))
  
  # Filter for long-term grab sample data
  df_grab <-
    df_sub %>% 
    filter(!is.na(!!as.symbol(grab_choice))) %>% 
    filter(!is.na(q)) %>% 
    rename(c = !!as.symbol(grab_choice)) 
  
# Reduce sensor data to the various sampling times (ie, daily, weekly, monthly) ----
  # Here we create 10 versions of each sampling frequency (e.g., df_15min_daily_1 through df_15min_daily_10)  
  
  # Daily (sample always occurs sometime between 8am and 5pm)
  df_15min_daily_list <- list()
  for(i in 1:10){
    df_15min_daily_list[[i]] <-
      df_15min %>% 
      mutate(period = ifelse(lubridate::hour(datetime) %in% c(8:16), "sampling", "not_sampling")) %>% 
      filter(period == "sampling") %>% 
      group_by(date) %>% 
      # This was randomly sample a row from each group (here, date)
      slice_sample(n = 1) %>% 
      mutate(df = "daily") %>% 
      ungroup()
  }

  # Weekly (sample always occurs M-F sometime between 8am and 5pm)
  df_15min_weekly_list <- list()
  for(i in 1:10){
    df_15min_weekly_list[[i]] <-
      df_15min %>% 
      mutate(period = ifelse(lubridate::hour(datetime) %in% c(8:16), "sampling", "not_sampling")) %>% 
      filter(period == "sampling") %>% 
      mutate(wday = lubridate::wday(date, label = TRUE)) %>% 
      filter(!wday %in% c("Sat", "Sun")) %>% 
      group_by(lubridate::year(date), lubridate::week(date)) %>% 
      # This was randomly sample a row from each group (here, date)
      slice_sample(n = 1) %>% 
      mutate(df = "weekly") %>% 
      ungroup() %>% 
      select(-c(`lubridate::year(date)`, `lubridate::week(date)`))
  }    

  # Monthly (sample always occurs M-F sometime between 8am and 5pm)
  df_15min_monthly_list <- list()
  for(i in 1:10){
    df_15min_monthly_list[[i]] <-
      df_15min %>% 
      mutate(period = ifelse(lubridate::hour(datetime) %in% c(8:16), "sampling", "not_sampling")) %>% 
      filter(period == "sampling") %>% 
      mutate(wday = lubridate::wday(date, label = TRUE)) %>% 
      filter(!wday %in% c("Sat", "Sun")) %>% 
      group_by(lubridate::year(date), lubridate::month(date)) %>% 
      # This was randomly sample a row from each group (here, date)
      slice_sample(n = 1) %>% 
      mutate(df = "monthly") %>% 
      ungroup() %>% 
      select(-c(`lubridate::year(date)`, `lubridate::month(date)`))
  }    
  

# Fit the models with each dataset  ----
  ## Create list of datasets ----
  df_list <- list(df_grab %>% mutate(df = "lt_grab"),
                  df_15min %>% mutate(df = "15min"),
                  df_15min_daily_list[[1]],
                  df_15min_weekly_list[[1]],
                  df_15min_monthly_list[[1]])
  
  ## Fit models, save residuals, and create plot df ----
    ### c0 = c0_bf ----
    c0_set <- c0_bf
  
    # i = 1
    resids_final <- tibble()
    bdf_all <- tibble()
    plot_df_all <- tibble()
    
    for(i in 1:length(df_list)) {
        
        print(as.character(df_list[[i]][1, "df"]))
      
        # Equation 0: Regress log10(C) on log(Q)
        eq0_fit <- lm(log10(c) ~ log10(q), data=df_list[[i]])
        eq0_resids <- 
          tibble(resids = as.numeric(residuals(eq0_fit)),
                 q = df_list[[i]]$q,
                 model = "eq0",
                 df = df_list[[i]]$df)
        print(sprintf("RMSE Eq 0: %s", RMSE(summary(eq0_fit)$residuals)))    
        
        # Equation 1/Model 1: Single mixing volume, constant load, inflow has zero concentration
          # Force n
          eq1 <- function(q, a, n) (a*q**(-1/n))
          eq1_fit <- nlsLM(c ~ eq1(q, a, n=rec_n), data=df_list[[i]], start=list(a=380))
          eq1_resids <- 
            tibble(resids = as.numeric(residuals(eq1_fit)),
                   q = df_list[[i]]$q,
                   model = "eq1",
                   param_a = eq1_fit$m$getPars()[[1]],
                   n_source = "set",
                   param_n = rec_n,
                   df = df_list[[i]]$df)
          print(sprintf("RMSE Eq 1-1: %s", RMSE(summary(eq1_fit)$residuals)))
          
          # Don't force n
          eq1_nFit_fit <- nlsLM(c ~ eq1(q, a, n), data=df_list[[i]], start=list(a=380, n=0.1))
          eq1_nFit_resids <- 
            tibble(resids = as.numeric(residuals(eq1_nFit_fit)),
                   q = df_list[[i]]$q,
                   model = "eq1_nFit",
                   param_a = eq1_nFit_fit$m$getPars()[[1]],
                   n_source = "fit",
                   param_n = eq1_nFit_fit$m$getPars()[[2]],
                   df = df_list[[i]]$df)
          print(sprintf("RMSE Eq 1-2: %s", RMSE(summary(eq1_nFit_fit)$residuals)))          
        
        # Equation 2/Model 2: Single mixing volume, constant load, inflow has constant concentration
          # Force n
          eq2 <- function(q, a, n, c0) (a*q**(-1/n)+c0)
          eq2_fit <- nlsLM(c ~ eq2(q, a, n=rec_n, c0=c0_set), data=df_list[[i]], start=list(a=380))
          eq2_resids <- 
            tibble(resids = as.numeric(residuals(eq2_fit)),
                   q = df_list[[i]]$q,
                   model = "eq2",
                   param_a = eq2_fit$m$getPars()[[1]],
                   n_source = "set",
                   param_n = rec_n,
                   df = df_list[[i]]$df)  
          print(sprintf("RMSE Eq 2-1: %s", RMSE(summary(eq2_fit)$residuals)))
          
          # Don't force n
          eq2_nFit_fit <- nlsLM(c ~ eq2(q, a, n, c0=c0_set), data=df_list[[i]], start=list(a=380, n=1))
          eq2_nFit_resids <- 
            tibble(resids = as.numeric(residuals(eq2_nFit_fit)),
                   q = df_list[[i]]$q,
                   model = "eq2_nFit",
                   param_a = eq2_nFit_fit$m$getPars()[[1]],
                   n_source = "fit",
                   param_n = eq2_nFit_fit$m$getPars()[[2]],
                   df = df_list[[i]]$df)  
          print(sprintf("RMSE Eq 2-1: %s", RMSE(summary(eq2_nFit_fit)$residuals)))          
        
        # Equation 3/Model 3a: Same as models 1 or 2 except one assumes concentration changes little as compared to volume
        eq3 <- function(q, f,e) (f-e*log(q))
        eq3_fit <- nlsLM(c ~ eq3(q, f, e), data=df_list[[i]], start=list(f=100, e=100))
        eq3_resids <- 
          tibble(resids = as.numeric(residuals(eq3_fit)),
                 q = df_list[[i]]$q,
                 model = "eq3",
                 param_f = eq3_fit$m$getPars()[[1]],
                 param_e = eq3_fit$m$getPars()[[2]],
                 df = df_list[[i]]$df)   
        print(sprintf("RMSE Eq 3: %s", RMSE(summary(eq3_fit)$residuals)))
        
        # Equations 4,4a/Model 3b: Same as model 1 except one assumes volume changes little as compared to concentration
          # Force n
          eq4 <- function(q, h, g, n) (h*exp(-g*q**(1/n)))
          eq4_fit <- nlsLM(c ~ eq4(q, h, g, n=rec_n), data=df_list[[i]], start=list(h=175, g=0.00001))
          eq4_resids <- 
            tibble(resids = as.numeric(residuals(eq4_fit)),
                   q = df_list[[i]]$q,
                   model = "eq4",
                   param_h = eq4_fit$m$getPars()[[1]],
                   param_g = eq4_fit$m$getPars()[[2]],
                   n_source = "set",
                   param_n = rec_n,
                   df = df_list[[i]]$df)   
          print(sprintf("RMSE Eq 4-1: %s", RMSE(summary(eq4_fit)$residuals)))
          
          # Don't force n
          eq4_nFit_fit <- nlsLM(c ~ eq4(q, h, g, n), data=df_list[[i]], start=list(h=175, g=0.00001, n=0.1))
          eq4_nFit_resids <- 
            tibble(resids = as.numeric(residuals(eq4_nFit_fit)),
                   q = df_list[[i]]$q,
                   model = "eq4_nFit",
                   param_h = eq4_nFit_fit$m$getPars()[[1]],
                   param_g = eq4_nFit_fit$m$getPars()[[2]],
                   n_source = "fit",
                   param_n = eq4_nFit_fit$m$getPars()[[3]],
                   df = df_list[[i]]$df)   
          print(sprintf("RMSE Eq 4-2: %s", RMSE(summary(eq4_nFit_fit)$residuals)))          
        
        # Equations 5, 5a/Model 3c: Same as model 2 except one assumes volume changes little as compared to concentration
          # Force n
          eq5 <- function(q, h, g, n, c0) (h*exp(-g*q**(1/n))+c0)
          eq5_fit <- nlsLM(c ~ eq5(q, h, g, n=rec_n, c0=c0_set), data=df_list[[i]], start=list(h=160, g=0.00001))
          eq5_resids <- 
            tibble(resids = as.numeric(residuals(eq5_fit)),
                   q = df_list[[i]]$q,
                   model = "eq5",
                   param_h = eq5_fit$m$getPars()[[1]],
                   param_g = eq5_fit$m$getPars()[[2]],
                   n_source = "set",
                   param_n = rec_n,
                   df = df_list[[i]]$df)   
          print(sprintf("RMSE Eq 5-1: %s", RMSE(summary(eq5_fit)$residuals)))
          
          # Don't force n
          eq5_nFit_fit <- nlsLM(c ~ eq5(q, h, g, n, c0=c0_set), data=df_list[[i]], start=list(h=160, g=0.00001, n=1))
          eq5_nFit_resids <- 
            tibble(resids = as.numeric(residuals(eq5_nFit_fit)),
                   q = df_list[[i]]$q,
                   model = "eq5_nFit",
                   param_h = eq5_nFit_fit$m$getPars()[[1]],
                   param_g = eq5_nFit_fit$m$getPars()[[2]],
                   n_source = "fit",
                   param_n = eq5_nFit_fit$m$getPars()[[3]],
                   df = df_list[[i]]$df)   
          print(sprintf("RMSE Eq 5-2: %s", RMSE(summary(eq5_nFit_fit)$residuals)))          
        
        # Equation 6, model 3d: same as 1 and 2, neither volume nor concentration changes much
          # Force n
          eq6 <- function(q, m, j, n) (m-j*q**(1/n))
          eq6_fit <- nlsLM(c ~ eq6(q, m, j, n=rec_n), data=df_list[[i]], start=list(m=100, j=100))
          eq6_resids <- 
            tibble(resids = as.numeric(residuals(eq6_fit)),
                   q = df_list[[i]]$q,
                   model = "eq6",
                   param_m = eq6_fit$m$getPars()[[1]],
                   param_j = eq6_fit$m$getPars()[[2]],
                   n_source = "set",
                   param_n = rec_n,
                   df = df_list[[i]]$df) 
          print(sprintf("RMSE Eq 6-1: %s", RMSE(summary(eq6_fit)$residuals)))
          
          # Don't force n
          eq6_nFit_fit <- nlsLM(c ~ eq6(q, m, j, n), data=df_list[[i]], start=list(m=100, j=100, n = 0.1))
          eq6_nFit_resids <- 
            tibble(resids = as.numeric(residuals(eq6_nFit_fit)),
                   q = df_list[[i]]$q,
                   model = "eq6_nFit",
                   param_m = eq6_nFit_fit$m$getPars()[[1]],
                   param_j = eq6_nFit_fit$m$getPars()[[2]],
                   n_source = "fit",
                   param_n = eq6_nFit_fit$m$getPars()[[3]],
                   df = df_list[[i]]$df) 
          print(sprintf("RMSE Eq 6-2: %s", RMSE(summary(eq6_nFit_fit)$residuals)))          
        
        # Equation 7/Model 4: Single mixing volume; constant load; all outflow from inflow, and inflow has zero concentration
          # Don't force n
          eq7 <- function(q, s, b, n) (s/(1+(b*q**(1/n))))
          eq7_fit <- nlsLM(c ~ eq7(q, s, b, n=rec_n), data=df_list[[i]],start=list(s=175, b=0.00001))
          eq7_resids <- 
            tibble(resids = as.numeric(residuals(eq7_fit)),
                   q = df_list[[i]]$q,
                   model = "eq7",
                   param_s = eq7_fit$m$getPars()[[1]],
                   param_b = eq7_fit$m$getPars()[[2]],
                   n_source = "set",
                   param_n = rec_n,
                   df = df_list[[i]]$df) 
          print(sprintf("RMSE Eq 7-1: %s", RMSE(summary(eq7_fit)$residuals)))
          
          # Don't force n
          eq7_nFit_fit <- nlsLM(c ~ eq7(q, s, b, n), data=df_list[[i]],start=list(s=175, b=0.00001, n = 0.1))
          eq7_nFit_resids <- 
            tibble(resids = as.numeric(residuals(eq7_nFit_fit)),
                   q = df_list[[i]]$q,
                   model = "eq7_nFit",
                   param_s = eq7_nFit_fit$m$getPars()[[1]],
                   param_b = eq7_nFit_fit$m$getPars()[[2]],
                   n_source = "fit",
                   param_n = eq7_nFit_fit$m$getPars()[[3]],
                   df = df_list[[i]]$df) 
          print(sprintf("RMSE Eq 7-2: %s", RMSE(summary(eq7_nFit_fit)$residuals)))          
        
        # Equation 8/Model 5: single mixing volume, inflow = outflow, some part of mixing volume remains constant (V0). Inflow has constant concentration (C0)
          # Force n
          eq8 <- function(q, s, b, n, c0) ((s-c0)/(1+(b*q**(1/n)))+c0)
          eq8_fit <- nlsLM(c ~ eq8(q, s, b, n=rec_n, c0=c0_set), data=df_list[[i]],start=list(s=175, b=0.00001))
          eq8_resids <- 
            tibble(resids = as.numeric(residuals(eq8_fit)),
                   q = df_list[[i]]$q,
                   model = "eq8",
                   param_s = eq8_fit$m$getPars()[[1]],
                   param_b = eq8_fit$m$getPars()[[2]],
                   n_source = "set",
                   param_n = rec_n,
                   df = df_list[[i]]$df)   
          print(sprintf("RMSE Eq 8-1: %s", RMSE(summary(eq8_fit)$residuals)))
          
          # Don't force n
          eq8_nFit_fit <- nlsLM(c ~ eq8(q, s, b, n, c0=c0_set), data=df_list[[i]],start=list(s=175, b=0.00001, n=2))
          eq8_nFit_resids <- 
            tibble(resids = as.numeric(residuals(eq8_nFit_fit)),
                   q = df_list[[i]]$q,
                   model = "eq8_nFit",
                   param_s = eq8_nFit_fit$m$getPars()[[1]],
                   param_b = eq8_nFit_fit$m$getPars()[[2]],
                   n_source = "fit",
                   param_n = eq8_nFit_fit$m$getPars()[[3]],
                   df = df_list[[i]]$df)   
          print(sprintf("RMSE Eq 8-2: %s", RMSE(summary(eq8_nFit_fit)$residuals)))          
        
          # Combine all residual dfs
          resids_final <- 
            bind_rows(resids_final,
                      eq0_resids,
                      eq1_resids,
                      eq1_nFit_resids,
                      eq2_resids,
                      eq2_nFit_resids,
                      eq3_resids, 
                      eq4_resids,
                      eq4_nFit_resids,
                      eq5_resids,
                      eq5_nFit_resids,
                      eq6_resids,
                      eq6_nFit_resids,
                      eq7_resids,
                      eq7_nFit_resids,
                      eq8_resids,
                      eq8_nFit_resids) %>% 
            mutate(c0_source = "baseflow",
                   c0_value = c0_bf) %>% 
            select(df, model, c0_source, c0_value, q, resids, n_source, param_n, param_a, 
                   param_b, param_e, param_f, param_g, param_h, param_j, param_m, param_s)
        
          # Plot model fits
            df2 = df_list[[i]][order(df_list[[i]]$q,decreasing=TRUE),]
            
            binBoundaries = KirchnerBinning(df2, min_per_bin=20)
            
            bqs = c()
            bcs = c()
            sigmas = c()
            for (j in 1:(length(binBoundaries)-1)){
              bqs = c(bqs,mean(df2$q[binBoundaries[j]:binBoundaries[j+1]],na.rm=TRUE))
              bcs = c(bcs,mean(df2$c[binBoundaries[j]:binBoundaries[j+1]],na.rm=TRUE))
              sigmas = c(sigmas,sd(df2$c[binBoundaries[j]:binBoundaries[j+1]],na.rm=TRUE)/sqrt(binBoundaries[j+1]-binBoundaries[j]))
            }  
            
            # get dfs ready for ggplot
            bdf <- data.frame(bqs,bcs,sigmas) %>% 
              mutate(df = as.character(df_list[[i]][1, "df"]))
            plotdf = data.frame(df_list[[i]]$q, predict(eq0_fit), predict(eq1_fit), predict(eq1_nFit_fit),
                                predict(eq2_fit), predict(eq2_nFit_fit), predict(eq3_fit), predict(eq4_fit), 
                                predict(eq4_nFit_fit), predict(eq5_fit), predict(eq5_nFit_fit), predict(eq6_fit),
                                predict(eq6_nFit_fit), predict(eq7_fit), predict(eq7_nFit_fit), predict(eq8_fit), predict(eq8_nFit_fit))
            colnames(plotdf) = c('q','eq0','eq1','eq1_nFit','eq2','eq2_nFit','eq3','eq4','eq4_nFit',
                                 'eq5','eq5_nFit','eq6','eq6_nFit','eq7','eq7_nFit','eq8','eq8_nFit')
            
            plotdf <- melt(plotdf,  id.vars = 'q', variable.name = 'model') %>% 
              mutate(df = as.character(df_list[[i]][1, "df"]))
            
            # Aggregate dfs
            bdf_all <- bind_rows(bdf_all, bdf) %>% 
              mutate(c0_source = "baseflow",
                     c0_value = c0_bf)
            plot_df_all <- bind_rows(plot_df_all, plotdf) %>% 
              mutate(c0_source = "baseflow",
                     c0_value = c0_bf)
        
    }
    
      # Remove objects
      rm(eq0_fit, eq1_fit, eq1_nFit_fit, eq2_fit, eq2_nFit_fit, eq3_fit, eq4_fit,
         eq4_nFit_fit, eq5_fit, eq5_nFit_fit, eq6_fit, eq6_nFit_fit, eq7_fit,
         eq7_nFit_fit, eq8_fit, eq8_nFit_fit,
         eq0_resids, eq1_resids, eq1_nFit_resids, eq2_resids, eq2_nFit_resids, eq3_resids, eq4_resids,
         eq4_nFit_resids, eq5_resids, eq5_nFit_resids, eq6_resids, eq6_nFit_resids, eq7_resids,
         eq7_nFit_resids, eq8_resids, eq8_nFit_resids)    
  
    
    ### REPEAT where c0 = c0_wd ----
    c0_set <- c0_wd
    
    # i = 1
    resids_final_2 <- tibble()
    bdf_all_2 <- tibble()
    plot_df_all_2 <- tibble()
    
    for(i in 1:length(df_list)) {
      
      print(as.character(df_list[[i]][1, "df"]))

      # Equation 2/Model 2: Single mixing volume, constant load, inflow has constant concentration
      # Force n
      eq2 <- function(q, a, n, c0) (a*q**(-1/n)+c0)
      eq2_fit <- nlsLM(c ~ eq2(q, a, n=rec_n, c0=c0_set), data=df_list[[i]], start=list(a=380))
      eq2_resids <- 
        tibble(resids = as.numeric(residuals(eq2_fit)),
               q = df_list[[i]]$q,
               model = "eq2",
               param_a = eq2_fit$m$getPars()[[1]],
               n_source = "set",
               param_n = rec_n,
               df = df_list[[i]]$df)  
      print(sprintf("RMSE Eq 2-1: %s", RMSE(summary(eq2_fit)$residuals)))
      
      # Don't force n
      eq2_nFit_fit <- nlsLM(c ~ eq2(q, a, n, c0=c0_set), data=df_list[[i]], start=list(a=380, n=1))
      eq2_nFit_resids <- 
        tibble(resids = as.numeric(residuals(eq2_nFit_fit)),
               q = df_list[[i]]$q,
               model = "eq2_nFit",
               param_a = eq2_nFit_fit$m$getPars()[[1]],
               n_source = "fit",
               param_n = eq2_nFit_fit$m$getPars()[[2]],
               df = df_list[[i]]$df)  
      print(sprintf("RMSE Eq 2-1: %s", RMSE(summary(eq2_nFit_fit)$residuals)))          
      
      # Equations 5, 5a/Model 3c: Same as model 2 except one assumes volume changes little as compared to concentration
      # Force n
      eq5 <- function(q, h, g, n, c0) (h*exp(-g*q**(1/n))+c0)
      eq5_fit <- nlsLM(c ~ eq5(q, h, g, n=rec_n, c0=c0_set), data=df_list[[i]], start=list(h=160, g=0.00001))
      eq5_resids <- 
        tibble(resids = as.numeric(residuals(eq5_fit)),
               q = df_list[[i]]$q,
               model = "eq5",
               param_h = eq5_fit$m$getPars()[[1]],
               param_g = eq5_fit$m$getPars()[[2]],
               n_source = "set",
               param_n = rec_n,
               df = df_list[[i]]$df)   
      print(sprintf("RMSE Eq 5-1: %s", RMSE(summary(eq5_fit)$residuals)))
      
      # Don't force n
      eq5_nFit_fit <- nlsLM(c ~ eq5(q, h, g, n, c0=c0_set), data=df_list[[i]], start=list(h=160, g=0.00001, n=1))
      eq5_nFit_resids <- 
        tibble(resids = as.numeric(residuals(eq5_nFit_fit)),
               q = df_list[[i]]$q,
               model = "eq5_nFit",
               param_h = eq5_nFit_fit$m$getPars()[[1]],
               param_g = eq5_nFit_fit$m$getPars()[[2]],
               n_source = "fit",
               param_n = eq5_nFit_fit$m$getPars()[[3]],
               df = df_list[[i]]$df)   
      print(sprintf("RMSE Eq 5-2: %s", RMSE(summary(eq5_nFit_fit)$residuals)))          
      
      # Equation 8/Model 5: single mixing volume, inflow = outflow, some part of mixing volume remains constant (V0). Inflow has constant concentration (C0)
      # Force n
      eq8 <- function(q, s, b, n, c0) ((s-c0)/(1+(b*q**(1/n)))+c0)
      eq8_fit <- nlsLM(c ~ eq8(q, s, b, n=rec_n, c0=c0_set), data=df_list[[i]],start=list(s=175, b=0.00001))
      eq8_resids <- 
        tibble(resids = as.numeric(residuals(eq8_fit)),
               q = df_list[[i]]$q,
               model = "eq8",
               param_s = eq8_fit$m$getPars()[[1]],
               param_b = eq8_fit$m$getPars()[[2]],
               n_source = "set",
               param_n = rec_n,
               df = df_list[[i]]$df)   
      print(sprintf("RMSE Eq 8-1: %s", RMSE(summary(eq8_fit)$residuals)))
      
      # Don't force n
      eq8_nFit_fit <- nlsLM(c ~ eq8(q, s, b, n, c0=c0_set), data=df_list[[i]],start=list(s=175, b=0.00001, n=2))
      eq8_nFit_resids <- 
        tibble(resids = as.numeric(residuals(eq8_nFit_fit)),
               q = df_list[[i]]$q,
               model = "eq8_nFit",
               param_s = eq8_nFit_fit$m$getPars()[[1]],
               param_b = eq8_nFit_fit$m$getPars()[[2]],
               n_source = "fit",
               param_n = eq8_nFit_fit$m$getPars()[[3]],
               df = df_list[[i]]$df)   
      print(sprintf("RMSE Eq 8-2: %s", RMSE(summary(eq8_nFit_fit)$residuals)))          
      
      # Combine all residual dfs
      resids_final_2 <- 
        bind_rows(resids_final_2,
                  eq2_resids,
                  eq2_nFit_resids,
                  eq5_resids,
                  eq5_nFit_resids,
                  eq8_resids,
                  eq8_nFit_resids) %>% 
        mutate(c0_source = "wet_deposition",
               c0_value = c0_wd) %>% 
        select(df, model, c0_source, c0_value, q, resids, n_source, param_n, param_a, param_b, param_g, param_h, param_s)
      
      # Plot model fits
      df2 = df_list[[i]][order(df_list[[i]]$q,decreasing=TRUE),]
      
      binBoundaries = KirchnerBinning(df2, min_per_bin=20)
      
      bqs = c()
      bcs = c()
      sigmas = c()
      for (j in 1:(length(binBoundaries)-1)){
        bqs = c(bqs,mean(df2$q[binBoundaries[j]:binBoundaries[j+1]],na.rm=TRUE))
        bcs = c(bcs,mean(df2$c[binBoundaries[j]:binBoundaries[j+1]],na.rm=TRUE))
        sigmas = c(sigmas,sd(df2$c[binBoundaries[j]:binBoundaries[j+1]],na.rm=TRUE)/sqrt(binBoundaries[j+1]-binBoundaries[j]))
      }  
      
      # get dfs ready for ggplot
      bdf <- data.frame(bqs,bcs,sigmas) %>% 
        mutate(df = as.character(df_list[[i]][1, "df"]))
      plotdf = data.frame(df_list[[i]]$q, predict(eq0_fit), predict(eq1_fit), predict(eq1_nFit_fit),
                          predict(eq2_fit), predict(eq2_nFit_fit), predict(eq3_fit), predict(eq4_fit), 
                          predict(eq4_nFit_fit), predict(eq5_fit), predict(eq5_nFit_fit), predict(eq6_fit),
                          predict(eq6_nFit_fit), predict(eq7_fit), predict(eq7_nFit_fit), predict(eq8_fit), predict(eq8_nFit_fit))
      colnames(plotdf) = c('q','eq0','eq1','eq1_nFit','eq2','eq2_nFit','eq3','eq4','eq4_nFit',
                           'eq5','eq5_nFit','eq6','eq6_nFit','eq7','eq7_nFit','eq8','eq8_nFit')
      
      plotdf <- melt(plotdf,  id.vars = 'q', variable.name = 'model') %>% 
        mutate(df = as.character(df_list[[i]][1, "df"]))
      
      # Aggregate dfs
      bdf_all_2 <- bind_rows(bdf_all_2, bdf) %>% 
        mutate(c0_source = "wet_deposition",
               c0_value = c0_wd)
      plot_df_all_2 <- bind_rows(plot_df_all_2, plotdf) %>% 
        mutate(c0_source = "wet_deposition",
               c0_value = c0_wd)      
    }        
  
  
  ## Evaluate model performance using repeated K/V-fold cross-validation ----    
    ### Create function to fit and assess models
      # All models
      cv_mods <- function(split, ...){
        ### fit models with the analysis partition
        eq1_fit <- nlsLM(c ~ eq1(q, a, n=rec_n), data=rsample::analysis(split), start=list(a=380))
        eq2_fit <- nlsLM(c ~ eq2(q, a, n=rec_n, c0=c0_set), data=rsample::analysis(split), start=list(a=380))
        eq3_fit <- nlsLM(c ~ eq3(q, f, e), data=rsample::analysis(split), start=list(f=100, e=100))
        eq4_fit <- nlsLM(c ~ eq4(q, h, g, n=rec_n), data=rsample::analysis(split), start=list(h=175, g=0.00001))
        eq5_fit <- nlsLM(c ~ eq5(q, h, g, n=rec_n, c0=c0_set), data=rsample::analysis(split), start=list(h=160, g=0.00001))
        eq6_fit <- nlsLM(c ~ eq6(q, m, j, n=rec_n), data=rsample::analysis(split), start=list(m=100, j=100))
        eq7_fit <- nlsLM(c ~ eq7(q, s, b, n=rec_n), data=rsample::analysis(split),start=list(s=175, b=0.00001))
        eq8_fit <- nlsLM(c ~ eq8(q, s, b, n=rec_n, c0=c0_set), data=rsample::analysis(split),start=list(s=175, b=0.00001))
        ### take the dependent variable from the assessment partition
        y <- rsample::assessment(split)$c
        ### create the residuals  getting model predictions from the assessment partition
        e1 <- (y - predict(eq1_fit, newdata=rsample::assessment(split)))
        e2 <- (y - predict(eq2_fit, newdata=rsample::assessment(split)))
        e3 <- (y - predict(eq3_fit, newdata=rsample::assessment(split)))
        e4 <- (y - predict(eq4_fit, newdata=rsample::assessment(split)))
        e5 <- (y - predict(eq5_fit, newdata=rsample::assessment(split)))
        e6 <- (y - predict(eq6_fit, newdata=rsample::assessment(split)))
        e7 <- (y - predict(eq7_fit, newdata=rsample::assessment(split)))
        e8 <- (y - predict(eq8_fit, newdata=rsample::assessment(split)))
        ## return the cross-validated residuals from both models as 
        ## a data frame. 
        data.frame(
          e1 = e1,
          e2 = e2,
          e3 = e3,
          e4 = e4,
          e5 = e5,
          e6 = e6,
          e7 = e7,
          e8 = e8
        )
      }
      
      # Just models w/ c0 value
      cv_mods_c0 <- function(split, ...){
        ### fit models with the analysis partition
        eq2_fit <- nlsLM(c ~ eq2(q, a, n=rec_n, c0=c0_set), data=rsample::analysis(split), start=list(a=380))
        eq5_fit <- nlsLM(c ~ eq5(q, h, g, n=rec_n, c0=c0_set), data=rsample::analysis(split), start=list(h=160, g=0.00001))
        eq8_fit <- nlsLM(c ~ eq8(q, s, b, n=rec_n, c0=c0_set), data=rsample::analysis(split),start=list(s=175, b=0.00001))
        ### take the dependent variable from the assessment partition
        y <- rsample::assessment(split)$c
        ### create the residuals  getting model predictions from the assessment partition
        e2 <- (y - predict(eq2_fit, newdata=rsample::assessment(split)))
        e5 <- (y - predict(eq5_fit, newdata=rsample::assessment(split)))
        e8 <- (y - predict(eq8_fit, newdata=rsample::assessment(split)))
        ## return the cross-validated residuals from both models as 
        ## a data frame. 
        data.frame(
          e2 = e2,
          e5 = e5,
          e8 = e8
        )
      }    
  
    ### The vfold_cv function sets up the cross-validation partitions

    ### Create list of datasets/dataset lists
    df_list_cv1 <- list(df_grab %>% mutate(df = "lt_grab"),
                        df_15min %>% mutate(df = "15min"))
    
    df_list_cv2 <- list(df_15min_daily_list,
                        df_15min_weekly_list,
                        df_15min_monthly_list)
                    
  
    # USE THIS FOR df_grab AND df_15min (i.e., the dfs where we use all the data and don't do 10 versions of each)
      # First time with one c0 value
      # Set c0
      c0_set <- c0_bf
      
      out1_c0_1 <- tibble()
      for(i in 1:length(df_list_cv1)) {
        
        # Run cross-validation  
        doParallel::registerDoParallel()
        out <- vfold_cv(df_list_cv1[[i]],
                        v = 10,
                        repeats = 10) %>% 
          ## estimate the cv on all of the partitions
          mutate(err = map(splits,
                           cv_mods)) %>% 
          ## drop the splits column (THIS WAS NECESSARY FOR THE LARGER DATASETS; OTHERWISE CODE GOT STUCK HERE)
          select(-splits)
        
        # Unnest err column in out (I separated this function for larger datasets that get hung up in previous step)
        out_unnest <-
          out %>% 
          unnest(err)
        
        # Get the average cross-validation error for each model: 
        out2 <-
          out_unnest %>%
          group_by(id) %>%
          ## sum up the sums of squared errors across partitions
          summarise(across(e1:e8, ~sd(.x))) %>%
          ## calculate average CV error:
          # summarise(across(e1:e2, ~mean(.x)))
          pivot_longer(cols = e1:e8, names_to = "eq", values_to = "sse") %>%
          group_by(eq) %>%
          summarize(rmse_mean = mean(sse),
                    rmse_n = n(),
                    rmse_sd = sd(sse),
                    rmse_se = rmse_sd/sqrt(rmse_n)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se,
                 ci_lower = rmse_mean - margin_err,
                 ci_upper = rmse_mean + margin_err) %>% 
          mutate(df = as.character(df_list_cv1[[i]][1, "df"]),
                 c0_source = "baseflow",
                 c0_value = c0_set)
        
          # Aggregate the outputs
          out1_c0_1 <- bind_rows(out1_c0_1, out2)
      }
      
      # Second time with second c0 value
      # Set c0
      c0_set <- c0_wd
      
      out1_c0_2 <- tibble()
      for(i in 1:length(df_list_cv1)) {
        
        # Run cross-validation  
        doParallel::registerDoParallel()
        out <- vfold_cv(df_list_cv1[[i]],
                        v = 10,
                        repeats = 10) %>% 
          ## estimate the cv on all of the partitions
          mutate(err = map(splits,
                           cv_mods_c0)) %>% 
          ## drop the splits column (THIS WAS NECESSARY FOR THE LARGER DATASETS; OTHERWISE CODE GOT STUCK HERE)
          select(-splits)
        
        # Unnest err column in out (I separated this function for larger datasets that get hung up in previous step)
        out_unnest <-
          out %>% 
          unnest(err)
        
        # Get the average cross-validation error for each model: 
        out2 <-
          out_unnest %>%
          group_by(id) %>%
          ## sum up the sums of squared errors across partitions
          summarise(across(e2:e8, ~sd(.x))) %>%
          ## calculate average CV error:
          # summarise(across(e1:e2, ~mean(.x)))
          pivot_longer(cols = e2:e8, names_to = "eq", values_to = "sse") %>%
          group_by(eq) %>%
          summarize(rmse_mean = mean(sse),
                    rmse_n = n(),
                    rmse_sd = sd(sse),
                    rmse_se = rmse_sd/sqrt(rmse_n)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se,
                 ci_lower = rmse_mean - margin_err,
                 ci_upper = rmse_mean + margin_err) %>% 
          mutate(df = as.character(df_list_cv1[[i]][1, "df"]),
                 c0_source = "wet_deposition",
                 c0_value = c0_set)
        
        # Aggregate the outputs
        out1_c0_2 <- bind_rows(out1_c0_2, out2)
      }      
  
  
  
    # USE THIS FOR THE REDUCED REPS OF THE 15-min DATA (e.g., df_15min_daily, etc.)
      # First time with one c0 value
      # Set c0
      c0_set <- c0_bf
      
      out2_c0_1 <- tibble()
      for(i in 1:length(df_list_cv2)) {  
        
        # Loop through the lists of dataframes    
        out_sub <- tibble()
        for (j in 1:length(df_list_cv2[[i]])) {
          out <- vfold_cv(df_list_cv2[[i]][[j]],
                          v = 10,
                          repeats = 10) %>% 
            ## estimate the cv on all of the partitions
            mutate(err = map(splits,
                             cv_mods)) %>% 
            ## drop the splits column (THIS WAS NECESSARY FOR THE LARGER DATASETS; OTHERWISE CODE GOT STUCK HERE)
            select(-splits)
          
          # Unnest err column in out (I separated this function for larger datasets that get hung up in previous step)
          out_unnest <-
            out %>% 
            unnest(err) %>% 
            mutate(df_rep = j)
          
          out_sub <- bind_rows(out_sub, out_unnest)
        }
        
        # get the average cross-validation error for each model: 
        out2 <-
          out_sub %>%
          group_by(df_rep, id) %>%
          ## sum up the sums of squared errors across partitions
          summarise(across(e1:e8, ~sd(.x))) %>%
          ## calculate average CV error:
          # summarise(across(e1:e2, ~mean(.x)))
          pivot_longer(cols = e1:e8, names_to = "eq", values_to = "sse") %>%
          group_by(df_rep, eq) %>%
          summarize(rmse_mean = mean(sse),
                    rmse_n = n(),
                    rmse_sd = sd(sse),
                    rmse_se = rmse_sd/sqrt(rmse_n)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se,
                 ci_lower = rmse_mean - margin_err,
                 ci_upper = rmse_mean + margin_err) %>%
          group_by(eq) %>%
          summarize(rmse_mean_2 = mean(rmse_mean),
                    rmse_n_2 = n(),
                    rmse_sd_2 = sd(rmse_mean),
                    rmse_se_2 = rmse_sd_2/sqrt(rmse_n_2)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n_2 - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se_2,
                 ci_lower = rmse_mean_2 - margin_err,
                 ci_upper = rmse_mean_2 + margin_err) %>%
          rename(rmse_mean = rmse_mean_2,
                 rmse_n = rmse_n_2,
                 rmse_sd = rmse_sd_2,
                 rmse_se = rmse_se_2) %>%
          mutate(method = "reps_separate") %>% 
          mutate(df = as.character(df_list_cv2[[i]][[1]][1, "df"]),
                 c0_source = "baseflow",
                 c0_value = c0_set)
        
        # Aggreate outputs
        out2_c0_1 <- bind_rows(out2_c0_1, out2)
      }
      
      # Second time with second c0 value
      # Set c0
      c0_set <- c0_wd
      
      out2_c0_2 <- tibble()
      for(i in 1:length(df_list_cv2)) {  
        
        # Loop through the lists of dataframes    
        out_sub <- tibble()
        for (j in 1:length(df_list_cv2[[i]])) {
          out <- vfold_cv(df_list_cv2[[i]][[j]],
                          v = 10,
                          repeats = 10) %>% 
            ## estimate the cv on all of the partitions
            mutate(err = map(splits,
                             cv_mods_c0)) %>% 
            ## drop the splits column (THIS WAS NECESSARY FOR THE LARGER DATASETS; OTHERWISE CODE GOT STUCK HERE)
            select(-splits)
          
          # Unnest err column in out (I separated this function for larger datasets that get hung up in previous step)
          out_unnest <-
            out %>% 
            unnest(err) %>% 
            mutate(df_rep = j)
          
          out_sub <- bind_rows(out_sub, out_unnest)
        }
        
        # get the average cross-validation error for each model: 
        out2 <-
          out_sub %>%
          group_by(df_rep, id) %>%
          ## sum up the sums of squared errors across partitions
          summarise(across(e2:e8, ~sd(.x))) %>%
          ## calculate average CV error:
          # summarise(across(e1:e2, ~mean(.x)))
          pivot_longer(cols = e2:e8, names_to = "eq", values_to = "sse") %>%
          group_by(df_rep, eq) %>%
          summarize(rmse_mean = mean(sse),
                    rmse_n = n(),
                    rmse_sd = sd(sse),
                    rmse_se = rmse_sd/sqrt(rmse_n)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se,
                 ci_lower = rmse_mean - margin_err,
                 ci_upper = rmse_mean + margin_err) %>%
          group_by(eq) %>%
          summarize(rmse_mean_2 = mean(rmse_mean),
                    rmse_n_2 = n(),
                    rmse_sd_2 = sd(rmse_mean),
                    rmse_se_2 = rmse_sd_2/sqrt(rmse_n_2)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n_2 - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se_2,
                 ci_lower = rmse_mean_2 - margin_err,
                 ci_upper = rmse_mean_2 + margin_err) %>%
          rename(rmse_mean = rmse_mean_2,
                 rmse_n = rmse_n_2,
                 rmse_sd = rmse_sd_2,
                 rmse_se = rmse_se_2) %>%
          mutate(method = "reps_separate") %>% 
          mutate(df = as.character(df_list_cv2[[i]][[1]][1, "df"]),
                 c0_source = "wet_deposition",
                 c0_value = c0_set)
        
        # Aggreate outputs
        out2_c0_2 <- bind_rows(out2_c0_2, out2)
      }      

      
# END code that requires repeating for each dataset; you can move on from here if you've run all the datasets you want ----
      
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
  