# Fit the CQ data to Hall's models & do cross-validation

# Load libraries ----
  library("tidyverse")    # general workhorse
  library("here")         # takes the guesswork out of dealing with file paths
  library("data.table")   # read in large datasets faster
  library("reshape2")     # manipulate dataframes
  library("rsample")      # for model cross-validation
  library("doParallel")   # for faster computing
  library("minpack.lm")   # alternative to 'nls'; much more flexible!
  
# Which site and solute would you like to model?
  site_choice <- "HBF"
  sol_choice <- "doc"
  grab_choice <- "doc_grab"

# Options ----
  set.seed(1217)
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
  
  rm(df_sub)
  
# Reduce sensor data to the various sampling times (ie, daily, weekly, monthly) ----
  # Here we create 20 versions of each sampling frequency (e.g., df_15min_daily_1 through df_15min_daily_10)  
  
  # Daily (sample always occurs sometime between 8am and 5pm)
  df_15min_daily_list <- list()
  for(i in 1:30){
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
  for(i in 1:30){
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
  for(i in 1:30){
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
  
  # Write these data to CSV for future use
  df_grab %>% 
    mutate(df = "lt_grab") %>% 
    bind_rows(df_15min %>% 
                mutate(df = "15min"),
              df_15min_daily_list[[1]] %>% 
                mutate(df = "daily"),
              df_15min_weekly_list[[1]] %>% 
                mutate(df = "weekly"),
              df_15min_monthly_list[[1]] %>% 
                mutate(df = "monthly")) %>% 
    write_csv(here("data", "out", paste0("cq_all_datasets", "_", sol_choice, "_", site_choice, ".csv")))  
  
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
                 param_b = eq0_fit$coefficients[[2]],
                 df = df_list[[i]]$df)
        # Here I recalculate the residual using back-transformed predicted values 
        resids_eq0 <- df_list[[i]]$c - 10^(predict(eq0_fit))
        print(sprintf("RMSE Eq 0: %s", RMSE(resids_eq0)))
        # print(sprintf("RMSE Eq 0: %s", RMSE(summary(eq0_fit)$residuals)))    
        
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
          eq4_nFit_fit <- nlsLM(c ~ eq4(q, h, g, n), data=df_list[[i]], start=list(h=175, g=0.00001, n=0.5))
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
          eq6_nFit_fit <- nlsLM(c ~ eq6(q, m, j, n), data=df_list[[i]], start=list(m=100, j=100, n = 0.5))
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
    
    # Write to CSV for future use
    resids_final %>% 
      write_csv(here("data", "out",  paste0("cq_all_resids", "_", "baseflow", "_", sol_choice, "_", site_choice, ".csv")))
    bdf_all %>% 
      write_csv(here("data", "out", paste0("cq_all_bdf", "_", "baseflow", "_", sol_choice, "_", site_choice, ".csv")))
    plot_df_all %>% 
      write_csv(here("data", "out",  paste0("cq_all_fits", "_", "baseflow", "_", sol_choice, "_", site_choice, ".csv")))    
    
    # Remove objects
    rm(eq0_fit, eq1_fit, eq1_nFit_fit, eq2_fit, eq2_nFit_fit, eq3_fit, eq4_fit,
       eq4_nFit_fit, eq5_fit, eq5_nFit_fit, eq6_fit, eq6_nFit_fit, eq7_fit,
       eq7_nFit_fit, eq8_fit, eq8_nFit_fit,
       eq0_resids, eq1_resids, eq1_nFit_resids, eq2_resids, eq2_nFit_resids, eq3_resids, eq4_resids,
       eq4_nFit_resids, eq5_resids, eq5_nFit_resids, eq6_resids, eq6_nFit_resids, eq7_resids,
       eq7_nFit_resids, eq8_resids, eq8_nFit_resids)    
  
    
    ### REPEAT where c0 = c0_wd ----
    ### Only need to re-run models with C0 as input
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
      plotdf = data.frame(df_list[[i]]$q, 
                          predict(eq2_fit), predict(eq2_nFit_fit), 
                          predict(eq5_fit), predict(eq5_nFit_fit), 
                          predict(eq8_fit), predict(eq8_nFit_fit))
      colnames(plotdf) = c('q','eq2','eq2_nFit','eq5','eq5_nFit','eq8','eq8_nFit')
      
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
    
    # Write to CSV for future use
    resids_final_2 %>% 
      write_csv(here("data", "out",  paste0("cq_all_resids", "_", "wetdep", "_", sol_choice, "_", site_choice, ".csv")))
    bdf_all_2 %>% 
      write_csv(here("data", "out", paste0("cq_all_bdf", "_", "wetdep", "_", sol_choice, "_", site_choice, ".csv")))
    plot_df_all_2 %>% 
      write_csv(here("data", "out",  paste0("cq_all_fits", "_", "wetdep", "_", sol_choice, "_", site_choice, ".csv")))    
    
    # Remove objects
    rm(eq2_fit, eq2_nFit_fit, eq5_fit, eq5_nFit_fit, eq8_fit, eq8_nFit_fit,
       eq2_resids, eq2_nFit_resids, eq5_resids, eq5_nFit_resids, eq8_resids, eq8_nFit_resids)      
  
  
# Evaluate model performance using repeated K/V-fold cross-validation ----    
  ## Create function to fit and assess models ----
    ### All models ----
    cv_mods <- function(split, ...){
      ### fit models with the analysis partition
      eq0_fit <- lm(log10(c) ~ log10(q), data=rsample::analysis(split))
      eq1_fit <- nlsLM(c ~ eq1(q, a, n=rec_n), data=rsample::analysis(split), start=list(a=380))
      # eq1_n_fit <- nlsLM(c ~ eq1(q, a, n), data=rsample::analysis(split), start=list(a=380, n=0.1))
      eq1_n_fit <- nlsLM(c ~ eq1(q, a, n), data=rsample::analysis(split), start=list(a=380, n=1))
      eq2_fit <- nlsLM(c ~ eq2(q, a, n=rec_n, c0=c0_set), data=rsample::analysis(split), start=list(a=380))
      eq2_n_fit <- nlsLM(c ~ eq2(q, a, n, c0=c0_set), data=rsample::analysis(split), start=list(a=380, n=1))
      eq3_fit <- nlsLM(c ~ eq3(q, f, e), data=rsample::analysis(split), start=list(f=100, e=100))
      eq4_fit <- nlsLM(c ~ eq4(q, h, g, n=rec_n), data=rsample::analysis(split), start=list(h=175, g=0.00001))
      eq4_n_fit <- nlsLM(c ~ eq4(q, h, g, n), data=rsample::analysis(split), start=list(h=175, g=0.00001, n=0.5))
      eq5_fit <- nlsLM(c ~ eq5(q, h, g, n=rec_n, c0=c0_set), data=rsample::analysis(split), start=list(h=160, g=0.00001))
      eq5_n_fit <- nlsLM(c ~ eq5(q, h, g, n, c0=c0_set), data=rsample::analysis(split), start=list(h=160, g=0.00001, n=1))
      eq6_fit <- nlsLM(c ~ eq6(q, m, j, n=rec_n), data=rsample::analysis(split), start=list(m=100, j=100))
      eq6_n_fit <- nlsLM(c ~ eq6(q, m, j, n), data=rsample::analysis(split), start=list(m=100, j=100, n = 0.5))
      eq7_fit <- nlsLM(c ~ eq7(q, s, b, n=rec_n), data=rsample::analysis(split), start=list(s=175, b=0.00001))
      eq7_n_fit <- nlsLM(c ~ eq7(q, s, b, n), data=rsample::analysis(split), start=list(s=175, b=0.00001, n = 0.1))
      eq8_fit <- nlsLM(c ~ eq8(q, s, b, n=rec_n, c0=c0_set), data=rsample::analysis(split), start=list(s=175, b=0.00001))
      eq8_n_fit <- nlsLM(c ~ eq8(q, s, b, n, c0=c0_set), data=rsample::analysis(split), start=list(s=175, b=0.00001, n=1))
      ### take the dependent variable from the assessment partition
      y <- rsample::assessment(split)$c
      ### create the residuals  getting model predictions from the assessment partition
      e0 <- (y - 10^(predict(eq0_fit, newdata=rsample::assessment(split))))
      e1 <- (y - predict(eq1_fit, newdata=rsample::assessment(split)))
      e1_nFit <- (y - predict(eq1_n_fit, newdata=rsample::assessment(split)))
      e2 <- (y - predict(eq2_fit, newdata=rsample::assessment(split)))
      e2_nFit <- (y - predict(eq2_n_fit, newdata=rsample::assessment(split)))
      e3 <- (y - predict(eq3_fit, newdata=rsample::assessment(split)))
      e4 <- (y - predict(eq4_fit, newdata=rsample::assessment(split)))
      e4_nFit <- (y - predict(eq4_n_fit, newdata=rsample::assessment(split)))
      e5 <- (y - predict(eq5_fit, newdata=rsample::assessment(split)))
      e5_nFit <- (y - predict(eq5_n_fit, newdata=rsample::assessment(split)))
      e6 <- (y - predict(eq6_fit, newdata=rsample::assessment(split)))
      e6_nFit <- (y - predict(eq6_n_fit, newdata=rsample::assessment(split)))
      e7 <- (y - predict(eq7_fit, newdata=rsample::assessment(split)))
      e7_nFit <- (y - predict(eq7_n_fit, newdata=rsample::assessment(split)))
      e8 <- (y - predict(eq8_fit, newdata=rsample::assessment(split)))
      e8_nFit <- (y - predict(eq8_n_fit, newdata=rsample::assessment(split)))
      ## return the cross-validated residuals from both models as 
      ## a data frame. 
      data.frame(
        e0 = e0,
        e1 = e1,
        e1_nFit = e1_nFit,
        e2 = e2,
        e2_nFit = e2_nFit,
        e3 = e3,
        e4 = e4,
        e4_nFit = e4_nFit,
        e5 = e5,
        e5_nFit = e5_nFit,
        e6 = e6,
        e6_nFit = e6_nFit,
        e7 = e7,
        e7_nFit = e7_nFit,
        e8 = e8,
        e8_nFit = e8_nFit
      )
    }
    
    ### Define possibly model
    possibly_cv_mods = possibly(.f = cv_mods, otherwise = NULL)
    
    ### Just models w/ c0 value ----
    cv_mods_c0 <- function(split, ...){
      ### fit models with the analysis partition
      eq2_fit <- nlsLM(c ~ eq2(q, a, n=rec_n, c0=c0_set), data=rsample::analysis(split), start=list(a=380))
      eq2_n_fit <- nlsLM(c ~ eq2(q, a, n, c0=c0_set), data=rsample::analysis(split), start=list(a=380, n=1))
      eq5_fit <- nlsLM(c ~ eq5(q, h, g, n=rec_n, c0=c0_set), data=rsample::analysis(split), start=list(h=160, g=0.00001))
      eq5_n_fit <- nlsLM(c ~ eq5(q, h, g, n, c0=c0_set), data=rsample::analysis(split), start=list(h=160, g=0.00001, n=1))
      eq8_fit <- nlsLM(c ~ eq8(q, s, b, n=rec_n, c0=c0_set), data=rsample::analysis(split),start=list(s=175, b=0.00001))
      eq8_n_fit <- nlsLM(c ~ eq8(q, s, b, n, c0=c0_set), data=rsample::analysis(split), start=list(s=175, b=0.00001, n=1))
      ### take the dependent variable from the assessment partition
      y <- rsample::assessment(split)$c
      ### create the residuals  getting model predictions from the assessment partition
      e2 <- (y - predict(eq2_fit, newdata=rsample::assessment(split)))
      e2_nFit <- (y - predict(eq2_n_fit, newdata=rsample::assessment(split)))
      e5 <- (y - predict(eq5_fit, newdata=rsample::assessment(split)))
      e5_nFit <- (y - predict(eq5_n_fit, newdata=rsample::assessment(split)))
      e8 <- (y - predict(eq8_fit, newdata=rsample::assessment(split)))
      e8_nFit <- (y - predict(eq8_n_fit, newdata=rsample::assessment(split)))
      ## return the cross-validated residuals from both models as 
      ## a data frame. 
      data.frame(
        e2 = e2,
        e2_nFit = e2_nFit,
        e5 = e5,
        e5_nFit = e5_nFit,
        e8 = e8,
        e8_nFit = e8_nFit
      )
    }
    
    ### Define possibly model
    possibly_cv_mods_c0 = possibly(.f = cv_mods_c0, otherwise = NULL)
  
  ## Create list of datasets/dataset lists ----
  df_list_cv1 <- list(df_grab %>% mutate(df = "lt_grab"),
                      df_15min %>% mutate(df = "15min"))
  
  df_list_cv2 <- list(df_15min_daily_list,
                      df_15min_weekly_list,
                      df_15min_monthly_list)
 
  ## Do cross-validation ----   
    # This method assumes that the time series data are sorted by timestamp ahead of time
    # Helpful websites:
      # https://www.tidymodels.org/learn/models/time-series/
      # https://rsample.tidymodels.org/reference/rolling_origin.html
      # https://www.tmwr.org/resampling#rolling
      # https://openforecast.org/adam/rollingOrigin.html
    
    ### DF_GRAB AND DF_15MIN ----
    # USE THIS FOR df_grab AND df_15min (i.e., the dfs where we use all the data and don't do 10 versions of each)
    
      #### First time with one c0 value ----
      # Set c0
      c0_set <- c0_bf
      
      
      out1_c0_1 <- tibble()
      for(i in 1:length(df_list_cv1)) {
        
        # Run cross-validation  
        doParallel::registerDoParallel()
        
        # Split smaller grab sample df slightly differently than 15-min data set
        if(df_list_cv1[[i]]$df[[1]] == "lt_grab"){
          
          # Set sample sizes of analysis and assessment splits; 
          # We will make the minimum sample size for analysis = 20
          initial_n = ceiling(nrow(df_list_cv1[[i]])/2)
          if(initial_n < 20){
            intial_n = 20
          }
          assess_n = ceiling(initial_n/3)
          skip_n = 1
          
          # Create splits
          splits <- rsample::rolling_origin(df_list_cv1[[i]],
                                            initial = initial_n,
                                            assess = assess_n,
                                            skip = skip_n)
          
          # Pull min and max for each split for NRMSE later
          maxmin_c <-
            splits %>% 
            mutate(assessments = map(splits, assessment)) %>%
            unnest(assessments) %>%
            summarise(min_c = min(c), 
                      max_c = max(c), 
                      .by = id)
        
        # Split larger 15-min dataset   
        } else {
          
          # 1 day = 24 * 60 / 15 = 96 rows
          # 14 days = 96 * 14 = 1344
          # 30 days = 96 * 30 = 2880
          # 3 months-ish or 300 days = 2880 * 3
          initial_n = ceiling(nrow(df_list_cv1[[i]])/2)
          assess_n = ceiling(initial_n/3)
          skip_n = 2880
          
          splits <- rsample::rolling_origin(df_list_cv1[[i]],
                                            initial = initial_n,
                                            assess = assess_n,
                                            skip = skip_n)
          
          # Pull min and max for each split for NRMSE later
          maxmin_c <-
            splits %>% 
            mutate(assessments = map(splits, assessment)) %>%
            unnest(assessments) %>%
            summarise(min_c = min(c), 
                      max_c = max(c), 
                      .by = id)
        }
        
        # Estimate the cv on all splits
        out <-
          splits %>% 
          mutate(err = map(splits,
                           possibly_cv_mods)) %>% 
          # Remove models with errors (typically singular gradient errors)
          compact() %>% 
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
          summarise(across(e0:e8_nFit, ~sd(.x))) %>%
          ## calculate average CV error:
          # summarise(across(e1:e2, ~mean(.x)))
          pivot_longer(cols = e0:e8_nFit, names_to = "eq", values_to = "sse") %>%
          full_join(maxmin_c, by = join_by(id)) %>% 
          mutate(n_sse = sse/(max_c - min_c)) %>% 
          group_by(eq) %>%
          summarize(rmse_mean = mean(sse),
                    rmse_n = n(),
                    rmse_sd = sd(sse),
                    rmse_se = rmse_sd/sqrt(rmse_n),
                    nrmse_mean = mean(n_sse),
                    nrmse_sd = sd(n_sse),
                    nrmse_se = nrmse_sd/sqrt(rmse_n)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se,
                 rmse_ci_lower = rmse_mean - margin_err,
                 rmse_ci_upper = rmse_mean + margin_err,
                 nrmse_ci_lower = nrmse_mean - margin_err,
                 nrmse_ci_upper = nrmse_mean + margin_err) %>% 
          mutate(df = as.character(df_list_cv1[[i]][1, "df"]),
                 c0_source = "baseflow",
                 c0_value = c0_set)
        
        # Aggregate the outputs
        out1_c0_1 <- bind_rows(out1_c0_1, out2)
      }
      
      #### Second time with second c0 value ----
      # Set c0
      c0_set <- c0_wd
      
      out1_c0_2 <- tibble()
      for(i in 1:length(df_list_cv1)) {
        
        # Run cross-validation  
        doParallel::registerDoParallel()
        
        # Split smaller grab sample df slightly differently than 15-min data set
        if(df_list_cv1[[i]]$df[[1]] == "lt_grab"){
          
          # Set sample sizes of analysis and assessment splits; 
          # We will make the minimum sample size for analysis = 20
          initial_n = ceiling(nrow(df_list_cv1[[i]])/2)
          if(initial_n < 20){
            intial_n = 20
          }
          assess_n = ceiling(initial_n/3)
          skip_n = 1
          
          # Create splits
          splits <- rsample::rolling_origin(df_list_cv1[[i]],
                                            initial = initial_n,
                                            assess = assess_n,
                                            skip = skip_n)
          
          # Pull min and max for each split for NRMSE later
          maxmin_c <-
            splits %>% 
            mutate(assessments = map(splits, assessment)) %>%
            unnest(assessments) %>%
            summarise(min_c = min(c), 
                      max_c = max(c), 
                      .by = id)
          
          # Split larger 15-min dataset   
        } else {
          
          # 1 day = 24 * 60 / 15 = 96 rows
          # 14 days = 96 * 14 = 1344
          # 30 days = 96 * 30 = 2880
          # 3 months-ish or 300 days = 2880 * 3
          initial_n = ceiling(nrow(df_list_cv1[[i]])/2)
          assess_n = ceiling(initial_n/3)
          skip_n = 2880
          
          splits <- rsample::rolling_origin(df_list_cv1[[i]],
                                            initial = initial_n,
                                            assess = assess_n,
                                            skip = skip_n)
          
          # Pull min and max for each split for NRMSE later
          maxmin_c <-
            splits %>% 
            mutate(assessments = map(splits, assessment)) %>%
            unnest(assessments) %>%
            summarise(min_c = min(c), 
                      max_c = max(c), 
                      .by = id)
        }
        
        # Estimate the cv on all splits
        out <-
          splits %>% 
          mutate(err = map(splits,
                           possibly_cv_mods_c0)) %>% 
          # Remove models with errors (typically singular gradient errors)
          compact() %>% 
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
          summarise(across(e2:e8_nFit, ~sd(.x))) %>%
          ## calculate average CV error:
          # summarise(across(e1:e2, ~mean(.x)))
          pivot_longer(cols = e2:e8_nFit, names_to = "eq", values_to = "sse") %>%
          full_join(maxmin_c, by = join_by(id)) %>% 
          mutate(n_sse = sse/(max_c - min_c)) %>% 
          group_by(eq) %>%
          summarize(rmse_mean = mean(sse),
                    rmse_n = n(),
                    rmse_sd = sd(sse),
                    rmse_se = rmse_sd/sqrt(rmse_n),
                    nrmse_mean = mean(n_sse),
                    nrmse_sd = sd(n_sse),
                    nrmse_se = nrmse_sd/sqrt(rmse_n)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se,
                 rmse_ci_lower = rmse_mean - margin_err,
                 rmse_ci_upper = rmse_mean + margin_err,
                 nrmse_ci_lower = nrmse_mean - margin_err,
                 nrmse_ci_upper = nrmse_mean + margin_err) %>% 
          mutate(df = as.character(df_list_cv1[[i]][1, "df"]),
                 c0_source = "wet_deposition",
                 c0_value = c0_set)
        
        # Aggregate the outputs
        out1_c0_2 <- bind_rows(out1_c0_2, out2)
      }            
    
    ### REDUCED REPS OF THE 15-min DATA ----
    # USE THIS FOR THE REDUCED REPS OF THE 15-min DATA (e.g., df_15min_daily, etc.)
      
      #### First time with one c0 value ----
      # Set c0
      c0_set <- c0_bf
      
      out2_c0_1 <- tibble()
      for(i in 1:length(df_list_cv2)) {  
  
        # Run cross-validation  
        doParallel::registerDoParallel()              
        
        # Loop through the lists of dataframes    
        out_sub <- tibble()
        maxmin_c_sub <- tibble()
        for (j in 1:length(df_list_cv2[[i]])) {
          
          # # Split larger daily sample df slightly differently than smaller dataset
          if(df_list_cv2[[i]][[j]]$df[[1]] == "daily") {
            
            # Set sample sizes of analysis and assessment splits; 
            # We will make the minimum sample size for analysis = 20
            initial_n = ceiling(nrow(df_list_cv2[[i]][[j]])/2)
            assess_n = ceiling(initial_n/3)
            skip_n = 30
            
            splits <- rsample::rolling_origin(df_list_cv2[[i]][[j]],
                                              initial = initial_n,
                                              assess = assess_n,
                                              skip = skip_n)
            
            # Pull min and max for each split for NRMSE later
            maxmin_c <-
              splits %>% 
              mutate(assessments = map(splits, assessment)) %>%
              unnest(assessments) %>%
              summarise(min_c = min(c), 
                        max_c = max(c), 
                        .by = id) %>% 
              mutate(df_rep = j)
          
          # Split weekly datasets  
          } else if(df_list_cv2[[i]][[j]]$df[[1]] == "weekly"){
            
            initial_n = ceiling(nrow(df_list_cv2[[i]][[j]])/2)
            assess_n = ceiling(initial_n/3)
            skip_n = 4
            
            splits <- rsample::rolling_origin(df_list_cv2[[i]][[j]],
                                              initial = initial_n,
                                              assess = assess_n,
                                              skip = skip_n)
            
            # Pull min and max for each split for NRMSE later
            maxmin_c <-
              splits %>% 
              mutate(assessments = map(splits, assessment)) %>%
              unnest(assessments) %>%
              summarise(min_c = min(c), 
                        max_c = max(c), 
                        .by = id) %>% 
              mutate(df_rep = j) 
            
          # Split monthly datasets  
          } else {
            
            initial_n = ceiling(nrow(df_list_cv2[[i]][[j]])/2)
            if(initial_n < 20){
              intial_n = 20
            }
            assess_n = ceiling(initial_n/3)
            skip_n = 0
            
            splits <- rsample::rolling_origin(df_list_cv2[[i]][[j]],
                                              initial = initial_n,
                                              assess = assess_n,
                                              skip = skip_n)
            
            # Pull min and max for each split for NRMSE later
            maxmin_c <-
              splits %>% 
              mutate(assessments = map(splits, assessment)) %>%
              unnest(assessments) %>%
              summarise(min_c = min(c), 
                        max_c = max(c), 
                        .by = id) %>% 
              mutate(df_rep = j)             
          }
          
          out <- 
            splits %>% 
            ## estimate the cv on all of the partitions
            mutate(err = map(splits,
                             possibly_cv_mods)) %>% 
            # Remove models with errors (typically singular gradient errors)
            compact() %>% 
            ## drop the splits column (THIS WAS NECESSARY FOR THE LARGER DATASETS; OTHERWISE CODE GOT STUCK HERE)
            select(-splits)
          
          # Unnest err column in out (I separated this function for larger datasets that get hung up in previous step)
          out_unnest <-
            out %>% 
            unnest(err) %>% 
            mutate(df_rep = j)
          
          out_sub <- bind_rows(out_sub, out_unnest)
          maxmin_c_sub <- bind_rows(maxmin_c_sub, maxmin_c)
        }
        
        # get the average cross-validation error for each model: 
        out2 <-
          out_sub %>%
          group_by(df_rep, id) %>%
          ## sum up the sums of squared errors across partitions
          summarise(across(e0:e8_nFit, ~sd(.x))) %>%
          ## calculate average CV error:
          # summarise(across(e1:e2, ~mean(.x)))
          pivot_longer(cols = e0:e8_nFit, names_to = "eq", values_to = "sse") %>%
          full_join(maxmin_c_sub, by = join_by(df_rep, id)) %>% 
          mutate(n_sse = sse/(max_c - min_c)) %>%
          group_by(df_rep, eq) %>%
          summarize(rmse_mean = mean(sse),
                    rmse_n = n(),
                    rmse_sd = sd(sse),
                    rmse_se = rmse_sd/sqrt(rmse_n),
                    nrmse_mean = mean(n_sse),
                    nrmse_sd = sd(n_sse),
                    nrmse_se = nrmse_sd/sqrt(rmse_n)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se,
                 rmse_ci_lower = rmse_mean - margin_err,
                 rmse_ci_upper = rmse_mean + margin_err,
                 nrmse_ci_lower = nrmse_mean - margin_err,
                 nrmse_ci_upper = nrmse_mean + margin_err) %>%
          group_by(eq) %>%
          summarize(rmse_mean_2 = mean(rmse_mean),
                    rmse_n_2 = n(),
                    rmse_sd_2 = sd(rmse_mean),
                    rmse_se_2 =  rmse_sd_2/sqrt(rmse_n_2),
                    nrmse_mean_2 = mean(nrmse_mean),
                    nrmse_sd_2 = sd(nrmse_mean),
                    nrmse_se_2 = nrmse_sd_2/sqrt(rmse_n_2)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n_2 - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se_2,
                 rmse_ci_lower = rmse_mean_2 - margin_err,
                 rmse_ci_upper = rmse_mean_2 + margin_err,
                 nrmse_ci_lower = nrmse_mean_2 - margin_err,
                 nrmse_ci_upper = nrmse_mean_2 + margin_err) %>%
          rename(rmse_mean = rmse_mean_2,
                 rmse_n = rmse_n_2,
                 rmse_sd = rmse_sd_2,
                 rmse_se = rmse_se_2,
                 nrmse_mean = nrmse_mean_2,
                 nrmse_sd = nrmse_sd_2,
                 nrmse_se = nrmse_se_2) %>%
          mutate(method = "reps_separate") %>% 
          mutate(df = as.character(df_list_cv2[[i]][[1]][1, "df"]),
                 c0_source = "baseflow",
                 c0_value = c0_set)
        
        # Aggreate outputs
        out2_c0_1 <- bind_rows(out2_c0_1, out2)
      }
      
      #### Second time with second c0 value ----
      # Set c0
      c0_set <- c0_wd
      
      out2_c0_2 <- tibble()
      for(i in 1:length(df_list_cv2)) {  
        
        # Loop through the lists of dataframes    
        out_sub <- tibble()
        maxmin_c_sub <- tibble()
        for (j in 1:length(df_list_cv2[[i]])) {
          
          # Split larger daily sample df slightly differently than smaller dataset
          if(df_list_cv2[[i]][[j]]$df[[1]] == "daily") {
            
            # Set sample sizes of analysis and assessment splits; 
            # We will make the minimum sample size for analysis = 20
            initial_n = ceiling(nrow(df_list_cv2[[i]][[j]])/2)
            assess_n = ceiling(initial_n/3)
            skip_n = 30
            
            splits <- rsample::rolling_origin(df_list_cv2[[i]][[j]],
                                              initial = initial_n,
                                              assess = assess_n,
                                              skip = skip_n)
            
            # Pull min and max for each split for NRMSE later
            maxmin_c <-
              splits %>% 
              mutate(assessments = map(splits, assessment)) %>%
              unnest(assessments) %>%
              summarise(min_c = min(c), 
                        max_c = max(c), 
                        .by = id) %>% 
              mutate(df_rep = j)
            
            # Split weekly datasets  
          } else if(df_list_cv2[[i]][[j]]$df[[1]] == "weekly"){
            
            initial_n = ceiling(nrow(df_list_cv2[[i]][[j]])/2)
            assess_n = ceiling(initial_n/3)
            skip_n = 4
            
            splits <- rsample::rolling_origin(df_list_cv2[[i]][[j]],
                                              initial = initial_n,
                                              assess = assess_n,
                                              skip = skip_n)
            
            # Pull min and max for each split for NRMSE later
            maxmin_c <-
              splits %>% 
              mutate(assessments = map(splits, assessment)) %>%
              unnest(assessments) %>%
              summarise(min_c = min(c), 
                        max_c = max(c), 
                        .by = id) %>% 
              mutate(df_rep = j) 
            
            # Split monthly datasets  
          } else {
            
            initial_n = ceiling(nrow(df_list_cv2[[i]][[j]])/2)
            if(initial_n < 20){
              intial_n = 20
            }
            assess_n = ceiling(initial_n/3)
            skip_n = 0
            
            splits <- rsample::rolling_origin(df_list_cv2[[i]][[j]],
                                              initial = initial_n,
                                              assess = assess_n,
                                              skip = skip_n)
            
            # Pull min and max for each split for NRMSE later
            maxmin_c <-
              splits %>% 
              mutate(assessments = map(splits, assessment)) %>%
              unnest(assessments) %>%
              summarise(min_c = min(c), 
                        max_c = max(c), 
                        .by = id) %>% 
              mutate(df_rep = j)             
          }
          
          out <- 
            splits %>% 
            ## estimate the cv on all of the partitions
            mutate(err = map(splits,
                             possibly_cv_mods_c0)) %>% 
            # Remove models with errors (typically singular gradient errors)
            compact() %>% 
            ## drop the splits column (THIS WAS NECESSARY FOR THE LARGER DATASETS; OTHERWISE CODE GOT STUCK HERE)
            select(-splits)
          
          # Unnest err column in out (I separated this function for larger datasets that get hung up in previous step)
          out_unnest <-
            out %>% 
            unnest(err) %>% 
            mutate(df_rep = j)
          
          out_sub <- bind_rows(out_sub, out_unnest)
          maxmin_c_sub <- bind_rows(maxmin_c_sub, maxmin_c)
        }
        
        # get the average cross-validation error for each model: 
        out2 <-
          out_sub %>%
          group_by(df_rep, id) %>%
          ## sum up the sums of squared errors across partitions
          summarise(across(e2:e8_nFit, ~sd(.x))) %>%
          ## calculate average CV error:
          # summarise(across(e1:e2, ~mean(.x)))
          pivot_longer(cols = e2:e8_nFit, names_to = "eq", values_to = "sse") %>%
          full_join(maxmin_c_sub, by = join_by(df_rep, id)) %>% 
          mutate(n_sse = sse/(max_c - min_c)) %>%
          group_by(df_rep, eq) %>%
          summarize(rmse_mean = mean(sse),
                    rmse_n = n(),
                    rmse_sd = sd(sse),
                    rmse_se = rmse_sd/sqrt(rmse_n),
                    nrmse_mean = mean(n_sse),
                    nrmse_sd = sd(n_sse),
                    nrmse_se = nrmse_sd/sqrt(rmse_n)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se,
                 rmse_ci_lower = rmse_mean - margin_err,
                 rmse_ci_upper = rmse_mean + margin_err,
                 nrmse_ci_lower = nrmse_mean - margin_err,
                 nrmse_ci_upper = nrmse_mean + margin_err) %>%
          group_by(eq) %>%
          summarize(rmse_mean_2 = mean(rmse_mean),
                    rmse_n_2 = n(),
                    rmse_sd_2 = sd(rmse_mean),
                    rmse_se_2 =  rmse_sd_2/sqrt(rmse_n_2),
                    nrmse_mean_2 = mean(nrmse_mean),
                    nrmse_sd_2 = sd(nrmse_mean),
                    nrmse_se_2 = nrmse_sd_2/sqrt(rmse_n_2)) %>%
          mutate(alpha = 0.05,
                 deg_fr = rmse_n_2 - 1,
                 t_score = qt(p = alpha/2, df = deg_fr, lower.tail = F),
                 margin_err = t_score * rmse_se_2,
                 rmse_ci_lower = rmse_mean_2 - margin_err,
                 rmse_ci_upper = rmse_mean_2 + margin_err,
                 nrmse_ci_lower = nrmse_mean_2 - margin_err,
                 nrmse_ci_upper = nrmse_mean_2 + margin_err) %>%
          rename(rmse_mean = rmse_mean_2,
                 rmse_n = rmse_n_2,
                 rmse_sd = rmse_sd_2,
                 rmse_se = rmse_se_2,
                 nrmse_mean = nrmse_mean_2,
                 nrmse_sd = nrmse_sd_2,
                 nrmse_se = nrmse_se_2) %>%
          mutate(method = "reps_separate") %>% 
          mutate(df = as.character(df_list_cv2[[i]][[1]][1, "df"]),
                 c0_source = "wet_deposition",
                 c0_value = c0_set)
        
        # Aggreate outputs
        out2_c0_2 <- bind_rows(out2_c0_2, out2)
      }      

      #### Write to CSV for future use ----
      out1_c0_1 %>%
        bind_rows(out1_c0_2,
                  out2_c0_1,
                  out2_c0_2) %>% 
        write_csv(here("data", "out",  paste0("cq_all_cv", "_", sol_choice, "_", site_choice, ".csv")))    
      
