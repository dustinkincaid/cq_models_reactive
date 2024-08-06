# Workflow for cq_models_reactive

#### Download aggregated data from Google Drive folder
Use **1_download_data_from_googledrive.R** to download the aggregated C-Q data (nh_cq_data_filtered.csv) for all 
NH watersheds that is assembled using **0_prep_cq_data.R** and used in most subsequent scripts. **0_prep_cq_data.R**
requires raw data files that have not been uploaded.

#### Estimate C0 from wet deposition data and baseflow concentrations
If you need to re-estimate C0 for either wet deposition or baseflow, you can use **2_summarize_wet_deposition.R**
or **3_estimate_c0_from_baseflow_concs.R**, respectively. Otherwise, you'll use the output files from these, i.e., 
wet_dep_mean_concs.csv and c0_baseflow_estimates.csv, respectively.

#### Create figures of raw data
Use **4_create_figure_comparing_grabs_to_sensor.R** and **5_create_figure_logC_vs_logQ.R** to plot grab samples vs 
sensor data and log(C)~log(Q) plots, respectively.

#### Model C-Q data with power-law
Use **6_model_cq_power_law_XX_XXX.R** to model C-Q data with the power-law using nonlinear least squares regression.
