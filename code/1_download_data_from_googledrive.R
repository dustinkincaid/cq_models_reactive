# Download the aggregated data created in 0_prep_cq_data.R

if(!require(googledrive)) install.packages("googledrive")
library(googledrive)

# Shared link to file (with embedded file ID)
# https://drive.google.com/file/d/1aJcTg1f_9ks7oyxUItZBa9iHXss3eUob/view?usp=sharing

# Put googledrive into a de-authorized state
drive_deauth()

# Pull embedded file ID from 
public_file <-  drive_get(as_id("1aJcTg1f_9ks7oyxUItZBa9iHXss3eUob"))
drive_download(file = public_file, 
               path = "data/nh_cq_data_filtered.csv",
               overwrite = TRUE)
