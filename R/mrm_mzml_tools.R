library(MSnbase)
library(glue)
library(tidyverse)
library(readr)
library(here)

# Functions to read and convert MRM chromatograms from mzML files
# ---------------------------------------------------------------


# Bo Burla, Singapore Lipidomics Incubator
# 2021


# Save chromatograms extracted from mzML files as CSV files: 
# - One CSV file per MRM feature
# - Each CSV file: Column 1: Retention time, Following columns: intensities of samples
convert_mzml_to_csv <- function(path_mzml, path_output){

  file <- dir(path = path_mzml,
              full.name = TRUE,
              pattern = "mzML$")
  
  mrm <- readSRMData(file)
  
  file_names <- str_replace(string = fs::path_file(file), pattern = ".mzML","")
  feature_count <- nrow(mrm)
  
  for (i in seq(feature_count)) {
    
    df_rt <- enframe( x = rtime(mrm[i,1]), value =  "RT", name = "Timepoint")
    df_rt <- df_rt[,-1]
    
    feat <- mrm[i,]     
    df_int <- map_dfc(feat, "intensity")
  
    df <- bind_cols(df_rt, df_int)
    colnames(df) <- c("RT", file_names)
    
    feature_name <- glue("Feat-{formatC(i, width=3, flag='0')}_Q1_{round(precursorMz(feat),1)[1]}-Q3_{round(productMz(feat)[1],1)}")
  
    write_csv(x = df, file = paste0(path_output, "/", feature_name, ".csv"))
    
    return(mrm)
  }
}


