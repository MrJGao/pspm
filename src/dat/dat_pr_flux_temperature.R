#'******************************************************************************
#' Clean downloaded FLUXNET2015 data to get temperature records for phenology
#' modeling. The data were downloaded directly from the FLUXNET2015 website. The
#' downloaded raw data are stored in HPC under `Urban_pheno/Data/FLUXNET2015`.
#'******************************************************************************

library(data.table)
library(lubridate)

source("src/base.R")

# in: 
# fluxnest2015
data_dir <- "data/FLUXNET2015/"
zip_data_files <- list.files(
    data_dir, 
    pattern = "FLX.*FLUXNET2015_SUBSET.*.zip$"
)


# out: flux_measured_daily_temperature.csv
flux_temp_file <-  "data/flux_measured_daily_temperature.csv"


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt <- lapply(zip_data_files, function(zf) {
    # just for debugging:
    # zipfile <- file.path(data_dir, 
    #     "FLX_US-NR1_FLUXNET2015_SUBSET_1998-2014_1-4.zip")
    
    zipfile <- file.path(data_dir, zf)
    ziplist <- unzip(zipfile, list = TRUE) # get the file list
    target_file <- ziplist[grep(".*._SUBSET_DD_.*.", ziplist$Name), "Name"]
    
    # get site ID
    siteID <- strsplit(target_file, "_")[[1]][2]
    
    # read data and subset temperatre related columns only
    dt <- fread(cmd = paste("unzip -p", zipfile, target_file))
    dt <- dt[, .(siteID = siteID, TIMESTAMP, TA_F, TA_F_QC)]
    return(dt)
})
dt <- do.call(rbind, dt)

# Convert to Date format
dt[, Date := as_date(as.character(TIMESTAMP))]

dt <- dt[, .(siteID, Date, Temp = TA_F, Temp_QC = TA_F_QC)]


# Write out flux site-measured daily temperature
fwrite(dt, file = flux_temp_file)
