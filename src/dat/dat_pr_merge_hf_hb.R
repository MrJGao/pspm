#'******************************************************************************
#' Description: Merge HF and HB phenological and temperature data together.
#'******************************************************************************

library(data.table)
library(lubridate)

source("src/base.R")


# in:
# HF data
hf_pheno_dt_file <- "data/hf_grd_temperature_pheno.csv"
# HB data
hb_pheno_dt_file <- "data/hb_grd_temperature_pheno.csv"

# out: hb_hf_pheno_temp.csv
hb_hf_pheno_temp_file <- "data/hb_hf_pheno_temp.csv"

# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# HF
hf_pheno_dt <- fread(hf_pheno_dt_file)
hf_pheno_dt[, Date := as_date(Date)]

# HB
hb_pheno_dt <- fread(hb_pheno_dt_file)
hb_pheno_dt[, Date := as_date(Date)]

# Combine them together
hb_hf_dt <- rbind(
    hf_pheno_dt[, .(siteID, Date,
        Tmean = Temp, Tmax = NA, Tmin = NA, SOS,
        PhenoYear
    )],
    hb_pheno_dt[, .(siteID, Date, Tmean, Tmax, Tmin, SOS, PhenoYear)]
)

# Export the combined data table
fwrite(hb_hf_dt, hb_hf_pheno_temp_file)


