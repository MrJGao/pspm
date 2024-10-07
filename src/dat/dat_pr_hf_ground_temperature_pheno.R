#'******************************************************************************
#' Clean and extract the Harvard Forest ground measured temperature and ground
#' observed phenology.
#'******************************************************************************
rm(list=ls())

source("src/base.R")
library(data.table)
library(lubridate)


# in: 
# Ground phenology from Bayesian LSP project
hf_pheno_file <- paste0(
    "Q:/My Drive/Research/Landsat_LSP/Data/",
    "HF/Harvard_Forest_groud/HF_ground.Rds"
)
# flux tower measured temperature
flux_temp_file <- "data/flux_measured_daily_temperature.csv"

# out: hf_grd_temperature_pheno.csv
grd_pheno_file <- "data/hf_grd_temperature_pheno.csv"


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hf_pheno <- readRDS(hf_pheno_file)
hf_pheno <- hf_pheno[, 
    .(siteID = "US-Ha1", SOS = round(spring), PhenoYear = year)
]

flux_temp <- fread(flux_temp_file)
flux_temp[, Date := as_date(Date)]
flux_temp <- flux_temp[siteID == "US-Ha1", ]


# Clean and format
grd_pheno_dt <- data.table()
for (i in 1:nrow(hf_pheno)) {
    cur_sy <- hf_pheno[i, ]
    # NOTE: I extract temperature from DOY -110 to 255
    # R uses a 0 based index for dates
    end_date <- as_date(
        255 - 1, 
        origin = paste0(cur_sy[["PhenoYear"]], "-01-01")
    )
    start_date <- end_date - 364
    if (leap_year(as_date(paste0(
        as.numeric(cur_sy[["PhenoYear"]]) - 1, "-01-01")))
    ) {
        start_date <- start_date - 1
    }
    
    cur_temp <- flux_temp[siteID == cur_sy$siteID & 
        between(Date, start_date, end_date), ]
    # Filter out site years with less than 365 days temperature
    if (nrow(cur_temp) != 365 | any(is.na(cur_temp$Temp))) {
        next
    }

    cur_temp$SOS <- cur_sy$SOS
    cur_temp$PhenoYear <- cur_sy$PhenoYear

    grd_pheno_dt <- rbind(grd_pheno_dt, cur_temp)
}


# hf ground pheno + ground temperature
fwrite(grd_pheno_dt, file = grd_pheno_file)
