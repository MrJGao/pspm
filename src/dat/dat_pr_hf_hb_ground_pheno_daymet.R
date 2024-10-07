#'******************************************************************************
#' Description: Clean ground phenological observations at Harvard Forest and 
#' Hubbard Brook Experimental Forest as well as the corresponding Daymet 
#' temperatures.
#'******************************************************************************
rm(list=ls())

source("src/base.R")
library(raster)
library(data.table)
library(lubridate)


# in: Daymet from AppEEARS
dm_file_1 <- file.path(
    gdir, "Data",
    "HF_HB_daymet/Daymet-HF-HB-1-DAYMET-004-results.csv"
)
dm_file_2 <- file.path(
    gdir, "Data",
    "HF_HB_daymet/Daymet-HF-HB-2-DAYMET-004-results.csv"
)
# HF ground phenology from Bayesian LSP project
hf_pheno_file <- paste0(
    "Q:/My Drive/Research/Landsat_LSP/Data/",
    "HF/Harvard_Forest_groud/HF_ground.Rds"
)
# HB ground phenology from Bayesian LSP project
hb_pheno_file <- "Q:/My Drive/Research/Landsat_LSP/Pipeline/phn_phenos.Rds"

# out: hf_hb_pheno_daymet.csv
grd_pheno_file <- file.path(gdir, "Data/ARD/hf_hb_pheno_daymet.csv")


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Daymet data
dm <- rbind(fread(dm_file_1), fread(dm_file_2))
dm <- dm[, .(
    siteID = ID, Lat = Latitude, Lon = Longitude,
    Date = as_date(Date), Tmax = DAYMET_004_tmax, Tmin = DAYMET_004_tmin
)]

# ~ Harvard Forest
hf_pheno <- readRDS(hf_pheno_file)
hf_pheno <- hf_pheno[, 
    .(siteID = "US-Ha1", SOS = round(spring), PhenoYear = year)
]


# ~ Hubbard Brook
hb_pheno <- readRDS(hb_pheno_file)
# Only need midgup
hb_pheno <- hb_pheno[, .(siteID = site, SOS = yday(midgup), 
    PhenoYear = year(midgup))]
hb_pheno <- na.omit(hb_pheno)
# Correct siteID name
hb_pheno$siteID <- gsub("X", "", hb_pheno$siteID)


# ~ Clean and format
pheno <- rbind(hf_pheno, hb_pheno)
# For each pheno year, find the corresponding Daymet temperatures
grd_pheno_dt <- data.table()
for (i in 1:nrow(pheno)) {
    cur_sy <- pheno[i, ]
    # NOTE: The Daymet calendar is based on a standard calendar year.
    # All Daymet years, including leap years, have 1 - 365 days.
    # For leap years, the Daymet database includes leap day (February 29) and 
    # values for December 31 are discarded from leap years to maintain a 365-day 
    # per year.
    
    # R uses a 0 based index for dates
    end_date <- as_date(
        255 - 1, 
        origin = paste0(cur_sy[["PhenoYear"]], "-01-01")
    ) 
    start_date <- end_date - 364
    if (leap_year(as_date(paste0(as.numeric(cur_sy[["PhenoYear"]]) - 1, 
        "-01-01")))) {
        start_date <- start_date - 1
    }
    cur_dm <- dm[
        siteID == cur_sy$siteID & between(Date, start_date, end_date), 
    ]
    cur_dm$SOS <- cur_sy$SOS
    cur_dm$PhenoYear <- cur_sy$PhenoYear

    grd_pheno_dt <- rbind(grd_pheno_dt, cur_dm)
}

# Calculate daily mean temperature
grd_pheno_dt[, Tmean := (Tmax + Tmin) / 2]


# ground phenology daymet dt
fwrite(grd_pheno_dt, file = grd_pheno_file)





