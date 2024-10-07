#'******************************************************************************
#' Description: Clean and format spring phneology data and Daymet temperature 
#' data for flux sites. Data were downloaded by 
#' `src/dat_dl_pheno_for_flux_sites.R`.
#'******************************************************************************
rm(list=ls())

source("src/base.R")
library(data.table)
library(lubridate)


# in: 
# MCD12Q2 at flux sites
mcd_flux_file_1 <- file.path(gdir, 
    "Data/Fluxsites/Pheno-Daymet-for-Fluxsites-MCD12Q2-006-results.csv"
)
mcd_flux_file_2 <- file.path(gdir, 
    "Data/Fluxsites/Pheno-Daymet-for-Fluxsites-2-MCD12Q2-006-results.csv"
)
# Daymet at flux sites
dm_flux_file_1 <- file.path(gdir,
    "Data/Fluxsites/Pheno-Daymet-for-Fluxsites-DAYMET-004-results.csv"
)
dm_flux_file_2 <- file.path(gdir,
    "Data/Fluxsites/Pheno-Daymet-for-Fluxsites-2-DAYMET-004-results.csv"
)


# out: fluxsite_mcd12q2.csv
mcd_flux_ard_file <- file.path(gdir, "Data/ARD/fluxsite_mcd12q2.csv")
# out: fluxsite_daymet.csv
dm_flux_ard_file <- file.path(gdir, "Data/ARD/fluxsite_daymet.csv")
# out: flux_mcd12q2_daymet_dt.csv
flux_mcd_dm_file <- file.path(gdir,
    "Data/ARD/flux_mcd12q2_daymet_dt.csv"
)


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~ Data clean to ARD ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCD12Q2 at flux sites
flux_q2 <- rbind(
    fread(mcd_flux_file_1),
    fread(mcd_flux_file_2)
)
flux_q2 <- flux_q2[, .(siteID = ID, Latitude, Longitude, IGBP = Category, 
    Year = year(Date), Midgup = MCD12Q2_006_MidGreenup_0, 
    MidgupQA = MCD12Q2_006_QA_Detailed_0_MidGreenup, 
    MidgupQA_dec = MCD12Q2_006_QA_Detailed_0_MidGreenup_Description)]

# MCD12Q2 for flux sites
fwrite(flux_q2, file = mcd_flux_ard_file)


# Daymet at flux sites
flux_dm <- rbind(
    fread(dm_flux_file_1),
    fread(dm_flux_file_2)
)
flux_dm <- flux_dm[, .(siteID = ID, Latitude, Longitude, IGBP = Category, Date, 
    Tmax = DAYMET_004_tmax, Tmin = DAYMET_004_tmin)]
flux_dm[, Tmean := (Tmax + Tmin) / 2] # Calculate daily mean temperature

# Daymet for flux sites
fwrite(flux_dm, file = dm_flux_ard_file)


# ~ Format before modeling ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Phenology obvervations from MCD12Q2
pheno <- fread(mcd_flux_ard_file)

# Temperature observations from Daymet
temp <- fread(dm_flux_ard_file)
temp[, Date := as_date(Date)]


# ~ Format data
# ~~~~~~~~~~~~~~~~~
# convert to DOY
pheno[, SOS := yday(as_date(Midgup))]
# only use `Best` observations
pheno <- pheno[MidgupQA_dec == "Best", ]

# Since Daymet only covers North America
pheno <- pheno[siteID %in% unique(temp$siteID), ]

# NOTE: I wantted to format data like `phenor` does, but I figured that it's a 
# bit hard to view the result data as it's a list rather than a data table. 
# So, I format my data as a normal data table.

# For each site, retrieve temperature from this year's DOY 255 day to 365 days 
# back. This way we don't need to deal with leap years. 
# This can be changed later.
sy_temp <- apply(pheno[, .(siteID, Year, SOS)], 1, function(sy) {
    # NOTE: The Daymet calendar is based on a standard calendar year.
    # All Daymet years, including leap years, have 1 - 365 days.
    # For leap years, the Daymet database includes leap day (February 29) and 
    # values for December 31 are discarded from leap years to maintain a 
    # 365-day year.

    # R uses a 0 based index for dates
    end_date <- as_date(255 - 1, origin = paste0(sy[["Year"]], "-01-01")) 
    start_date <- end_date - 364
    if (leap_year(as_date(paste0(as.numeric(sy[["Year"]]) - 1, "-01-01")))) {
        start_date <- start_date - 1
    }
    sy_temp <- temp[
        siteID == sy[["siteID"]] & between(Date, start_date, end_date), 
        .(siteID, Date, Tmax, Tmin, Tmean)
    ]
    sy_temp$PhenoYear <- sy[["Year"]]
    sy_temp$SOS <- sy[["SOS"]]
    return(sy_temp)
})
flux_mcd12q2_daymet_dt <- do.call(rbind, sy_temp)


# MCD12Q2 and Daymet for flux sites
fwrite(flux_mcd12q2_daymet_dt, file = flux_mcd_dm_file)
