#'******************************************************************************
#' Model comparion for MCD12Q2 + Daymet at flux sites but using 
#' site years that have flux measured temperature so to test whether flux 
#' measured temperature increased the model fit over Daymet temperature.
#'******************************************************************************

# Since the number of site-years data are limited for flux measured temperature, 
# do pooled analysis only.

source("src/mod_compare_base.R")


# ~ Pooled ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read data in
flux_dt <- fread("data/flux_mcd12q2_daymet_dt.csv")
flux_dt[, Date := as_date(Date)]
# merge site meta
site_meta <- fread("data/raw/sites_fluxnet2015_Tier1.csv")
flux_dt <- merge(flux_dt, site_meta[, .(siteID, IGBP)], by = "siteID")
# filter out SOS > 255
flux_dt <- flux_dt[SOS < 255 & IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), ]

# flux_dt[, .N / 365, siteID]


# in: flux temperature
flux_temp <- fread("data/flux_measured_daily_temperature.csv")
flux_temp[, Date := as_date(Date)]


# Merge the two
flux_dt <- merge(flux_dt, flux_temp, by = c("siteID", "Date"), all.x = TRUE)

# Filter out site-years with measured temperature data that have missing values
flux_dm_temp_dt <- by(flux_dt, flux_dt[, .(siteID, PhenoYear)], function(sy) {
    if (any(is.na(sy$Temp)) == FALSE) {
        return(sy)
    }
})
flux_dm_temp_dt <- do.call(rbind, flux_dm_temp_dt)
# flux_dm_temp_dt[, .N / 365, by = siteID]

# Filter out non-forest sites but use Daymet temperature
flux_dm_temp_dt <- flux_dm_temp_dt[
    IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), 
    .(siteID, Date, Tmax = NA, Tmin = NA, Tmean, SOS, PhenoYear)
]



print("Fit MCD12Q2 + Daymet pooled for sites that have flux temperature....")

flux_dm_temp_li <- FormatDataForPhenoModel(flux_dm_temp_dt)
# out: save flux_dm_temp_li for later use
saveRDS(flux_dm_temp_li, "pipe/flux_dm_temp_li.Rds")

# Read the saved Rds
flux_dm_temp_li <- readRDS("pipe/flux_dm_temp_li.Rds")

# Fit pooled models
mod_mcd_dm_flux_temp_pooled <- FitCompareModels(flux_dm_temp_li)
# out: model fit for MCD12Q2 + Daymet at flux sites
saveRDS(mod_mcd_dm_flux_temp_pooled, 
    "pipe/goodness-of-fit/mod_mcd_dm_flux-temp-sites_pooled.Rds"
)

