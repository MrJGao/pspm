#'******************************************************************************
#' Description: Model comparion for MCD12Q2 + flux_temperature at flux sites. 
#' B/c flux measured temperature data for most sites have limited site years, 
#' this combination can only fit pooled model.
#'******************************************************************************
source("src/mod_compare_base.R")

outdir <- "pipe/goodness-of-fit"
if (!dir.exists(outdir)) {
    dir.create(outdir)
}

# ~ Pooled ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # in: flux temperature
# flux_temp <- fread("data/flux_measured_daily_temperature.csv")
# flux_temp[, Date := as_date(Date)]

# # in: flux mcd12q2 phenology
# flux_dt <- fread("data/flux_mcd12q2_daymet_dt.csv")
# flux_dt[, Date := as_date(Date)]
# # merge site meta
# site_meta <- fread("data/raw/sites_fluxnet2015_Tier1.csv")
# flux_dt <- merge(flux_dt, site_meta[, .(siteID, IGBP)], by = "siteID")
# # filter out SOS > 255
# flux_dt <- flux_dt[SOS < 255 & IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), ]


# # Merge the two
# flux_dt <- merge(flux_dt, flux_temp, by = c("siteID", "Date"), all.x = TRUE)

# # Filter out site-years with measured temperature data that have missing values
# flux_temp_dt <- by(flux_dt, flux_dt[, .(siteID, PhenoYear)], function(sy) {
#     if (any(is.na(sy$Temp)) == FALSE) {
#         return(sy)
#     }
# })
# flux_temp_dt <- do.call(rbind, flux_temp_dt)
# flux_temp_dt[, .N / 365, by = siteID]

# # Filter out non-forest sites
# flux_temp_dt <- flux_temp_dt[IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), .(
#     siteID, Date, Tmax = NA, Tmin = NA, Tmean = Temp, SOS, PhenoYear
# )]

# print("Fit MCD12Q2 + Flux_temperature....")
# flux_temp_li <- FormatDataForPhenoModel(flux_temp_dt)

# # out: save for later use
# saveRDS(flux_temp_li, file.path("pipe", "flux_temp_li.Rds"))


flux_temp_li <- readRDS(file.path("pipe", "flux_temp_li.Rds"))


# Fit pooled models
mod_mcd_flux_temp_pooled <- FitCompareModels(flux_temp_li)

# out: model fit for MCD12Q2 + Daymet at flux sites
saveRDS(mod_mcd_flux_temp_pooled, file.path(
    "pipe", "goodness-of-fit", "mcd_flux_temp_pool",
    "mod_mcd_flux_temp_pooled.Rds"
))