#'******************************************************************************
#' Description: Linear regress SOS by preseason monthly mean temperature to act 
#' as a bench mark for spring phenology model comparison.
#'******************************************************************************

# The data in this linear regression analysis includes:
#   - MCD12Q2 + Daymet at flux sites
#   - MCD12Q2 + flux measured temperature
#   - HF & HB pheno + ground temperature
#   - PEP725 + E-OBS

# Follow Basler et al 2016, the linear regression models conducted here will try 
# different single monthly mean temperatures from Jan to May, and some 
# combinations of the mean temperatures of several months, the best linear 
# regression result determined by RMSE will be returned.

source("src/base.R")

library(data.table)
library(lubridate)



LmFit <- function(dt) {
    # Convert to site anomalies
    Convert2Anomaly <- function(month_mean) {
        month_mean[, Tmean := scale(Tmean, scale = FALSE), by = "siteID"]
        month_mean$SOS <- as.double(month_mean$SOS)
        month_mean[, SOS := scale(SOS, scale = FALSE), by = "siteID"]
        
        return(month_mean)
    }

    # Try single month means
    mod_fits <- list()
    for (i in 1:5) {
        month_mean <- unique(
            dt[
                month(Date) == i, 
                .(Tmean = mean(Tmean), SOS = SOS),
                by = .(siteID, PhenoYear)
            ]
        )
        
        month_mean <- Convert2Anomaly(month_mean)
        
        # plot(month_mean[siteID == "US-Ha1", .(SOS, Tmean)])
        # plot(month_mean[, .(SOS, Tmean)], 
        #     col = factor(month_mean$siteID), 
        #     pch = 16
        # )
        lm_fit <- lm(SOS ~ Tmean, data = month_mean)
        # summary(lm_fit)
        rmse <- sqrt(mean(lm_fit$residuals^2))

        mod_fits <- append(mod_fits, list( 
            list(mod = lm_fit, rmse = rmse)))
    }

    # Try Jan-March mean
    month_mean <- unique(dt[month(Date) %in% c(1, 2, 3), 
        .(Tmean = mean(Tmean), SOS = SOS),
        by = .(siteID, PhenoYear)
    ])
    month_mean <- Convert2Anomaly(month_mean)

    # plot(month_mean[, .(SOS, Tmean)])
    lm_fit <- lm(SOS ~ Tmean, data = month_mean)
    # summary(lm_fit)
    rmse <- sqrt(mean(lm_fit$residuals^2))

    mod_fits <- append(mod_fits, list(
        list(mod = lm_fit, rmse = rmse)
    ))

    # Try March-May mean
    month_mean <- unique(dt[month(Date) %in% c(3, 4, 5),
        .(Tmean = mean(Tmean), SOS = SOS),
        by = .(siteID, PhenoYear)
    ])
    month_mean <- Convert2Anomaly(month_mean)

    # plot(month_mean[, .(SOS, Tmean)])
    lm_fit <- lm(SOS ~ Tmean, data = month_mean)
    # summary(lm_fit)
    rmse <- sqrt(mean(lm_fit$residuals^2))

    mod_fits <- append(mod_fits, list(
        list(mod = lm_fit, rmse = rmse)
    ))


    # Get the best single month fit
    best_mod_fit <- list(mod = NULL, rmse = 9999)
    for (i in seq_along(mod_fits)) {
        rmse <- mod_fits[[i]]$rmse
        if (best_mod_fit$rmse > rmse) {
            best_mod_fit$mod <- mod_fits[[i]]$mod
            best_mod_fit$rmse <- rmse
        }
    }

    return(best_mod_fit)
}



# ~ MCD12Q2 + Daymet at flux sites ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read data in
flux_dt <- fread(file.path(gdir, "Data/ARD/flux_mcd12q2_daymet_dt.csv"))
flux_dt[, Date := as_date(Date)]
# merge site meta
site_meta <- fread(file.path(gdir, "Data/ARD/sites_fluxnet2015_Tier1.csv"))
flux_dt <- merge(flux_dt, site_meta[, .(siteID, IGBP)], by = "siteID")
# filter out SOS > 255
flux_dt <- flux_dt[SOS < 255 & IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), ]

# Pooled linear regression
mod_mcd_dm_flux_lm_pooled <- LmFit(flux_dt)
# out: mod_mcd_dm_flux_lm_pooled.Rds
saveRDS(mod_mcd_dm_flux_lm_pooled, file.path(gdir, "Pipeline", 
    "mod_mcd_dm_flux_lm_pooled.Rds"))

# Site-specific linear regression
mod_mcd_dm_flux_lm_site <- by(flux_dt, flux_dt$siteID, LmFit)
# out: mod_mcd_dm_flux_lm_site.Rds
saveRDS(mod_mcd_dm_flux_lm_site, file.path(gdir, "Pipeline", 
    "mod_mcd_dm_flux_lm_site.Rds"))




# ~ MCD12Q2 + flux_temperature ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#! Some sites have too few site-years, so no site-specific linear regression for
#!   this combination

# in: flux temperature
flux_temp <- fread(file.path(
    gdir,
    "Data/ARD/flux_measured_daily_temperature.csv"
))
flux_temp[, Date := as_date(Date)]

# in: flux mcd12q2 phenology
flux_dt <- fread(file.path(gdir, "Data/ARD/flux_mcd12q2_daymet_dt.csv"))
flux_dt[, Date := as_date(Date)]
# merge site meta
site_meta <- fread(file.path(gdir, "Data/ARD/sites_fluxnet2015_Tier1.csv"))
flux_dt <- merge(flux_dt, site_meta[, .(siteID, IGBP)], by = "siteID")
# filter out SOS > 255
flux_dt <- flux_dt[SOS < 255 & IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), ]


# Merge the two
flux_dt <- merge(flux_dt, flux_temp, by = c("siteID", "Date"), all.x = TRUE)

# Filter out site-years with measured temperature data that have missing values
flux_temp_dt <- by(flux_dt, flux_dt[, .(siteID, PhenoYear)], function(sy) {
    if (any(is.na(sy$Temp)) == FALSE) {
        return(sy)
    }
})
flux_temp_dt <- do.call(rbind, flux_temp_dt)
# flux_temp_dt[, .N / 365, by = siteID]

# Filter out non-forest sites
flux_temp_dt <- flux_temp_dt[IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), .(
    siteID, Date,
    Tmax = NA, Tmin = NA, Tmean = Temp, SOS, PhenoYear
)]



# Pooled linear regression
mod_mcd_flux_temp_lm_pooled <- LmFit(flux_temp_dt)
# out: mod_mcd_flux_temp_lm_pooled.Rds
saveRDS(mod_mcd_flux_temp_lm_pooled, file.path(gdir, "Pipeline", 
    "mod_mcd_flux_temp_lm_pooled.Rds"))



# ~ MCD12Q2 + Daymet for site years that also have flux temperature ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flux_dm_temp_dt <- by(flux_dt, flux_dt[, .(siteID, PhenoYear)], function(sy) {
    if (any(is.na(sy$Temp)) == FALSE) {
        return(sy)
    }
})
flux_dm_temp_dt <- do.call(rbind, flux_dm_temp_dt)

# Pooled linear regression
mod_mcd_dm_flux_temp_pooled <- LmFit(flux_dm_temp_dt)
# out: mod_mcd_dm_flux_temp_pooled.Rds
saveRDS(mod_mcd_dm_flux_temp_pooled, file.path(gdir, "Pipeline", 
    "mod_mcd_dm_flux_temp_pooled.Rds"))


# ~ HF & HB ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HF
hf_pheno_dt <- fread(file.path(gdir, "Data/ARD/hf_grd_temperature_pheno.csv"))
hf_pheno_dt[, Date := as_date(Date)]

# HB
hb_pheno_dt <- fread(file.path(gdir, "Data/ARD/hb_grd_temperature_pheno.csv"))
hb_pheno_dt[, Date := as_date(Date)]

hb_hf_dt <- rbind(
    hf_pheno_dt[, .(siteID, Date,
        Tmean = Temp, Tmax = NA, Tmin = NA, SOS,
        PhenoYear
    )],
    hb_pheno_dt[, .(siteID, Date, Tmean, Tmax, Tmin, SOS, PhenoYear)]
)


mod_hb_hf_lm_pooled <- LmFit(hb_hf_dt)
# out: mod_hb_hf_lm_pooled.Rds
saveRDS(mod_hb_hf_lm_pooled, file.path(gdir, "Pipeline", 
    "mod_hb_hf_lm_pooled.Rds"))

mod_hb_hf_lm_site <- by(hb_hf_dt, hb_hf_dt$siteID, LmFit)
# out: mod_hb_hf_lm_site.Rds
saveRDS(mod_hb_hf_lm_site, file.path(gdir, "Pipeline", 
    "mod_hb_hf_lm_site.Rds"))



# ~ PEP725 + E-OBS ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# There are so many sites, so use parallel processing
library(parallel)

pep725_species <- c("aesculus", "alnus", "betula", "fagus", 
    "fraxinus", "quercus"
)

# species_name <- "aesculus" # debug
for (species_name in pep725_species) {
    # Read site-specific caches
    caches <- list.files(
        file.path(hpc_local_dir,"Pipeline/tmp"),
        paste0(species_name, ".*.Rds$"),
        full.names = TRUE
    )

    # Make cluster
    cl <- makeCluster(detectCores() - 1)
    calls <- clusterCall(cl, function() {
        wd <- getwd()
        setwd(wd)

        suppressWarnings({
            source("Code/mod_compare_base.R")
        })
    })
    clusterExport(cl, c("LmFit"))

    # x <- caches[1] # debug
    mod_pep725_site <- clusterApply(cl = cl, caches, function(x) {
        cur_site <- readRDS(x)
        lm_fit <- LmFit(cur_site)
        return(lm_fit)
    })
    # Get site names
    mod_pep725_site_names <- sapply(caches, function(x) {
        siteID <- strsplit(tools::file_path_sans_ext(basename(x)), "_")[[1]][2]
        return(siteID)
    })
    names(mod_pep725_site) <- mod_pep725_site_names

    # out: mod_pep725_lm_site_*.Rds
    saveRDS(mod_pep725_site, file.path(gdir, "Pipeline", paste0(
        "mod_pep725_lm_site_", species_name, ".Rds"
    )))

    stopCluster(cl)
}


