#'******************************************************************************
#' Cross-validation for linear regression on multiple data sources.
#'******************************************************************************

# The data in this linear regression analysis includes:
#   - MCD12Q2 + Daymet at flux sites
#   - MCD12Q2 + flux measured temperature
#   - HF & HB pheno + ground temperature
#   - PEP725 + E-OBS

# To eliminate the dependence in data from the same site, all data will be 
# converted to site anomalies before spliting them into training and testing 
# sets. Follow Basler et al 2016, the linear regression models conducted here 
# will try different single monthly mean temperatures from Jan to May, and some 
# combinations of the mean temperatures of several months, the best linear 
# regression result determined by RMSE will be used to predict the test dataset.

# Linear regression is relatively fast, so this script only uses the local 
# machine to parallelly perform the job. If the speed is too slow in the future, 
# we will use HPC to conduct the job.


source("src/base.R")

library(data.table)
library(lubridate)
library(parallel)

# Convert to site anomalies
Convert2Anomaly <- function(month_mean) {
    month_mean[, Tmean := scale(Tmean, scale = FALSE), by = "siteID"]
    month_mean$SOS <- as.double(month_mean$SOS)
    month_mean[, SOS := scale(SOS, scale = FALSE), by = "siteID"]

    return(month_mean)
}


LeaveOneSiteOutCV <- function(data_dt) {
    LmSiteCv <- function(month_mean, site) {
        # Split to training and test sets
        cv_res <- by(month_mean, month_mean$siteID, function(test_site) {
            train_dt <- month_mean[siteID != unique(test_site$siteID)]
            
            lm_fit <- lm(SOS ~ Tmean, data = train_dt)
            pred <- predict(lm_fit, newdata = test_site)

            return(data.table(site = test_site$siteID, 
                obs = test_site$SOS, pred = pred))
        })
        cv_res <- do.call(rbind, cv_res)

        return(cv_res)
    }
    
    # Try single month means
    mod_fits <- list()
    for (i in 1:5) {
        month_mean <- unique(data_dt[month(Date) == i,
            .(Tmean = mean(Tmean), SOS = SOS),
            by = .(siteID, PhenoYear)
        ])

        month_mean <- Convert2Anomaly(month_mean)

        cv_res <- LmSiteCv(month_mean)
        rmse <- sqrt(mean((cv_res$pred - cv_res$obs)^2))

        mod_fits <- append(mod_fits, list(
            list(mod = cv_res, rmse = rmse, type = paste0("month_", i))
        ))
    }

    # Try Jan-March mean
    month_mean <- unique(data_dt[month(Date) %in% c(1, 2, 3),
        .(Tmean = mean(Tmean), SOS = SOS),
        by = .(siteID, PhenoYear)
    ])
    month_mean <- Convert2Anomaly(month_mean)

    cv_res <- LmSiteCv(month_mean)
    rmse <- sqrt(mean((cv_res$pred - cv_res$obs)^2))

    mod_fits <- append(mod_fits, list(
        list(mod = cv_res, rmse = rmse, type = paste0("Jan-Mar"))
    ))

    # Try March-May mean
    month_mean <- unique(data_dt[month(Date) %in% c(3, 4, 5),
        .(Tmean = mean(Tmean), SOS = SOS),
        by = .(siteID, PhenoYear)
    ])
    month_mean <- Convert2Anomaly(month_mean)

    cv_res <- LmSiteCv(month_mean)
    rmse <- sqrt(mean((cv_res$pred - cv_res$obs)^2))

    mod_fits <- append(mod_fits, list(
        list(mod = cv_res, rmse = rmse, type = paste0("Mar-May"))
    ))


    # Get the best single month fit
    best_mod_fit <- list(mod = NULL, rmse = 9999)
    for (i in seq_along(mod_fits)) {
        rmse <- mod_fits[[i]]$rmse
        if (best_mod_fit$rmse > rmse) {
            best_mod_fit$mod <- mod_fits[[i]]$mod
            best_mod_fit$rmse <- rmse
            best_mod_fit$type <- mod_fits[[i]]$type
        }
    }

    return(best_mod_fit)
}


LeaveOneYearOutCV <- function(site_dt) {
    if (uniqueN(site_dt[, .(PhenoYear)]) < 3) {
        return(NULL)
    }
    # site_dt <- flux_dt[siteID == "US-Ha1",]
    LmYearCV <- function(month_mean) {
        cv_res <- by(month_mean, month_mean$PhenoYear, function(test_sy) {
            train_dt <- month_mean[PhenoYear != unique(test_sy$PhenoYear)]
            
            lm_fit <- lm(SOS ~ Tmean, data = train_dt)
            pred <- predict(lm_fit, newdata = test_sy)

            return(data.table(site = test_sy$siteID, 
                obs = test_sy$SOS, pred = pred))
        })
        cv_res <- do.call(rbind, cv_res)

        return(cv_res)
    }
    
    # Try single month means
    mod_fits <- list()
    for (i in 1:5) {
        month_mean <- unique(site_dt[month(Date) == i,
            .(Tmean = mean(Tmean), SOS = SOS),
            by = .(siteID, PhenoYear)
        ])

        month_mean <- Convert2Anomaly(month_mean)

        cv_res <- LmYearCV(month_mean)
        rmse <- sqrt(mean((cv_res$pred - cv_res$obs)^2))

        mod_fits <- append(mod_fits, list(
            list(mod = cv_res, rmse = rmse, type = paste0("month_", i))
        ))
    }

    # Try Jan-March mean
    month_mean <- unique(site_dt[month(Date) %in% c(1, 2, 3),
        .(Tmean = mean(Tmean), SOS = SOS),
        by = .(siteID, PhenoYear)
    ])
    month_mean <- Convert2Anomaly(month_mean)

    cv_res <- LmYearCV(month_mean)
    rmse <- sqrt(mean((cv_res$pred - cv_res$obs)^2))

    mod_fits <- append(mod_fits, list(
        list(mod = cv_res, rmse = rmse, type = paste0("Jan-Mar"))
    ))

    # Try March-May mean
    month_mean <- unique(site_dt[month(Date) %in% c(3, 4, 5),
        .(Tmean = mean(Tmean), SOS = SOS),
        by = .(siteID, PhenoYear)
    ])
    month_mean <- Convert2Anomaly(month_mean)

    cv_res <- LmYearCV(month_mean)
    rmse <- sqrt(mean((cv_res$pred - cv_res$obs)^2))

    mod_fits <- append(mod_fits, list(
        list(mod = cv_res, rmse = rmse, type = paste0("Mar-May"))
    ))


    # Get the best single month fit
    best_mod_fit <- list(mod = NULL, rmse = 9999)
    for (i in seq_along(mod_fits)) {
        rmse <- mod_fits[[i]]$rmse
        if (best_mod_fit$rmse > rmse) {
            best_mod_fit$mod <- mod_fits[[i]]$mod
            best_mod_fit$rmse <- rmse
            best_mod_fit$type <- mod_fits[[i]]$type
        }
    }

    return(best_mod_fit)
}



outdir <- "pipe/mod_cv"
if (!dir.exists(outdir)) {
    dir.create(outdir)
}


# ~ MCD12Q2 + Daymet at flux sites ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print("Do MCD12Q2 + Daymet...")

# Read data in
flux_dt <- fread("data/flux_mcd12q2_daymet_dt.csv")
flux_dt[, Date := as_date(Date)]
# merge site meta
site_meta <- fread("data/raw/sites_fluxnet2015_Tier1.csv")
flux_dt <- merge(flux_dt, site_meta[, .(siteID, IGBP)], by = "siteID")
# filter out SOS > 255
flux_dt <- flux_dt[SOS < 255 & IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), ]

# Pooled linear regression
mod_mcd_dm_flux_lm_pooled_cv <- LeaveOneSiteOutCV(flux_dt)
# out: mod_mcd_dm_flux_lm_pooled_cv.Rds
saveRDS(mod_mcd_dm_flux_lm_pooled_cv, file.path(
    outdir,
    "mod_mcd_dm_flux_lm_pooled_cv.Rds"
))


# Site-specific linear regression
mod_mcd_dm_flux_lm_site_cv <- by(flux_dt, flux_dt$siteID, LeaveOneYearOutCV)


# out: mod_mcd_dm_flux_lm_site_cv.Rds
saveRDS(mod_mcd_dm_flux_lm_site_cv, file.path(
    outdir,
    "mod_mcd_dm_flux_lm_site_cv.Rds"
))




# ~ MCD12Q2 + flux_temperature ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#! Some sites have too few site-years, so no site-specific linear regression for
#!   this combination

print("Do MCD12Q2 + flux_temperature...")

# in: flux temperature
flux_temp <- fread("data/flux_measured_daily_temperature.csv")
flux_temp[, Date := as_date(Date)]

# in: flux mcd12q2 phenology
flux_dt <- fread(file.path(gdir, "data/flux_mcd12q2_daymet_dt.csv"))
flux_dt[, Date := as_date(Date)]
# merge site meta
site_meta <- fread(file.path(gdir, "data/raw/sites_fluxnet2015_Tier1.csv"))
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
flux_temp_dt[, .N / 365, by = siteID]

# Filter out non-forest sites
flux_temp_dt <- flux_temp_dt[IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), .(
    siteID, Date,
    Tmax = NA, Tmin = NA, Tmean = Temp, SOS, PhenoYear
)]



# Pooled linear regression
mod_mcd_flux_temp_lm_pooled_cv <- LeaveOneSiteOutCV(flux_temp_dt)
# out: mod_mcd_flux_temp_lm_pooled_cv.Rds
saveRDS(mod_mcd_flux_temp_lm_pooled_cv, file.path(
    outdir, 
    "mod_mcd_flux_temp_lm_pooled_cv.Rds"
))




# ~ HF & HB ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print("Do HF & HB...")

# HF
hf_pheno_dt <- fread(file.path(gdir, "data/hf_grd_temperature_pheno.csv"))
hf_pheno_dt[, Date := as_date(Date)]

# HB
hb_pheno_dt <- fread(file.path(gdir, "data/hb_grd_temperature_pheno.csv"))
hb_pheno_dt[, Date := as_date(Date)]

hb_hf_dt <- rbind(
    hf_pheno_dt[, .(siteID, Date,
        Tmean = Temp, Tmax = NA, Tmin = NA, SOS,
        PhenoYear
    )],
    hb_pheno_dt[, .(siteID, Date, Tmean, Tmax, Tmin, SOS, PhenoYear)]
)


mod_hb_hf_lm_pooled_cv <- LeaveOneSiteOutCV(hb_hf_dt)
# out: mod_hb_hf_lm_pooled_cv.Rds
saveRDS(mod_hb_hf_lm_pooled_cv, file.path(
    outdir,
    "mod_hb_hf_lm_pooled_cv.Rds"
))

mod_hb_hf_lm_site_cv <- by(hb_hf_dt, hb_hf_dt$siteID, LeaveOneYearOutCV)
# out: mod_hb_hf_lm_site_cv.Rds
saveRDS(mod_hb_hf_lm_site_cv, file.path(
    outdir,
    "mod_hb_hf_lm_site_cv.Rds"
))



# ~ PEP725 + E-OBS ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print("Do PEP725 + E-OBS...")

# There are so many sites, so use parallel processing
library(parallel)

pep725_species <- c("aesculus", "alnus", "betula", "fagus", 
    "fraxinus", "quercus"
)

# species_name <- "aesculus" # debug
for (species_name in pep725_species) {
    # Read site-specific caches
    caches <- list.files(
        file.path("pipe/tmp"),
        paste0(species_name, ".*.Rds$"),
        full.names = TRUE
    )

    # Make cluster
    cl <- makeCluster(detectCores())
    calls <- clusterCall(cl, function() {
        suppressWarnings({
            source("src/mod_compare_base.R")
        })
    })
    clusterExport(cl, c("LeaveOneYearOutCV", "Convert2Anomaly"))

    # x <- caches[1] # debug
    mod_pep725_site <- clusterApply(cl = cl, caches, function(x) {
        cur_site <- readRDS(x)
        lm_fit <- LeaveOneYearOutCV(cur_site)
        return(lm_fit)
    })
    # Get site names
    mod_pep725_site_names <- sapply(caches, function(x) {
        siteID <- strsplit(tools::file_path_sans_ext(basename(x)), "_")[[1]][2]
        return(siteID)
    })
    names(mod_pep725_site) <- mod_pep725_site_names

    # out: mod_pep725_lm_site_*.Rds
    saveRDS(mod_pep725_site, file.path(outdir, paste0(
        "mod_pep725_lm_site_", species_name, "_cv.Rds"
    )))

    stopCluster(cl)
}


