#'******************************************************************************
#' Description: For pooled analysis, do external validation to test spatial 
#' extrapolation.
#'******************************************************************************
source("src/base.R")
source("src/mod_compare_base.R")


SampleHalfSitesOut <- function(data_li, num_test_sites) {
    sites <- unique(data_li$site)
    test_sites <- sample(sites, num_test_sites)
    sy_idx_test <- (1:length(data_li$site))[data_li$site %in% test_sites]
    sy_idx_train <- (1:length(data_li$site))[-sy_idx_test]
    testing_li <- GetDataLiIndex(data_li, sy_idx_test)
    training_li <- GetDataLiIndex(data_li, sy_idx_train)

    res_li <- list(
        test_sites = test_sites, 
        train = training_li, 
        test = testing_li
    )
    return(res_li)
}


set.seed(8964)

# ~ MCD12Q2 + Daymet at flux sites ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flux_li <- readRDS(file.path(gdir, "Pipeline", "flux_li.Rds"))
# Sample 50% of sites out
flux_split_li <- SampleHalfSitesOut(flux_li, 13)
ext_valid_res <- CvCompareModels(flux_split_li)

saveRDS(ext_valid_res, file.path(
    gdir, "Pipeline", "mod_mcd_dm_flux_ext_valid.Rds"
))



# ~ HF & HB ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hf_hb_li <- readRDS(file.path(gdir, "Pipeline", "hb_hf_li.Rds"))
# Sample 50% of sites out
hf_hb_split_li <- SampleHalfSitesOut(hf_hb_li, 3)
ext_valid_res <- CvCompareModels(hf_hb_split_li)

saveRDS(ext_valid_res, file.path(
    gdir, "Pipeline", "mod_hf_hb_ext_valid.Rds"
))



# ~ MCD12Q2 + flux_temperature ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flux_temp_li <- readRDS(file.path(gdir, "Pipeline", "flux_temp_li.Rds"))
# Sample 50% of sites out
flux_temp_split_li <- SampleHalfSitesOut(flux_temp_li, 11)
ext_valid_res <- CvCompareModels(flux_temp_split_li)

saveRDS(ext_valid_res, file.path(
    gdir, "Pipeline", "mod_mcd_flux_temp_ext_valid.Rds"
))



# ~ Linear regression ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert to site anomalies
Convert2Anomaly <- function(month_mean) {
    month_mean[, Tmean := scale(Tmean, scale = FALSE), by = "siteID"]
    month_mean$SOS <- as.double(month_mean$SOS)
    month_mean[, SOS := scale(SOS, scale = FALSE), by = "siteID"]

    return(month_mean)
}

LmExtValid <- function(data_dt, test_sites) {
    LmSiteExtValid <- function(month_mean, test_sites) {
        # Split to training and test sets
        train_dt <- month_mean[siteID %in% test_sites, ]
        test_dt <- month_mean[!siteID %in% test_sites, ]

        lm_fit <- lm(SOS ~ Tmean, data = train_dt)
        pred <- predict(lm_fit, newdata = test_dt)

        return(data.table(site = test_dt$siteID, 
            obs = test_dt$SOS, pred = pred))
    }
    
    # Try single month means
    mod_fits <- list()
    for (i in 1:5) {
        month_mean <- unique(data_dt[month(Date) == i,
            .(Tmean = mean(Tmean), SOS = SOS),
            by = .(siteID, PhenoYear)
        ])

        month_mean <- Convert2Anomaly(month_mean)

        ext_valid_res <- LmSiteExtValid(month_mean, test_sites)
        rmse <- sqrt(mean((ext_valid_res$pred - ext_valid_res$obs)^2))

        mod_fits <- append(mod_fits, list(
            list(mod = ext_valid_res, rmse = rmse, type = paste0("month_", i))
        ))
    }

    # Try Jan-March mean
    month_mean <- unique(data_dt[month(Date) %in% c(1, 2, 3),
        .(Tmean = mean(Tmean), SOS = SOS),
        by = .(siteID, PhenoYear)
    ])
    month_mean <- Convert2Anomaly(month_mean)

    ext_valid_res <- LmSiteExtValid(month_mean, test_sites)
    rmse <- sqrt(mean((ext_valid_res$pred - ext_valid_res$obs)^2))

    mod_fits <- append(mod_fits, list(
        list(mod = ext_valid_res, rmse = rmse, type = paste0("Jan-Mar"))
    ))

    # Try March-May mean
    month_mean <- unique(data_dt[month(Date) %in% c(3, 4, 5),
        .(Tmean = mean(Tmean), SOS = SOS),
        by = .(siteID, PhenoYear)
    ])
    month_mean <- Convert2Anomaly(month_mean)

    ext_valid_res <- LmSiteExtValid(month_mean, test_sites)
    rmse <- sqrt(mean((ext_valid_res$pred - ext_valid_res$obs)^2))

    mod_fits <- append(mod_fits, list(
        list(mod = ext_valid_res, rmse = rmse, type = paste0("Mar-May"))
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




# Read data in
flux_dt <- fread(file.path(gdir, "Data/ARD/flux_mcd12q2_daymet_dt.csv"))
flux_dt[, Date := as_date(Date)]
# merge site meta
site_meta <- fread(file.path(gdir, "Data/ARD/sites_fluxnet2015_Tier1.csv"))
flux_dt <- merge(flux_dt, site_meta[, .(siteID, IGBP)], by = "siteID")
# filter out SOS > 255
flux_dt <- flux_dt[SOS < 255 & IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), ]

# Pooled linear regression
mod_mcd_dm_flux_lm_ext_valid <- LmExtValid(flux_dt, flux_split_li$test_sites)
# out: mod_mcd_dm_flux_lm_pooled_cv.Rds
saveRDS(mod_mcd_dm_flux_lm_pooled_cv, file.path(
    gdir, "Pipeline",
    "mod_mcd_dm_flux_lm_pooled_cv.Rds"
))

