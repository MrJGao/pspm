#'******************************************************************************
#' Cross-validation for process-based spring phenology models using 
#' MCD12Q2 + Daymet at flux sites.
#'******************************************************************************

args = commandArgs(trailingOnly = TRUE)
idx <- as.numeric(args[1])


source("src/mod_compare_base.R")


# ~ Site-specific ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For site-specific models, do leave-one-year-out cv
print("Fit MCD12Q2 + Daymet site-specific....")

# Make cluster and fit pooled cv models
# cl <- makeCluster(detectCores())
# calls <- clusterCall(cl, function() {
#     suppressWarnings({
#         source("src/mod_compare_base.R")
#     })
# })
# clusterExport(cl, c("flux_li"))

DoSite <- function(x) {
    flux_li <- readRDS("pipe/flux_li.Rds")
    # x <- "CA-Man"
    site_idx <- grep(x, flux_li$site)
    # If less than 10 site years, skip this site
    if (length(site_idx) < 10) return(NULL)
    
    site_li <- GetDataLiIndex(flux_li, site_idx)
    cv_res <- data.table()
    for (i in 1:length(site_idx)) {
        site_split_li <- SplitTrainTest(data_li = site_li, idx = i, 
            folds = length(site_idx))
        mod <- CvCompareModels(site_split_li)
        cv_res <- rbind(cv_res, data.table(
            site = x, 
            obs = mod$obs,
            pred_TT = mod$TT,
            pred_PA = mod$PA,
            pred_SQ = mod$SQ,
            pred_AT = mod$AT,
            pred_UN = mod$UN
        ))
    }

    # Make unconverged results to NA
    cv_res[cv_res < 0 | cv_res > 255] <- NA

    # Convert predictions to anomalies
    cv_res[, ":="(
        obs = scale(obs, scale = FALSE),
        pred_TT = scale(pred_TT, scale = FALSE),
        pred_PA = scale(pred_PA, scale = FALSE),
        pred_SQ = scale(pred_SQ, scale = FALSE),
        pred_AT = scale(pred_AT, scale = FALSE),
        pred_UN = scale(pred_UN, scale = FALSE)
    )]

    return(cv_res)
}

flux_li <- readRDS("pipe/flux_li.Rds")
sites <- unique(flux_li$site)

res <- DoSite(sites[idx])

saveRDS(res, 
    file.path("pipe/mod_cv/mcd_dm_flux_site", paste0(sites[idx], ".Rds"))
)


# mod_mcd_dm_flux_site_cv <- lapply(
#     unique(flux_li$site), 
    
# )

# names(mod_mcd_dm_flux_site_cv) <- unique(flux_li$site)

# # out: cv for MCD12Q2 + Daymet at flux sites site-specific
# saveRDS(mod_mcd_dm_flux_site_cv, 
#     file.path("pipe/mod_cv", "mod_mcd_dm_flux_site_cv.Rds")
# )


# shut down cluster
# stopCluster(cl)
