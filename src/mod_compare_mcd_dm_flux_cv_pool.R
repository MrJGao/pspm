#'******************************************************************************
#' Cross-validation for process-based spring phenology models using 
#' MCD12Q2 + Daymet at flux sites. For pooled analysis
#'******************************************************************************


source("src/mod_compare_base.R")

args = commandArgs(trailingOnly = TRUE)
idx <- as.numeric(args[1])

flux_li <- readRDS("pipe/flux_li.Rds")


# ~ Pooled ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For pooled analysis, do leave-one-site-out cross validation
print("Fit MCD12Q2 + Daymet pooled cv....")

# Make cluster and fit pooled cv models
# cl <- makeCluster(46)
# calls <- clusterCall(cl, function() {
#     suppressWarnings({
#         source("src/mod_compare_base.R")
#     })
# })
# clusterExport(cl, c("flux_li"))

sites <- unique(flux_li$site)

x <- sites[idx]

split_li <- LeaveOneSiteOut(data_li = flux_li, site = x)
cv_res <- CvCompareModels(split_li)

saveRDS(cv_res, 
    file.path("pipe/mod_cv", "mcd_dm_flux_pool", 
        paste0(x, ".Rds")
    )
)

# mod_mcd_dm_flux_pooled_cv <- lapply(sites, function(x) {
#     split_li <- LeaveOneSiteOut(data_li = flux_li, site = x)
#     cv_res <- CvCompareModels(split_li)

#     return(cv_res)
# })
# names(mod_mcd_dm_flux_pooled_cv) <- paste0("test_", sites)

# # out: cv for MCD12Q2 + Daymet at flux sites pooled
# saveRDS(mod_mcd_dm_flux_pooled_cv, 
#     file.path("pipe/mod_cv", "mod_mcd_dm_flux_pooled_cv.Rds")
# )

# stopCluster(cl)

