#'******************************************************************************
#' Cross-validation for process-based spring phenology models using 
#' MCD12Q2 + Flux temperature at flux sites.
#'******************************************************************************


source("src/mod_compare_base.R")

args = commandArgs(trailingOnly = TRUE)
idx <- as.numeric(args[1])

outdir <- "pipe/mod_cv"
if (!dir.exists(outdir)) {
    dir.create(outdir)
}

# ~ Pooled ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For pooled analysis, do leave-one-site-out cross validation
print("Fit MCD12Q2 + flux_temp pooled cv....")

# Make cluster and fit pooled cv models
# cl <- makeCluster(detectCores(), type = "SOCK")
# calls <- clusterCall(cl, function() {
#     suppressWarnings({
#         source("src/mod_compare_base.R")
#     })
# })
# clusterExport(cl, c("flux_temp_li"))

flux_temp_li <- readRDS("pipe/flux_temp_li.Rds")
sites <- unique(flux_temp_li$site)

x <- sites[idx]

split_li <- LeaveOneSiteOut(data_li = flux_temp_li, site = x)
cv_res <- CvCompareModels(split_li)

saveRDS(cv_res,
    file.path(outdir, "mcd_flux_temp_pool", paste0(x, ".Rds"))
)

# mod_mcd_flux_temp_pooled_cv <- lapply(sites, function(x) {
#     split_li <- LeaveOneSiteOut(data_li = flux_temp_li, site = x)
#     cv_res <- CvCompareModels(split_li)

#     return(cv_res)
# })
# names(mod_mcd_flux_temp_pooled_cv) <- paste0("test_", sites)

# # out: cv for MCD12Q2 + flux_temp at flux sites pooled
# saveRDS(mod_mcd_flux_temp_pooled_cv, 
#     file.path(outdir, "mod_mcd_flux_temp_pooled_cv.Rds"))


# shut down cluster
# stopCluster(cl)