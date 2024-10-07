#'******************************************************************************
#' Cross validation for HF & HB sites.
#'******************************************************************************


rm(list = ls())

library(parallel)

source("src/mod_compare_base.R")

args = commandArgs(trailingOnly = TRUE)
idx <- as.numeric(args[1])



# ~ Site-specific ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For site-specific models, do leave-one-year-out cv
# Make cluster and fit pooled cv models
# cl <- makeCluster(detectCores(), type = "SOCK")
# calls <- clusterCall(cl, function() {
#     suppressWarnings({
#         source("src/mod_compare_base.R")
#     })
# })
# clusterExport(cl, c("hf_hb_li"))


hf_hb_li <- readRDS("pipe/hb_hf_li.Rds")

sites <- unique(hf_hb_li$site)
x <- sites[idx]

site_idx <- grep(x, hf_hb_li$site)
# If less than 10 site years, skip this site
if (length(site_idx) < 10) return(NULL)

site_li <- GetDataLiIndex(hf_hb_li, site_idx)
cv_res <- data.table()
for (i in 1:length(site_idx)) {
    site_split_li <- SplitTrainTest(
        data_li = site_li, 
        idx = i,
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


saveRDS(cv_res, file.path("pipe/mod_cv", "hf_hb_site", paste0(x, ".Rds")))



# mod_hf_hb_site_cv <- lapply(sites, function(x) {
#     site_idx <- grep(x, hf_hb_li$site)
#     # If less than 10 site years, skip this site
#     if (length(site_idx) < 10) return(NULL)

#     site_li <- GetDataLiIndex(hf_hb_li, site_idx)
#     cv_res <- data.table()
#     for (i in 1:length(site_idx)) {
#         site_split_li <- SplitTrainTest(data_li = site_li, idx = i,
#             folds = length(site_idx))
#         mod <- CvCompareModels(site_split_li)
#         cv_res <- rbind(cv_res, data.table(
#             site = x,
#             obs = mod$obs,
#             pred_TT = mod$TT,
#             pred_PA = mod$PA,
#             pred_SQ = mod$SQ,
#             pred_AT = mod$AT,
#             pred_UN = mod$UN
#         ))
#     }

#     # Make unconverged results to NA
#     cv_res[cv_res < 0 | cv_res > 255] <- NA

#     # Convert predictions to anomalies
#     cv_res[, ":="(
#         obs = scale(obs, scale = FALSE),
#         pred_TT = scale(pred_TT, scale = FALSE),
#         pred_PA = scale(pred_PA, scale = FALSE),
#         pred_SQ = scale(pred_SQ, scale = FALSE),
#         pred_AT = scale(pred_AT, scale = FALSE),
#         pred_UN = scale(pred_UN, scale = FALSE)
#     )]

#     return(cv_res)
# })

# names(mod_hf_hb_site_cv) <- unique(hf_hb_li$site)

# # out: cv for HFHB site-specific
# saveRDS(mod_hf_hb_site_cv, file.path("pipe/mod_cv", "mod_hf_hb_site_cv.Rds"))


# shut down cluster
# stopCluster(cl)
