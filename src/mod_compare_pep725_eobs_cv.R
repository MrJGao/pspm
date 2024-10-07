#'******************************************************************************
#' Description: Cross-validation for PEP725 + E-OBS and save 
#'              the result to a file that will be used later to make figures.
#'******************************************************************************

# bsub < Code/mod_compare_hpc_pep725_eobs_cv.csh

source("src/mod_compare_base.R")

args = commandArgs(trailingOnly = TRUE)
idx <- as.numeric(args[1])


# pep725_species <- c("aesculus", "alnus", "betula", "fagus", 
#     "fraxinus", "quercus"
# )
# species_name <- pep725_species[5] # for debug

# Set up cluster
# cl <- makeCluster(mpi.universe.size() - 1, type = "MPI")
# # cl <- makeCluster(24, type = "SOCK")
# calls <- clusterCall(cl, function() {
#     suppressWarnings({
#         source("Code/mod_compare_base.R")
#     })
# })


# Fit species-specific models
caches <- list.files("pipe/pep_tmp", full.names = TRUE)

outdir <- file.path("pipe/mod_cv", "pep")
if (!dir.exists(outdir)) {
    dir.create(outdir)
}

print(paste("CV models......", basename(caches[idx])))

cur_site <- readRDS(caches[idx])

# Leave one year out cv
cv_res <- data.table()
for (j in sample(unique(cur_site$PhenoYear), 10)) {
    system.time({
        cur_site_train <- cur_site[PhenoYear != j, ]
        cur_site_test <- cur_site[PhenoYear == j, ]

        cur_train_li <- FormatDataForPhenoModel(cur_site_train)
        cur_test_li <- FormatDataForPhenoModel(cur_site_test)
        split_li <- list(train = cur_train_li, test = cur_test_li)

        mod <- CvCompareModels(split_li, if_un_single = TRUE)

        cv_res <- rbind(cv_res, data.table(
            site = caches[idx],
            obs = mod$obs,
            pred_TT = mod$TT,
            pred_PA = mod$PA,
            pred_SQ = mod$SQ,
            pred_AT = mod$AT,
            pred_UN = mod$UN
        ))
    })
}
# Make unconverged results to NA
cv_res[cv_res < 0 | cv_res > 255] <- NA

saveRDS(cv_res, file.path(outdir, basename(caches[idx])))


# pep_eobs <- clusterApply(cl = cl, caches, function(x) {
#     cur_site <- readRDS(x)

#     # Leave one year out cv
#     cv_res <- data.table()
    
#     for (j in sample(unique(cur_site$PhenoYear), 10)) {
#         system.time({
#             cur_site_train <- cur_site[PhenoYear != j, ]
#             cur_site_test <- cur_site[PhenoYear == j, ]

#             cur_train_li <- FormatDataForPhenoModel(cur_site_train)
#             cur_test_li <- FormatDataForPhenoModel(cur_site_test)
#             split_li <- list(train = cur_train_li, test = cur_test_li)

#             mod <- CvCompareModels(split_li)

#             cv_res <- rbind(cv_res, data.table(
#                 site = x,
#                 obs = mod$obs,
#                 pred_TT = mod$TT,
#                 pred_PA = mod$PA,
#                 pred_SQ = mod$SQ,
#                 pred_AT = mod$AT,
#                 pred_UN = mod$UN
#             ))
#         })
#     }
#     # Make unconverged results to NA
#     cv_res[cv_res < 0 | cv_res > 255] <- NA
    
#     return(cv_res)
# })

# names(pep_eobs) <- lapply(caches, function(x) {
#     strsplit(basename(x), "\\.")[[1]][1]
# })

# out: model fit for PEP725 + E-OBS
# saveRDS(pep_eobs, file.path(
#     hpc_dir, "Pipeline",
#     paste0("mod_", species_name, "_pep725_ebos_cv", ".Rds")
# ))

# shut down cluster
# snow::stopCluster(cl)
# Rmpi::mpi.quit()
