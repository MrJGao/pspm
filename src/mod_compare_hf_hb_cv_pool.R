#'******************************************************************************
#' Cross validation for HF & HB sites.
#'******************************************************************************


rm(list = ls())

library(parallel)

source("src/mod_compare_base.R")

args = commandArgs(trailingOnly = TRUE)
idx <- as.numeric(args[1])

hf_hb_li <- readRDS("pipe/hb_hf_li.Rds")


# ~ Pooled ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For pooled analysis, do leave-one-site-out cross validation
print("Fit HF & HB pooled cv....")


# Make cluster and fit pooled cv models
# cl <- makeCluster(detectCores(), type = "SOCK")
# # cl <- makeCluster(mpi.universe.size() - 1, type = "MPI")
# calls <- clusterCall(cl, function() {
#     suppressWarnings({
#         source("src/mod_compare_base.R")
#     })
# })
# clusterExport(cl, c("hf_hb_li"))


sites <- unique(hf_hb_li$site)
x <- sites[idx]

split_li <- LeaveOneSiteOut(data_li = hf_hb_li, site = x)
cv_res <- CvCompareModels(split_li)

saveRDS(cv_res, 
    file.path("pipe/mod_cv", "hf_hb_site", paste0(x, ".Rds"))
)


