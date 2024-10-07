#'******************************************************************************
#' Description: Test cross-validation for the simulated data.
#'******************************************************************************
source("src/mod_hpc_pspm_spec_sim_base.R")

library(parallel)



# Reads in the arguments
args <- commandArgs(trailingOnly = TRUE)

# Throw error if missing or too many arguments
argslen <- length(args)
if (argslen != 1) stop("Error: number of arguments is incorrect")

# The data
data_li <- com_li[[as.numeric(args[1])]]
# The species name
spec <- names(com_li)[as.numeric(args[1])]


cl <- makeCluster(detectCores() - 1, type = "SOCK")
calls <- clusterCall(cl, function() {
    suppressWarnings({
        source("src/mod_hpc_pspm_spec_sim_base.R")
    })
})

# Number of cross-validation folds
num_folds <- length(data_li$site) / 1

clusterExport(cl, c("num_folds", "data_li", "spec"))

cv_res <- clusterApply(cl, x = 1:num_folds, function(i) {
    split_li <- SplitTrainTest(data_li, i, folds = num_folds)
    fit <- CvCompareModels2(split_li)
    return(list(spec = spec, fit = fit))
})

# out: pspm_spec_sim_cv_*.Rds
saveRDS(cv_res, file.path(hpc_dir, "Pipeline", 
    paste0("pspm_spec_sim_cv_", spec, ".Rds")
))


stopCluster(cl)



