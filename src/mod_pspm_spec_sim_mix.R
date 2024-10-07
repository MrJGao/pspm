#'******************************************************************************
#' Description: Test model performance on aggregated data.
#'******************************************************************************
source("src/mod_hpc_pspm_spec_sim_base.R")

library(parallel)


# Reads in the arguments
args <- commandArgs(trailingOnly = TRUE)

# Throw error if missing or too many arguments
argslen <- length(args)
if (argslen != 1) stop("Error: number of arguments is incorrect")


# ~ Add random noise to the simulated data ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_idx <- as.numeric(args[1])

# Set seed for reproducing
set.seed(run_idx)

for (i in seq_along(mix_li)) {
    sos <- mix_li[[i]]$transition_dates
    sos_noise <- sos + rnorm(length(sos), 0, 1)
    mix_li[[i]]$transition_dates <- sos_noise
}



# Set the cluster
cl <- makeCluster(detectCores() - 1, type = "SOCK")
calls <- clusterCall(cl, function() {
    suppressWarnings({
        source("src/mod_hpc_pspm_spec_sim_base.R")
    })
})



# ~ Goodness-of-fit ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If the job has already been run, skip it
gfit_output <- file.path(
    hpc_dir, "Pipeline",
    "pspm_mix_cv",
    paste0("pspm_spec_mix_fit_", run_idx, ".Rds")
)

if (!file.exists(gfit_output)) {
    gfit <- clusterApply(cl, seq_along(mix_li), function(i) {
        i_li <- mix_li[[i]]
        spec <- names(mix_li)[i]
        fit <- FitCompareModels(i_li)
        return(list(spec = spec, fit = fit))
    })

    # out: pspm_spec_mix_fit.Rds
    saveRDS(gfit, gfit_output)
}




# ~ Cross-validation ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Number of cross-validation folds
num_folds <- length(m1_li$site) / 1

clusterExport(cl, c("num_folds"))


for (i in seq_along(mix_li)) {
    data_li <- mix_li[[i]]
    spec <- names(mix_li)[i]
    
    clusterExport(cl, c("data_li", "spec"))

    cv_output <- file.path(
        hpc_dir, "Pipeline", "pspm_mix_cv",
        paste0("pspm_spec_mix_cv_", run_idx, "_", i, ".Rds")
    )

    if (!file.exists(cv_output)) {
        cv_res <- clusterApply(cl, x = 1:num_folds, function(i) {
            split_li <- SplitTrainTest(data_li, i, folds = num_folds)
            fit <- CvCompareModels2(split_li)
            return(list(spec = spec, fit = fit))
        })

        # out: pspm_spec_sim_cv_*.Rds
        saveRDS(cv_res, cv_output)
    }
}


stopCluster(cl)



