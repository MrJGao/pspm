#'******************************************************************************
#' Description: Test goodness of fit for the simulated data.
#'******************************************************************************
source("src/mod_hpc_pspm_spec_sim_base.R")


library(parallel)

cl <- makeCluster(detectCores() - 1, type = "SOCK")
calls <- clusterCall(cl, function() {
    suppressWarnings({
        source("src/mod_hpc_pspm_spec_sim_base.R")
    })
})


output <- clusterApply(cl, seq_along(com_li), function(i) {
    i_li <- com_li[[i]]
    spec <- names(com_li)[i]
    fit <- FitCompareModels(i_li)
    return(list(spec = spec, fit = fit))
})

# out: pspm_spec_sim_fit.Rds
saveRDS(output, file.path(hpc_dir, "Pipeline", "pspm_spec_sim_fit.Rds"))


stopCluster(cl)





