#'******************************************************************************
#' Description: Bootstrapping to test model parameter correlations.
#'******************************************************************************

# bsub < src/mod_hpc_para_cor.csh

source("src/mod_compare_base.R")

library(Rmpi)
library(parallel)
library(snow)



FitModelPara <- function(data_li) {
    # TT
    mod_TT <- FitTheModel(
        model = ThermalTimeModel,
        cost_fun = CostRMSE,
        data = data_li,
        control = list(maxit = 1e5)
    )

    # PA
    mod_PA <- FitTheModel(
        model = ParallelModel,
        cost_fun = CostRMSE,
        data = data_li,
        control = list(maxit = 1e5)
    )

    # SQ
    mod_SQ <- FitTheModel(
        model = SequentialModel,
        cost_fun = CostRMSE,
        data = data_li,
        control = list(maxit = 1e5)
    )

    # AT
    mod_AT <- FitTheModel(
        model = AlternatingModel,
        cost_fun = CostRMSE,
        data = data_li,
        control = list(maxit = 1e5)
    )

    # UN
    mod_UN <- FitTheModel(
        model = UnifiedModel,
        cost_fun = CostRMSE,
        data = data_li,
        control = list(maxit = 1e5)
    )

    return(list(
        TT = mod_TT$optim$par, PA = mod_PA$optim$par, SQ = mod_SQ$optim$par, 
        AT = mod_AT$optim$par, UN = mod_UN$optim$par
    ))
}



SampleWithRep <- function(data_li) {
    samp_idx <- sample(
        1:length(data_li$site), 
        length(data_li$site), 
        replace = TRUE
    )

    samp_li <- list()
    samp_li$site <- data_li$site[samp_idx]
    samp_li$year <- data_li$year[samp_idx]
    samp_li$transition_dates <- data_li$transition_dates[samp_idx]
    samp_li$doy <- data_li$doy
    samp_li$Ti <- data_li$Ti[, samp_idx]
    samp_li$T_min <- data_li$T_min[, samp_idx]
    samp_li$T_max <- data_li$T_max[, samp_idx]

    return(samp_li)
}

data_li <- readRDS(file.path("pipe", "hb_hf_li.Rds"))

# Make cluster
# cl <- makeCluster(detectCores() - 1, type = "SOCK")
cl <- makeCluster(mpi.universe.size() - 1, type = "MPI")
calls <- clusterCall(cl, function() {
    suppressWarnings({
        source("src/mod_compare_base.R")
    })
})
clusterExport(cl, c("SampleWithRep", "FitModelPara", "data_li"))


# Sample 500 times
# B/c I changed the maxit to 1e5, it takes a very long time to run. So, I 
# submitted the job as individual 100 iteration sub-jobs and reduced the total
# number of iteration to 500.
# system.time({
    output <- clusterApply(cl, x = 1:600, function(x) {
        temp_file <- file.path(hpc_dir, 
            "Pipeline", "bs_temp",
            paste0("iter_", x, ".Rds")
        )
        if (file.exists(temp_file)) {
            return(NULL)
        }
        # Resample with replacement
        samp_li <- SampleWithRep(data_li)
        # # Fit the 5 models
        fit <- FitModelPara(samp_li)

        # Save result to a temporary file
        saveRDS(fit, file = temp_file)
        return(fit)
    })
# })


# shut down cluster
snow::stopCluster(cl)
Rmpi::mpi.quit()



# ~ Clean the result ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This part should run on hpc.
# bsub -Is -n 1 -W 20 tcsh

if (FALSE) {
    # Gather all results
    source("src/base.R")

    output_files <- list.files(file.path(hpc_dir, "Pipeline", "bs_temp"),
        "*.Rds$",
        full.names = TRUE
    )

    output <- lapply(output_files, function(x) {
        fit <- readRDS(x)
        return(fit)
    })

    # Extract and format result
    TT_par <- PA_par <- SQ_par <- AT_par <- UN_par <- NULL
    res <- lapply(output, function(ot) {
        TT_par <<- rbind(TT_par, ot$TT)
        PA_par <<- rbind(PA_par, ot$PA)
        SQ_par <<- rbind(SQ_par, ot$SQ)
        AT_par <<- rbind(AT_par, ot$AT)
        UN_par <<- rbind(UN_par, ot$UN)
    })

    # Assign parameter names
    TT_par <- data.table(TT_par)
    names(TT_par) <- c("t0", "T_base", "F_crit")

    PA_par <- data.table(PA_par)
    names(PA_par) <- c(
        "t0", "t0_chil", "T_base", "T_opt", "T_min", "T_max", "C_min",
        "F_crit", "C_req"
    )

    SQ_par <- data.table(SQ_par)
    names(SQ_par) <- c(
        "t0", "t0_chil", "T_base", "T_opt", "T_min", "T_max",
        "F_crit", "C_req"
    )

    AT_par <- data.table(AT_par)
    names(AT_par) <- c("t0", "T_base", "a", "b", "c")

    UN_par <- data.table(UN_par)
    names(UN_par) <- c("tc", "a_c", "b_c", "c_c", "b_f", "c_f", "w", "k", "C_req")


    # out: model parameter bootstrapping
    saveRDS(
        list(TT = TT_par, PA = PA_par, SQ = SQ_par, AT = AT_par, UN = UN_par),
        file.path(hpc_dir, "Pipeline", "mod_par_cor_1e5.Rds")
    )

}









