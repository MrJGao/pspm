#'******************************************************************************
#' Description: Inspired by Hanninen et al 2019, I wonder if I simulate some 
#' spring phenology data by a "true" model, would other models be able to fit
#' the data as well?
#'******************************************************************************

# bsub -n 24 -W 36:00 -J "SimPheno" -oo Pipeline/hpc_job_log/para_var_out -eo Pipeline/hpc_job_log/para_var_err -R "rusage[mem=12GB]" "Rscript Code/mod_hpc_parallel_par_variation.R"

rm(list=ls())

library(data.table)
library(lubridate)
library(parallel)

source("src/base.R")
source("src/vis/vis_base.R")
source("src/hlp/hlp_fit_spm.R")
source("src/hlp/hlp_spring_pheno_models.R")
source("src/vis/vis_spring_pheno_model_diag_figs.R")
source("src/mod_compare_base.R")




# ~ Read in some temperature data observed from the HB ground ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hb_dt <- fread(file.path(hpc_dir, "Data/ARD", "hb_grd_temperature_pheno.csv"))
# unique(hb_dt[, .N / 365, by = .(siteID)])

# This site was selected b/c it has 27 site years
site_dt <- hb_dt[siteID == "1B"]

site_li <- FormatDataForPhenoModel(site_dt)


# ~ Simulation ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assume that the Paralle model is the true model that works in nature.
#   The PA model is selected b/c it has 9 parameters, more than other models.

# Parameter values
# Chilling-related
t0_chill <- 51
T_opt <- 2
T_min <- -5
T_max <- 10
C_min <- 0.2
C_req <- 100
# Forcing-related
t0 <- 140
T_base <- 5
F_crit <- 500

par_true <- c(t0, t0_chill, T_base, T_opt, T_min, T_max, C_min, F_crit, C_req)
# Simulate SOS data
sim_sos <- ParallelModel(par = par_true, data = site_li)


# Now, fit the similated sos data by the models
sim_li <- site_li
sim_li$transition_dates <- sim_sos

#' There's a concern that if the sos used to fit the model are exactly the same 
#' values as the simulated sos, it could cause numerical issues that lead to 
#' the unstable parameters. So, I tried to add some random noise to the 
#' simulated sos. But, I found the parameters are still unstable, which means
#' it's not caused by the numerical issues
#' 
# sim_li$transition_dates <- sim_sos + rnorm(length(sim_sos), 0, 1)


# # This step takes several minutes
# pheno_models <- FitCompareModels(sim_li)

# # compare fit
# par(mfrow = c(1, 5), mgp = c(1.5, 0.5, 0))
# par(mar = c(0, 0, 2, 0), oma = c(3, 3, 0, 1))
# ScatterXY(pheno_models$TT$fitted_doy, sim_sos,
#     range = c(150, 250),
#     xlab = "", ylab = "", main = "TT"
# )
# ScatterXY(pheno_models$PA$fitted_doy, sim_sos,
#     range = c(150, 250),
#     xlab = "", ylab = "", main = "PA", yaxt = "n"
# )
# ScatterXY(pheno_models$SQ$fitted_doy, sim_sos,
#     range = c(150, 250),
#     xlab = "", ylab = "", main = "SQ", yaxt = "n"
# )
# ScatterXY(pheno_models$AT$fitted_doy, sim_sos,
#     range = c(150, 250),
#     xlab = "", ylab = "", main = "AT", yaxt = "n"
# )
# ScatterXY(pheno_models$UN$fitted_doy, sim_sos,
#     range = c(150, 250),
#     xlab = "", ylab = "", main = "UN", yaxt = "n"
# )

# # estimated parameters
# par_PA <- pheno_models$PA$optim$par
# par_PA[1] <- site_li$doy[par_PA[1]]
# par_PA[2] <- site_li$doy[par_PA[2]]
# names(par_PA) <- c(
#     "t0", "t0chl", "T_base", "T_opt", "T_min", "T_max", "C_min",
#     "F_crit", "C_req"
# )

# par_true
# par_true[1] <- site_li$doy[par_true[1]] # convert to doy
# par_true[2] <- site_li$doy[par_true[2]] # convert to doy
# names(par_true) <- c(
#     "t0", "t0chl", "T_base", "T_opt", "T_min", "T_max", "C_min",
#     "F_crit", "C_req"
# )


# print(rbind(par_true, par_PA))


# Here, I fit the PA model multiple times but slightly change temperature to add 
# some randomness, I want a plot that shows the model parameters are very 
# unstable, and more importantly, they might not reveal the true parameters


cl <- makeCluster(detectCores() - 1, type = "SOCK")
calls <- clusterCall(cl, function() {
    suppressWarnings({
        library(data.table)
        library(lubridate)

        source("src/base.R")
        source("src/vis/vis_base.R")
        source("src/hlp/hlp_fit_spm.R")
        source("src/hlp/hlp_spring_pheno_models.R")
        source("src/vis/vis_spring_pheno_model_diag_figs.R")
        source("src/mod_compare_base.R")
    })
})
clusterExport(cl, c("site_li", "par_true"))


output <- clusterApply(cl, 1:1000, function(i) {
    # Simulate data
    sim_li <- site_li
    # Slightly change temperatures to add some randomness
    sim_li$Ti <- site_li$Ti + rnorm(length(site_li$Ti), 0, 0.01)
    # Simulate SOS
    sim_sos <- ParallelModel(par = par_true, data = sim_li)
    sim_li$transition_dates <- sim_sos

    # Fit the model using the simulated data
    mod_PA <- FitTheModel(
        model = ParallelModel,
        cost_fun = CostRMSE,
        data = sim_li
    )
    pars <- mod_PA$optim$par
    fit <- mod_PA$fitted_doy
    rmse <- mod_PA$optim$value

    return(list(est_par = c(rmse, pars), fit = fit))
})

saveRDS(output, 
    file.path(hpc_dir, 
        "Pipeline", 
        "parallel_prediction_variation.Rds"
    )
)



# set.seed(8964)
# est_par <- lapply(1:20, function(i) {
#     # Simulate data
#     sim_li <- site_li
#     # Slightly change temperatures to add some randomness
#     sim_li$Ti <- site_li$Ti + rnorm(length(site_li$Ti), 0, 0.01)
#     # Simulate SOS
#     sim_sos <- ParallelModel(par = par_true, data = sim_li)
#     sim_li$transition_dates <- sim_sos

#     # Fit the model using the simulated data
#     mod_PA <- FitTheModel(
#         model = ParallelModel,
#         cost_fun = CostRMSE,
#         data = sim_li,
#         control = list(maxit = 1e5)
#     )
#     pars <- mod_PA$optim$par
#     rmse <- mod_PA$optim$value
#     return(c(rmse, pars))
# })
# est_par <- do.call(rbind, est_par)
# est_par <- rbind(c(0, par_true), est_par)
# colnames(est_par) <- c(
#     "RMSE",
#     "t0", "t0chl", "T_base", "T_opt", "T_min", "T_max", "C_min",
#     "F_crit", "C_req"
# )

# pa_par_var <- list(par_true = par_true, est_par = est_par)

# # out: parallel_par_variation.Rds
# saveRDS(pa_par_var, 
#     file.path(hpc_dir, 
#         "Pipeline", 
#         "parallel_par_variation.Rds"
#     )
# )


