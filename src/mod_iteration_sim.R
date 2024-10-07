#'******************************************************************************
#' Description: Use simulation to explore the relationship between parameter 
#' values and the maximum number of iteration of the simulated annealing. 
#' The idea is that I simulate some SOS data using the Parellel model, and fit 
#' the data using the model again but with different maximum number of iteration
#' times. For each maximum number of iteration setting, the RMSE and parameter 
#' values will be stored. 
#'******************************************************************************

# bsub < src/mod_hpc_iteration_sim.csh

source("src/base.R")
source("src/hlp/hlp_fit_spm.R")
source("src/hlp/hlp_spring_pheno_models.R")
source("src/mod_compare_base.R")

library(data.table)

library(Rmpi)
library(parallel)
library(snow)



# in: 
# Ground phenology and temperature at HB
hb_pheno_temp_file <- file.path(
    hpc_dir, "Data", 
    "hb_grd_temperature_pheno.csv"
)

# out: PA_sim_iter_same_data.csv
sim_iter_same_file <- file.path(
    hpc_dir, "Pipeline", 
    "PA_sim_iter_same_data.csv"
)
# out: PA_sim_iter_noisy_data.csv
sim_iter_noisy_file <- file.path(
    hpc_dir, "Pipeline", 
    "PA_sim_iter_noisy_data.csv"
)
# out: PA_sim_par_true.csv
sim_par_true_file <- file.path(hpc_dir, "Pipeline", "PA_sim_par_true.csv")


# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the site data
hb_dt <- fread(hb_pheno_temp_file)
# Use site `1B` as it has 27 site years of data
site_dt <- hb_dt[siteID == "1B"]

# Format the data to person-period format
site_li <- FormatDataForPhenoModel(site_dt)


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

# Simulate SOS by the true parameters
par_true <- c(t0, t0_chill, T_base, T_opt, T_min, T_max, C_min, F_crit, C_req)
sim_sos <- ParallelModel(par = par_true, data = site_li)

# Export true parameters
par_true_dt <- data.table(par_true)
par_true_dt$name <- c(
    "t0", "t0_chil", "T_base", "T_opt", "T_min", "T_max", "C_min", 
    "F_crit", "C_req"
)
fwrite(par_true_dt, file = sim_par_true_file)


# Format the simulated data
sim_li <- site_li
sim_li$transition_dates <- sim_sos

# Make cluster
# cl <- makeCluster(detectCores() - 1, type = "SOCK")
cl <- makeCluster(mpi.universe.size() - 1, type = "MPI")
calls <- clusterCall(cl, function() {
    suppressWarnings({
        source("src/base.R")
        source("src/hlp/hlp_fit_spm.R")
        source("src/hlp/hlp_spring_pheno_models.R")
        source("src/mod_compare_base.R")
        library(data.table)
    })
})
clusterExport(cl, c("site_li", "par_true", "sim_li"))


# All maxit values want to try
# iter_vec <- seq(1000, 20000, by = 500)
iter_vec <- c(
    seq(1, 1000, by = 50), 
    seq(1e3, 1e4, by = 1e2), 
    seq(1e4, 1e5, by = 2e3), 
    seq(1e5, 1e6, by = 5e4)
)


# Do the job
# Same data but different iteration times and seeds
output <- clusterApply(cl, iter_vec, function(iter_val) {
    # Fit the model using the simulated data
    mod_PA <- FitTheModel(
        model = ParallelModel,
        cost_fun = CostRMSE,
        data = sim_li,
        control = list(maxit = iter_val, seed = iter_val)
    )
    pars <- mod_PA$optim$par
    
    res_dt <- cbind(iter_val, mod_PA$optim$value, mod_PA$optim$counts, t(pars))

    return(res_dt)
})
output <- do.call(rbind, output)
colnames(output) <- c(
    "Maxit", "RMSE", "Num_call", 
    "t0", "t0_chil", "T_base", "T_opt", "T_min", "T_max", "C_min", 
    "F_crit", "C_req"
)
output <- data.table(output)


# Export
fwrite(output, file = sim_iter_same_file)



# Add some tiny noise into the data. seeds are different as well
output2 <- clusterApply(cl, iter_vec, function(iter_val) {
    set.seed(iter_val)
    # Add noise
    sim_li$Ti <- site_li$Ti + rnorm(length(site_li$Ti), 0, 0.01)

    mod_PA <- FitTheModel(
        model = ParallelModel,
        cost_fun = CostRMSE,
        data = sim_li,
        control = list(maxit = iter_val, seed = iter_val)
    )
    pars <- mod_PA$optim$par

    res_dt <- cbind(iter_val, mod_PA$optim$value, mod_PA$optim$counts, t(pars))

    return(res_dt)
})
output2 <- do.call(rbind, output2)
colnames(output2) <- c(
    "Maxit", "RMSE", "Num_call",
    "t0", "t0_chil", "T_base", "T_opt", "T_min", "T_max", "C_min", 
    "F_crit", "C_req"
)
output2 <- data.table(output2)

# Export
fwrite(output2, file = sim_iter_noisy_file)



# shut down cluster
snow::stopCluster(cl)
Rmpi::mpi.quit()


