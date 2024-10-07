library(data.table)
library(lubridate)

source("src/base.R")
source("src/vis/vis_base.R")
source("src/hlp/hlp_fit_spm.R")
source("src/hlp/hlp_spring_pheno_models.R")
source("src/vis/vis_spring_pheno_model_diag_figs.R")
source("src/mod_compare_base.R")



hb_dt <- fread(file.path(hpc_local_dir, "Data/ARD", "hb_grd_temperature_pheno.csv"))
# unique(hb_dt[, .N / 365, by = .(siteID)])

# This site was selected b/c it has 27 site years
site_dt <- hb_dt[siteID == "1B"]

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

par_true <- c(t0, t0_chill, T_base, T_opt, T_min, T_max, C_min, F_crit, C_req)
# Simulate SOS data
sim_sos <- ParallelModel(par = par_true, data = site_li)


# Now, fit the similated sos data by the models
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
    data = sim_li,
    control = list(maxit = 1e5)
)
pars <- mod_PA$optim$par
rmse <- mod_PA$optim$value

est_par <- rbind(par_true, pars)
est_par



# The range of temperature noise
range(rnorm(length(site_li$Ti), 0, 0.01))
