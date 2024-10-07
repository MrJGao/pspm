# ******************************************************************************
# Simulations to explore how observational noise affects model intercomparison.
# 
# This file contains simulated data and base functions.
# ******************************************************************************
source("src/base.R")
source("src/hlp/hlp_fit_spm.R")
library(data.table)


# Use HB data as an example
hfhb_dt <- fread(file.path(hpc_dir, "Data/ARD/hb_hf_pheno_temp.csv"))
site_li <- FormatDataForPhenoModel(hfhb_dt[siteID == "1B"])


# TT model true parameters
a_tt_true <- c(
    t0 = 140,
    T_base = 5,
    F_crit = 300
)
a_sos <- do.call(ThermalTimeModel, list(
    par = a_tt_true, data = site_li
))

# PA model true parameters
b_pa_true <- c(
    t0 = 110,
    t0_chill = 90,
    T_base = 2,
    T_opt = 7,
    T_min = 0,
    T_max = 10,
    C_min = 0.1,
    F_crit = 200,
    C_req = 150
)
b_sos <- do.call(ParallelModel, list(
    par = b_pa_true, data = site_li
))

# SQ model true parameters
c_sq_true <- c(
    t0 = 100,
    t0_chill = 60,
    T_base = 5,
    T_opt = 2,
    T_min = -15,
    T_max = 20,
    F_crit = 400,
    C_req = 5
)
c_sos <- do.call(SequentialModel, list(
    par = c_sq_true, data = site_li
))

# AT model true parameters
d_at_true <- c(
    t0 = 150,
    T_base = 3,
    a = 50,
    b = 5000,
    c = -0.05
)
d_sos <- do.call(AlternatingModel, list(
    par = d_at_true, data = site_li
))


# Simulate setting
noise_level <- seq(1, 10, by = 1)
models <- c("TT", "PA", "SQ", "AT")
iter_idx <- 1:100

sim_setting <- expand.grid(noise = noise_level, model = models, iter = iter_idx)

