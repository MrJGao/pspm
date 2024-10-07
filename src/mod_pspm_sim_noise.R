# ******************************************************************************
# Simulations to explore how observational noise affects model intercomparison.
# ******************************************************************************
source("src/mod_hpc_pspm_sim_noise_base.R")
source("src/mod_compare_base.R")

# library(parallel)
library(magrittr)


# Reads in the arguments
args <- commandArgs(trailingOnly = TRUE)

# Throw error if missing or too many arguments
argslen <- length(args)
if (argslen != 1) stop("Error: number of arguments is incorrect")

# ~ Add random noise to the simulated data ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_idx <- as.numeric(args[1])

set.seed(run_idx)

# Retrieve simulation parameters
this_sim_setting <- sim_setting[run_idx, ]
model <- this_sim_setting[, "model"] %>% as.character()
noise <- this_sim_setting[, "noise"] %>% as.integer()
iter <- this_sim_setting[, "iter"] %>% as.integer()


sim_li <- copy(site_li)

sim_sos <- switch(model,
    "TT" = a_sos + rnorm(length(sim_li$transition_dates), 0, noise),
    "PA" = b_sos + rnorm(length(sim_li$transition_dates), 0, noise),
    "SQ" = c_sos + rnorm(length(sim_li$transition_dates), 0, noise),
    "AT" = d_sos + rnorm(length(sim_li$transition_dates), 0, noise)
)

sim_li$transition_dates <- sim_sos

mod_fit <- FitCompareModels(sim_li)

mod_fit$sim <- list(
    model = model,
    noise = noise,
    iter = iter
)

outdir <- file.path(hpc_dir, "Pipeline", "pspm_sim_noise")
if (dir.exists(outdir) == FALSE) {
    dir.create(outdir)
}

saveRDS(mod_fit, file = file.path(outdir, 
    paste0(model, "_", noise, "_", iter, ".Rds"))
)
