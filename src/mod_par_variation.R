#'******************************************************************************
#' Description: Show that parameters of the process-based models could change 
#' dramatically if the data has a slight change.
#'******************************************************************************
source("src/base.R")
source("src/mod_compare_base.R")

library(data.table)



# in: read data in
hf_hb_li <- readRDS("pipe/hb_hf_li.Rds")

# ~ Process-based models ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fit spring phenology models
pheno_models <- FitCompareModels(hf_hb_li)

# Interpret the parameters
par_TT <- pheno_models$TT$optim$par
par_PA <- pheno_models$PA$optim$par
par_SQ <- pheno_models$SQ$optim$par
par_AT <- pheno_models$AT$optim$par
par_UN <- pheno_models$UN$optim$par

# Convert to DOY
par_TT[1] <- hf_hb_li$doy[par_TT[1]]
names(par_TT) <- c("t0", "T_base", "F_crit")

par_PA[1] <- hf_hb_li$doy[par_PA[1]]
par_PA[2] <- hf_hb_li$doy[par_PA[2]]
names(par_PA) <- c("t0", "t0chl", "T_base", "T_opt", "T_min", "T_max", "C_min", 
    "F_crit", "C_req")

par_SQ[1] <- hf_hb_li$doy[par_SQ[1]]
par_SQ[2] <- hf_hb_li$doy[par_SQ[2]]
names(par_SQ) <- c("t0", "t0chl", "T_base", "T_opt", "T_min", "T_max", 
    "F_crit", "C_req")

par_AT[1] <- hf_hb_li$doy[par_AT[1]]
names(par_AT) <- c("t0", "T_base", "a", "b", "c")


par_UN[1] <- hf_hb_li$doy[par_UN[1]]
names(par_UN) <- c("tc", "a_c", "b_c", "c_c", "b_f", "c_f", "w", "k", "C_req")


par_li <- list(par_TT, par_PA, par_SQ, par_AT, par_UN)


# ~ If we remove one site-year ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sub_li <- GetDataLiIndex(hf_hb_li, 
    (1:length(hf_hb_li$site))[-sample(1:length(hf_hb_li$site), 1)])

# Fit spring phenology models
pheno_models_sub <- FitCompareModels(sub_li)

# Interpret the parameters
par_sub_TT <- pheno_models_sub$TT$optim$par
par_sub_PA <- pheno_models_sub$PA$optim$par
par_sub_SQ <- pheno_models_sub$SQ$optim$par
par_sub_AT <- pheno_models_sub$AT$optim$par
par_sub_UN <- pheno_models_sub$UN$optim$par

# Convert to DOY
par_sub_TT[1] <- hf_hb_li$doy[par_sub_TT[1]]
names(par_sub_TT) <- c("t0", "T_base", "F_crit")

par_sub_PA[1] <- hf_hb_li$doy[par_sub_PA[1]]
par_sub_PA[2] <- hf_hb_li$doy[par_sub_PA[2]]
names(par_sub_PA) <- c(
    "t0", "t0chl", "T_base", "T_opt", "T_min", "T_max",
    "C_min", "F_crit", "C_req"
)

par_sub_SQ[1] <- hf_hb_li$doy[par_sub_SQ[1]]
par_sub_SQ[2] <- hf_hb_li$doy[par_sub_SQ[2]]
names(par_sub_SQ) <- c(
    "t0", "t0chl", "T_base", "T_opt", "T_min", "T_max",
    "F_crit", "C_req"
)

par_sub_AT[1] <- hf_hb_li$doy[par_sub_AT[1]]
names(par_sub_AT) <- c("t0", "T_base", "a", "b", "c")


par_sub_UN[1] <- hf_hb_li$doy[par_sub_UN[1]]
names(par_sub_UN) <- c(
    "tc", "a_c", "b_c", "c_c", "b_f", "c_f", "w", "k", "C_req"
)


par_sub_li <- list(par_sub_TT, par_sub_PA, par_sub_SQ, par_sub_AT, par_sub_UN)





# ~ For a single site-year ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Would it also change dramatically if just for a single site year
samp_site <- sample(unique(hf_hb_li$site), 1)
samp_site <- "1B"
samp_idx <- (1:length(hf_hb_li$site))[hf_hb_li$site == samp_site]

samp_site_li <- GetDataLiIndex(hf_hb_li, samp_idx)


site_model <- FitCompareModels(samp_site_li)

# Interpret the parameters
par_TT <- site_model$TT$optim$par
par_PA <- site_model$PA$optim$par
par_SQ <- site_model$SQ$optim$par
par_AT <- site_model$AT$optim$par
par_UN <- site_model$UN$optim$par

# Convert to DOY
par_TT[1] <- hf_hb_li$doy[par_TT[1]]
names(par_TT) <- c("t0", "T_base", "F_crit")

par_PA[1] <- hf_hb_li$doy[par_PA[1]]
par_PA[2] <- hf_hb_li$doy[par_PA[2]]
names(par_PA) <- c(
    "t0", "t0chl", "T_base", "T_opt", "T_min", "T_max", "C_min",
    "F_crit", "C_req"
)

par_SQ[1] <- hf_hb_li$doy[par_SQ[1]]
par_SQ[2] <- hf_hb_li$doy[par_SQ[2]]
names(par_SQ) <- c(
    "t0", "t0chl", "T_base", "T_opt", "T_min", "T_max",
    "F_crit", "C_req"
)

par_AT[1] <- hf_hb_li$doy[par_AT[1]]
names(par_AT) <- c("t0", "T_base", "a", "b", "c")

par_UN[1] <- hf_hb_li$doy[par_UN[1]]
names(par_UN) <- c("tc", "a_c", "b_c", "c_c", "b_f", "c_f", "w", "k", "C_req")


par_site_li <- list(par_TT, par_PA, par_SQ, par_AT, par_UN)


par(mfrow = c(2, 5), mar = c(2, 2, 1, 2))

plot(site_model$obs, site_model$TT$fitted_doy, 
    xlim = c(135, 170), ylim = c(135, 170)
)
plot(site_model$obs, site_model$PA$fitted_doy, 
    xlim = c(135, 170), ylim = c(135, 170)
)
plot(site_model$obs, site_model$SQ$fitted_doy, 
    xlim = c(135, 170), ylim = c(135, 170)
)
plot(site_model$obs, site_model$AT$fitted_doy, 
    xlim = c(135, 170), ylim = c(135, 170)
)
plot(site_model$obs, site_model$UN$fitted_doy, 
    xlim = c(135, 170), ylim = c(135, 170)
)



# ~ If we remove one site-year ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sub_li <- GetDataLiIndex(
    samp_site_li,
    (1:length(samp_site_li$site))[-sample(1:length(samp_site_li$site), 1)]
)

# Fit spring phenology models
site_models_sub <- FitCompareModels(sub_li)

# Interpret the parameters
par_sub_TT <- site_models_sub$TT$optim$par
par_sub_PA <- site_models_sub$PA$optim$par
par_sub_SQ <- site_models_sub$SQ$optim$par
par_sub_AT <- site_models_sub$AT$optim$par
par_sub_UN <- site_models_sub$UN$optim$par

# Convert to DOY
par_sub_TT[1] <- hf_hb_li$doy[par_sub_TT[1]]
names(par_sub_TT) <- c("t0", "T_base", "F_crit")

par_sub_PA[1] <- hf_hb_li$doy[par_sub_PA[1]]
par_sub_PA[2] <- hf_hb_li$doy[par_sub_PA[2]]
names(par_sub_PA) <- c(
    "t0", "t0chl", "T_base", "T_opt", "T_min", "T_max", "C_min",
    "F_crit", "C_req"
)

par_sub_SQ[1] <- hf_hb_li$doy[par_sub_SQ[1]]
par_sub_SQ[2] <- hf_hb_li$doy[par_sub_SQ[2]]
names(par_sub_SQ) <- c(
    "t0", "t0chl", "T_base", "T_opt", "T_min", "T_max",
    "F_crit", "C_req"
)

par_sub_AT[1] <- hf_hb_li$doy[par_sub_AT[1]]
names(par_sub_AT) <- c("t0", "T_base", "a", "b", "c")


par_sub_UN[1] <- hf_hb_li$doy[par_sub_UN[1]]
names(par_sub_UN) <- c(
    "tc", "a_c", "b_c", "c_c", "b_f", "c_f", "w", "k", "C_req"
)


par_site_sub_li <- list(par_sub_TT, par_sub_PA, par_sub_SQ, 
    par_sub_AT, par_sub_UN
)


plot(site_models_sub$obs, site_models_sub$TT$fitted_doy, 
    xlim = c(135, 170), ylim = c(135, 170)
)
plot(site_models_sub$obs, site_models_sub$PA$fitted_doy, 
    xlim = c(135, 170), ylim = c(135, 170)
)
plot(site_models_sub$obs, site_models_sub$SQ$fitted_doy, 
    xlim = c(135, 170), ylim = c(135, 170)
)
plot(site_models_sub$obs, site_models_sub$AT$fitted_doy, 
    xlim = c(135, 170), ylim = c(135, 170)
)
plot(site_models_sub$obs, site_models_sub$UN$fitted_doy, 
    xlim = c(135, 170), ylim = c(135, 170)
)
