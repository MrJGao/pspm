#'******************************************************************************
#' Description: Simulate daily mean temperature values to test spring pheno 
#' models.
#'******************************************************************************

source("src/base.R")
source("src/mod_compare_base.R")

library(data.table)



# Use HF & HB data as an example.
# Fit different models using the same data, then simulate a years' daily 
# temperature to test the model results.


# ~ Fit models ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Process-based models ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# in: read data in
hf_hb_li <- readRDS(file.path(gdir, "Pipeline", "hb_hf_li.Rds"))
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
# in: read data in
hf_hb_li <- readRDS(file.path(gdir, "Pipeline", "hb_hf_li.Rds"))
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
    "t0", "t0chl", "T_base", "T_opt", "T_min", "T_max", "C_min",
    "F_crit", "C_req"
)

par_sub_SQ[1] <- hf_hb_li$doy[par_sub_SQ[1]]
par_sub_SQ[2] <- hf_hb_li$doy[par_sub_SQ[2]]
names(par_sub_SQ) <- c("t0", "t0chl", "T_base", "T_opt", "T_min", "T_max", 
    "F_crit", "C_req")

par_sub_AT[1] <- hf_hb_li$doy[par_sub_AT[1]]
names(par_sub_AT) <- c("t0", "T_base", "a", "b", "c")


par_sub_UN[1] <- hf_hb_li$doy[par_sub_UN[1]]
names(par_sub_UN) <- c(
    "tc", "a_c", "b_c", "c_c", "b_f", "c_f", "w", "k", "C_req"
)


par_sub_li <- list(par_sub_TT, par_sub_PA, par_sub_SQ, par_sub_AT, par_sub_UN)



# ~ Linear regression ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For this dataset, the monthly mean temperatures in May fit the best model
LmFitMay <- function(dt) {
    # Convert to site anomalies
    Convert2Anomaly <- function(month_mean) {
        month_mean[, Tmean := scale(Tmean, scale = FALSE), by = "siteID"]
        month_mean$SOS <- as.double(month_mean$SOS)
        month_mean[, SOS := scale(SOS, scale = FALSE), by = "siteID"]
        
        return(month_mean)
    }

    month_mean <- unique(dt[month(Date) == 5, .(Tmean = mean(Tmean), SOS = SOS),
        by = .(siteID, PhenoYear)
    ])
    
    month_mean <- Convert2Anomaly(month_mean)
    
    # plot(month_mean[siteID == "US-Ha1", .(SOS, Tmean)])
    # plot(month_mean[, .(SOS, Tmean)], 
    #     col = factor(month_mean$siteID), 
    #     pch = 16
    # )
    lm_fit <- lm(SOS ~ Tmean, data = month_mean)
    # summary(lm_fit)
    rmse <- sqrt(mean(lm_fit$residuals^2))

    return(lm_fit)
}

# HF
hf_pheno_dt <- fread(file.path(gdir, "Data/ARD/hf_grd_temperature_pheno.csv"))
hf_pheno_dt[, Date := as_date(Date)]

# HB
hb_pheno_dt <- fread(file.path(gdir, "Data/ARD/hb_grd_temperature_pheno.csv"))
hb_pheno_dt[, Date := as_date(Date)]

hb_hf_dt <- rbind(
    hf_pheno_dt[, .(siteID, Date,
        Tmean = Temp, Tmax = NA, Tmin = NA, SOS,
        PhenoYear
    )],
    hb_pheno_dt[, .(siteID, Date, Tmean, Tmax, Tmin, SOS, PhenoYear)]
)

# Fit linear regression
mod_lm_pooled <- LmFitMay(hb_hf_dt)
confint(mod_lm_pooled)




# ~ Simulation ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Randomly take one site-years' temperatures and add 2 degree C on each day
# uni_sy <- unique(hb_hf_dt[, .(siteID, PhenoYear)])
# idx <- sample(1:nrow(uni_sy), size = 1)
# samp_sy <- uni_sy[idx,]
# samp_dt <- hb_hf_dt[siteID == samp_sy$siteID, ]

# Calculate site mean temperatures
hb_hf_dt[, ":=" (temp_year = year(Date), temp_doy = yday(Date))]

sim_dt <- data.table()
for (site in unique(hb_hf_dt$siteID)) {
    temp_sim_dt <- unique(hb_hf_dt[siteID == site & temp_doy != 366, .(
        SOS = mean(SOS), Tmean = mean(Tmean), PhenoYear = 2100,
        siteID, Tmax = NA, Tmin = NA
    ), by = "temp_doy"])
    temp_sim_dt$Date <- as.Date("2100-01-01") + hf_hb_li$doy
    # simulate temperature
    temp_sim_dt$Tmean[1:110] <- temp_sim_dt$Tmean[1:110] + 
        rnorm(110, mean = 2, sd = 2)
    sim_dt <- rbind(sim_dt, temp_sim_dt)
}

samp_sim_li <- FormatDataForPhenoModel(sim_dt)


# ~ Process-based model predictions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pred_TT <- do.call(ThermalTimeModel, list(par = pheno_models$TT$optim$par, 
    data = samp_sim_li))
pred_PA <- do.call(ParallelModel, list(par = pheno_models$PA$optim$par, 
    data = samp_sim_li))
pred_SQ <- do.call(SequentialModel, list(par = pheno_models$SQ$optim$par, 
    data = samp_sim_li))
pred_AT <- do.call(AlternatingModel, list(par = pheno_models$AT$optim$par, 
    data = samp_sim_li))
pred_UN <- do.call(UnifiedModel, list(par = pheno_models$UN$optim$par, 
    data = samp_sim_li))

preds <- cbind(pred_TT, pred_PA, pred_SQ, pred_AT, pred_UN)
preds[preds > 300] <- NA
boxplot(preds, ylim = range(c(preds, sim_dt[, mean(SOS)]), na.rm = TRUE))
abline(h = sim_dt[, mean(SOS)])


# B/c linear regression performs on site anomalies, here we first need to get 
# the sampled site, calculate average value of the mean temperatures in May, 
# simulate temperature, and center it using the same average temperature value, 
# then we can retrieve the sample year out.

# Get site means
samp_dt_month_mean <- unique(samp_dt[month(Date) == 5, 
    .(Tmean = mean(Tmean), SOS = SOS),
    by = .(siteID, PhenoYear)])
# Centering values
temp_ctr <- mean(samp_dt_month_mean$Tmean)
sos_ctr <- mean(samp_dt_month_mean$SOS)
# Simulate temperature for the sampled year
sim_dt <- samp_dt_month_mean[PhenoYear == samp_sy$PhenoYear, 
    .(Tmean = Tmean + 2 - temp_ctr)]

pred_sos <- predict(mod_lm_pooled, 
    newdata = sim_dt, 
    interval = "predict"
) + sos_ctr
obs_sos <- samp_dt_month_mean[PhenoYear == samp_sy$PhenoYear, SOS]




