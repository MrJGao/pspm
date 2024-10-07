# Simulate some phenology data at HF using the TT model. Then, fix t0 and Tbase
#  and refit the model. Grab the residuals and see if they are correlated w/
#  temperature accumulated from t0 to SOS.

source("src/base.R")
SourceHelper()
library(data.table)
library(magrittr)

dt <- fread(file.path(gdir, "Data/ARD", "hb_hf_pheno_temp.csv"))


site_li <- FormatDataForPhenoModel(dt[siteID == "7T"])

# Simulate data
t0_true <- 80
T_base_true <- 0
F_crit_true <- 500

sim_pheno <- do.call(
    ThermalTimeModel, 
    list(par = c(t0_true, T_base_true, F_crit_true), data = site_li)
)
sim_li <- copy(site_li)
sim_li$transition_dates <- sim_pheno


# Now, fix t0 and Tbase
TTfix <- function(par, data, diagnose = FALSE) {
    # extract the parameter values from the
    # par argument for readability
    t0 <- 110
    T_base <- 5
    F_crit <- par

    # simple degree day sum setup
    Rf <- data$Ti - T_base
    Rf[Rf < 0] <- 0
    Rf[1:t0, ] <- 0

    # DOY of budburst criterium
    doy <- apply(Rf, 2, function(xt) {
        data$doy[which(cumsum(xt) >= F_crit)[1]]
    })

    doy <- StandardizeOutput(doy)

    if (diagnose == TRUE) {
        return(list(doy = doy, Rf = Rf))
    } else {
        return(doy)
    }
}

# Fit the model
est_par <- GenSA::GenSA(
    par = NULL,
    model = TTfix,
    fn = CostRMSE,
    data = sim_li,
    lower = c(0),
    upper = c(2000),
    control = NULL
)
est <- do.call(TTfix, list(par = est_par$par, data = sim_li))

err <- sim_li$transition_dates - est


# Now calculate accumulated temperature from t0 to SOS
acc_forc <- lapply(1:length(sim_li$transition_dates), function(i) {
    sos <- sim_li$transition_dates[i]
    temp <- sim_li$Ti[110:(sos + 110), i]
    temp[temp < 0] <- 0
    gdd <- sum(temp)
    return(gdd)
})
acc_forc <- do.call(c, acc_forc)


par(mfrow = c(1, 2))
CompareXY_J(est, sim_li$transition_dates, xlab = "Est.", ylab = "Sim.")
plot(acc_forc, err, pch = 16)
mod_lm <- lm(err ~ acc_forc)
LmLines(mod_lm)







