#'******************************************************************************
#' Description: Survival analysis for observed spring phenology data including
#' - MCD12Q2 + Daymet
#' - MCD12Q2 + flux_temperature
#' - HF & HB + ground_temperature
#' - PEP725 + E-OBS
#'******************************************************************************

rm(list = ls())

source("src/base.R")
source("src/hlp/hlp_chilling_calc.R")
source("src/hlp/hlp_model_cost_functions.R")
source("src/hlp/hlp_spring_pheno_models.R")
source("src/hlp/hlp_fit_spm.R")

library(survival)
library(data.table)
library(magrittr)
library(lubridate)


# Calculate chilling and forcing from a person-period dataset
CalChillForc <- function(pp_dt,
    Tc = 0, tc = -20,
    Tb = 5, t0 = 1
) {
    # Calculate chilling
    # Let's use chill days for now
    pp_dt[, ":=" (chill = 0, acc_chill = 0)]
    pp_dt[doy > tc,
        chill := .SD[, ChillFix(Tmean, "h")],
        by = .(id)
    ]
    # Accumulated chilling
    pp_dt[doy > tc,
        acc_chill := .SD[, cumsum(chill)],
        by = .(id)
    ]


    # Calculate forcing
    pp_dt[, ":="(forc = 0, acc_forc = 0)]
    pp_dt[doy > t0,
        forc := .SD[
            ,
            ifelse(Tmean - Tb > 0, Tmean - Tb, 0)
        ],
        by = .(id)
    ]
    # Accumulated forcing
    pp_dt[doy > t0,
        acc_forc := .SD[, cumsum(forc)],
        by = .(id)
    ]

    return(pp_dt)
}


#' Prepare survival analysis data.
#' 
#' @param sos_dt SOS data.table containing at least `siteID`, `Date`, 
#' `Tmax`, `Tmin`, `Tmean`, `PhenoYear`, and `SOS`.
#' @param chill_forc_fun Function used to calculate daily and accumulated 
#' chilling and forcing.
#' @param ... Additional arguments passed to `chill_forc_fun`.
#' @return An list object containing necessary survival analysis data.
PrepareSurvData <- function(sos_dt, chill_forc_fun, ...) {
    # Person-period full
    pp_dt_full <- ToPersonPeriod(sos_dt)
    # Person-level table format
    pl_dt <- ToPersonLevel(pp_dt_full)
    # Calculate chilling and forcing
    pp_dt_full <- do.call(chill_forc_fun, args = list(pp_dt_full, ...))
    # Person-period format for survival analysis
    pp_dt <- ToPersonPeriodForSurv(pp_dt_full)

    # Caculate survival object
    ts <- survfit(Surv(doy, event) ~ 1, data = pl_dt)

    return(list(
        pp_dt_full = pp_dt_full,
        pl_dt = pl_dt,
        pp_dt = pp_dt,
        start_time = ts$time[1]
    ))
}



# ~ MCD12Q2 + Daymet ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mcd_dm_dt <- fread(file.path(gdir, "Data/ARD", "flux_mcd12q2_daymet_dt.csv"))
mcd_dm_dt[, Date := as.Date(Date)]

# merge site meta
site_meta <- fread(file.path(gdir, "Data/ARD", "sites_fluxnet2015_Tier1.csv"))
mcd_dm_dt <- merge(mcd_dm_dt, site_meta[, .(siteID, IGBP)], by = "siteID")

# filter out SOS > 255
mcd_dm_dt <- mcd_dm_dt[SOS < 255 & 
    IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), 
]



pdf(file.path(gdir, "Pipeline", "surv_fit_flux_mcd.pdf"), 
    width = 10, height = 4
)
fit_dt <- lapply(mcd_dm_dt[, unique(siteID)], function(x) {
    # x <- "US-NR1"

    site_dt <- mcd_dm_dt[siteID == x, ]

    if (nrow(site_dt) / 365 < 5) {
        return(NA)
    }
    
    site_li <- FormatDataForPhenoModel(site_dt)

    # Fit a TT model
    mod_tt <- FitTheModel(
        model = ThermalTimeModel,
        cost_fun = CostRMSE,
        data = site_li
    )

    tt_pars <- mod_tt$optim$par

    surv_data <- PrepareSurvData(site_dt, CalChillForc,
        Tc = 10, tc = -60,
        Tb = 2, t0 = 100
    )

    # Survival model fit
    afit1 <- glm(event ~ acc_forc * doy,
        family = "binomial", 
        data = surv_data$pp_dt[doy >= surv_data$start_time, ]
    )
    afit2 <- glm(event ~ acc_forc * doy + acc_chill,
        family = "binomial", 
        data = surv_data$pp_dt[doy >= surv_data$start_time, ]
    )
    # summary(afit2)

    ano <- anova(afit1, afit2)
    dif <- ano$Deviance[2]

    par(mfrow = c(1, 3))
    CompareXY_J(mod_tt$fitted_doy, site_li$transition_dates,
        main = "TT",
        xlab = "Est.",
        ylab = "Obs."
    )
    pred_dt1 <- PlotSurvGoodFit(afit1, surv_data$pp_dt_full, main = x)
    pred_dt2 <- PlotSurvGoodFit(afit2, surv_data$pp_dt_full, main = x)

    pred_dt1 <- SurvPredict(afit1, surv_data$pp_dt_full)
    pred_dt2 <- SurvPredict(afit2, surv_data$pp_dt_full)

AICc(site_li$transition_dates, mod_tt$fitted_doy, k = 3)
AICc(pred_dt1$SOS, pred_dt1$pred_SOS, k = 3)
AICc(pred_dt2$SOS, pred_dt2$pred_SOS, k = 4)

    # row <- data.table(
    #     siteID = x,
    #     RMSE = RMSE(pred_dt2$SOS, pred_dt2$pred_SOS),
    #     ano_dif = dif
    # )

    # return(row)
})

dev.off()







# ~ MCD12Q2 + flux_temperature ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# in: flux temperature
flux_temp <- fread(file.path(gdir, "Data/ARD", 
    "flux_measured_daily_temperature.csv"
))
flux_temp[, Date := as_date(Date)]

# in: flux mcd12q2 phenology
flux_dt <- fread(file.path(gdir, "Data/ARD/flux_mcd12q2_daymet_dt.csv"))
flux_dt[, Date := as_date(Date)]
# merge site meta
site_meta <- fread(file.path(gdir, "Data/ARD/sites_fluxnet2015_Tier1.csv"))
flux_dt <- merge(flux_dt, site_meta[, .(siteID, IGBP)], by = "siteID")
# filter out SOS > 255
flux_dt <- flux_dt[SOS < 255 & IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), ]


# Merge the two
flux_dt <- merge(flux_dt, flux_temp, by = c("siteID", "Date"), all.x = TRUE)

# Filter out site-years with measured temperature data that have missing values
flux_temp_dt <- by(flux_dt, flux_dt[, .(siteID, PhenoYear)], function(sy) {
    if (any(is.na(sy$Temp)) == FALSE) {
        return(sy)
    }
})
flux_temp_dt <- do.call(rbind, flux_temp_dt)
flux_temp_dt[, .N / 365, by = siteID]

# Filter out non-forest sites
flux_temp_dt <- flux_temp_dt[IGBP %in% c("MF", "DBF", "ENF", "EBF", "DNF"), .(
    siteID, Date,
    Tmax = NA, Tmin = NA, Tmean = Temp, SOS, PhenoYear
)]



fit_dt <- lapply(flux_temp_dt[, unique(siteID)], function(x) {
    # site_dt <- flux_temp_dt[siteID == "US-Ne1", unique(SOS)]
    site_dt <- flux_temp_dt[siteID == x, ]

    site_li <- FormatDataForPhenoModel(site_dt)

    # Fit a TT model
    mod_tt <- FitTheModel(
        model = ThermalTimeModel,
        cost_fun = CostRMSE,
        data = site_li
    )

    tt_pars <- mod_tt$optim$par

    # Person-period full
    pp_dt_full <- ToPersonPeriod(site_dt)
    # Person-level table format
    pl_dt <- ToPersonLevel(pp_dt_full)
    # Calculate chilling and forcing
    pp_dt_full <- CalChillForc(pp_dt_full, 
        Tc = 0, tc = -80, 
        Tb = tt_pars[2], t0 = tt_pars[1] - 110
    )
    # Person-period format for survival analysis
    pp_dt <- ToPersonPeriodForSurv(pp_dt_full)

    # Caculate survival object
    ts <- survfit(Surv(doy, event) ~ 1, data = pl_dt)

    # Survival model fit
    afit1 <- glm(event ~ acc_forc * doy,
        family = "binomial", data = pp_dt[doy >= ts$time[1], ]
    )
    afit2 <- glm(event ~ acc_forc * doy + acc_chill,
        family = "binomial", data = pp_dt[doy >= ts$time[1], ]
    )

    ano <- anova(afit1, afit2)
    dif <- ano$Deviance[2]

    # Predict raw hazard
    pred_dt <- SurvPredict(afit2, pp_dt_full)
    # CompareXY_J(pred_dt$pred_SOS, pred_dt$SOS, main = x)

    # readline(prompt = "enter to continue...")

    row <- data.table(
        siteID = x,
        RMSE = RMSE(pred_dt$SOS, pred_dt$pred_SOS),
        ano_dif = dif
    )

    return(row)
})
fit_dt <- do.call(rbind, fit_dt)

fit_dt[RMSE < 10 & ano_dif > 5, ]






# ~ HF & HB + ground_temperature ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hfhb_dt <- fread(file.path(gdir, "Data/ARD", "hb_hf_pheno_temp.csv"))[,
    Date := as_date(Date)
]

fit_dt <- lapply(hfhb_dt[, unique(siteID)], function(x) {
    # site_dt <- hfhb_dt[siteID == "US-Ne1", unique(SOS)]
    site_dt <- hfhb_dt[siteID == x, ]

    # Fit a TT model
    mod_tt <- FitTheModel(
        model = ThermalTimeModel,
        cost_fun = CostRMSE,
        data = site_li
    )

    tt_pars <- mod_tt$optim$par

    # Person-period full
    pp_dt_full <- ToPersonPeriod(site_dt)
    # Person-level table format
    pl_dt <- ToPersonLevel(pp_dt_full)
    # Calculate chilling and forcing
    pp_dt_full <- CalChillForc(pp_dt_full, 
        Tc = 0, tc = -80, 
        Tb = tt_pars[2], t0 = tt_pars[1] - 110
    )
    # Person-period format for survival analysis
    pp_dt <- ToPersonPeriodForSurv(pp_dt_full)

    # Caculate survival object
    ts <- survfit(Surv(doy, event) ~ 1, data = pl_dt)

    # Survival model fit
    afit1 <- glm(event ~ acc_forc * doy,
        family = "binomial", data = pp_dt[doy >= ts$time[1], ]
    )
    afit2 <- glm(event ~ acc_forc * doy + acc_chill,
        family = "binomial", data = pp_dt[doy >= ts$time[1], ]
    )

    ano <- anova(afit1, afit2)
    dif <- ano$Deviance[2]

    # Predict raw hazard
    pred_dt <- SurvPredict(afit2, pp_dt_full)
    CompareXY_J(pred_dt$pred_SOS, pred_dt$SOS, main = x)

    readline(prompt = "enter to continue...")

    row <- data.table(
        siteID = x,
        RMSE = RMSE(pred_dt$SOS, pred_dt$pred_SOS),
        ano_dif = dif
    )

    return(row)
})
fit_dt <- do.call(rbind, fit_dt)

fit_dt[RMSE < 10 & ano_dif > 5, ]
