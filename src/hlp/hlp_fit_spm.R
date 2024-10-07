#'******************************************************************************
#' Description: Helper functions to quickly fit spring phenology models.
#'******************************************************************************
source("src/hlp/hlp_temperature_season.R")
source("src/hlp/hlp_model_cost_functions.R")
source("src/hlp/hlp_spring_pheno_models.R")
source("src/hlp/hlp_chilling_calc.R")



#' Get spring phenology model parameter ranges.
#' 
#' @param model_name Model name abbreviation (e.g., TT, UN)
#' @return A data.table contains all parameters initial values and ranges 
#' of the model
#' @export
GetModelParRange <- function(model_name = NULL) {
    if (is.null(model_name)) {
        stop("Must provide a model name to get its parameter ranges")
    }

    par_file <- jsonlite::fromJSON("src/par_range.json")
    if (!model_name %in% names(par_file)) {
        stop("Provided model name is not in `par_range.json`")
    }

    mod_par <- par_file[[model_name]]
    
    par_dt <- data.table()
    par_dt$name <- names(mod_par)
    par_dt$init <- as.numeric(sapply(mod_par, "[[", "init"))
    par_dt$lower <- as.numeric(sapply(mod_par, "[[", "lower"))
    par_dt$upper <- as.numeric(sapply(mod_par, "[[", "upper"))

    return(par_dt)
}


#' Fit a spring phenology mdoel. 
#' The optimization method is fixed to `GenSA` because I found it's better than 
#' other methods in finding the optimal parameter values.
#' 
#' @param model The model function or string.
#' @param cost_fun The cost function or string.
#' @param data The data needed.
#' @param par The model parameters init values.
#' @param lower Lower boundaries of parameters.
#' @param upper Upper boundaries of parameters.
#' @param control The optimization algorithm controls. For example, set `maxit` 
#' to 500 would make the algorithm iterate 500 times maximum.
#' @return The model fit object.
FitTheModel <- function(model, cost_fun, data, par = NULL, 
    lower = -Inf, upper = Inf, control = list()
) {
    if (is.null(par) & (all(is.infinite(lower)) | all(is.infinite(upper)))) {
        # Get par ranges from the json file
        par_range <- GetModelParRange(
            model_name = as.character(substitute(model))
        )
        par <- par_range$init
        lower <- par_range$lower
        upper <- par_range$upper
    }

    est_par <- GenSA::GenSA(
        par = par,
        model = model,
        fn = cost_fun,
        data = data,
        lower = lower,
        upper = upper,
        control = control
    )
    est <- do.call(model, list(par = est_par$par, data = data))
    return(list(optim = est_par, fitted_doy = est))
}


#' Fit the CDSOM model.
#' 
#' @param dt A data table contains person-period data records.
#' @param x_names The column names in the data table that specify the names of
#' the columns to be used as predictors.
#' @param t0 The start counting date in DOY format. Default is -30.
#' @return A list contains the MCMC sampling result object and the data.table 
#' with fitted DOYs.
#' @export
FitCDSOMModel <- function(dt, x_names, t0 = -30) {
    require(R2jags)

    # Format JAGS data
    X <- as.matrix(
        data.table(
            intercept = 1,
            dt[, x_names, with = FALSE]
        )
    )
    Y <- (dt$event > 0) * 1

    data <- list(
        X = X,
        Y = Y,
        head_nodes = which(dt$doy %in% min(dt$doy):t0),
        main_nodes = which(!(dt$doy %in% min(dt$doy):t0)),
        n = nrow(X),
        np = ncol(X),
        hmax = 100,
        lambda = 1
    )

    model <- jags(
        model.file = textConnection(CDSOMModel()),
        data = data,
        parameters.to.save = c("beta", "yp", "h", "kappa"),
        n.chains = 1,
        n.iter = 10000,
        n.burnin = 3000
    )

    # The MCMC sampling result
    gibbs <- model$BUGSoutput$sims.list

    # Model fitted values
    compare <- data.table(
        PhenoYear = dt$PhenoYear,
        siteID = dt$siteID,
        yp = apply(gibbs$yp, 2, median),
        y = Y
    )

    onset <- compare[, .(
        onset_pred = min(c(Inf, which(yp == 1))),
        onset = min(c(Inf, which(y == 1)))
    ), by = .(siteID, PhenoYear)]

    onset <- onset[is.infinite(onset) == FALSE]

    return(list(
        gibbs = gibbs,
        fitted = onset
    ))
}


FitTT <- function(dt, ifplot = FALSE) {
    dt_li <- FormatDataForPhenoModel(dt)
    mod_tt <- FitTheModel(
        model = ThermalTimeModel,
        cost_fun = CostRMSE,
        data = dt_li
    )

    if (ifplot) {
        CompareXY_J(mod_tt$fitted_doy, dt_li$transition_dates,
            main = "TT"
        )
    }

    return(list(
        optim = mod_tt$optim,
        fitted_doy = mod_tt$fitted_doy,
        data_doy = dt_li$transition_dates
    ))
}


#' Fit `TTfixed` model
#'
#'
#' @param dt Data containing spring phenophase transition dates and the
#'   corresponding temperature records
#' @return The model fit object.
#' @export
FitTTfixed <- function(dt, t0, T_base, ifplot = FALSE) {
    dt_li <- FormatDataForPhenoModel(dt)

    est_par <- GenSA::GenSA(
        par = c(500),
        model = TTfix,
        fn = CostRMSE,
        data = dt_li,
        lower = c(0),
        upper = c(2000),
        t0 = t0,
        T_base = T_base,
        control = NULL
    )
    est <- do.call(TTfix, 
        list(par = est_par$par, data = dt_li, t0 = t0, T_base = T_base)
    )

    if (ifplot) {
        CompareXY_J(est, dt_li$transition_dates, main = "TT_fixed")
    }

    return(list(optim = est_par, 
        fitted_doy = est, 
        data_doy = dt_li$transition_dates)
    )
}



#' Fit Survvial models to evaluate chilling
#' 
#' @description 
#' Fit two survival models, one w/ chilling and one w/o chilling, then compare
#' model deviance to determine whether chilling is important.
#' 
#' @param dt A data.table contains `siteID`, `Date`, `Tmean`, `SOS`, and 
#'   `PhenoYear`
#' @param chill_start_threshold The threshold to determine when to start 
#'   counting chilling as well as the chilling temprature.
#' @param forc_start_threshold The threshold to determine when to start
#'   counting foring as well as the forcing temprature.
#' @return A data.table contains model fitting results.
FitSaChillComp <- function(dt, 
    chill_start_threshold = 0.5, forc_start_threshold = 0.05,
    ifplot = FALSE, ...
) {
    surv_data <- PrepareSurvDataPerSite(dt, 
        chill_start_threshold, forc_start_threshold
    )

    # Survival model fit
    afit1 <- glm(event ~ acc_forc * doy,
        family = "binomial",
        data = surv_data$pp_dt_for_mod
    )

    afit2 <- glm(event ~ acc_forc * doy + acc_chill,
        family = "binomial",
        data = surv_data$pp_dt_for_mod
    )

    null_rmse <- RMSE(surv_data$pl_dt$SOS, mean(surv_data$pl_dt$SOS))
    afit1_pred <- SurvPredict(afit1, surv_data$pp_dt_for_pred)
    afit1_rmse <- RMSE(afit1_pred$SOS, afit1_pred$pred_SOS)
    afit2_pred <- SurvPredict(afit2, surv_data$pp_dt_for_pred)
    afit2_rmse <- RMSE(surv_data$pl_dt$SOS, afit2_pred$pred_SOS)

    if (ifplot) {
        null_1 <- PlotSurvGoodFit(afit1_pred, ...)
        null_2 <- PlotSurvGoodFit(afit2_pred, ...)
    }

    dev_diff <- anova(afit1, afit2)$Deviance[2]

    afit1_aicc <- AICc(afit1_pred$SOS, afit1_pred$pred_SOS, length(coef(afit1)))
    afit2_aicc <- AICc(afit2_pred$SOS, afit2_pred$pred_SOS, length(coef(afit2)))

    res_dt <- data.table(
        siteID = surv_data$siteID,
        NoChill_RMSE = afit1_rmse,
        Chill_rmse = afit2_rmse,
        NoChill_AICc = afit1_aicc,
        Chill_AICc = afit2_aicc,
        Null_RMSE = null_rmse,
        Diff_dev = dev_diff,
        Chill_th = chill_start_threshold,
        Forc_th = forc_start_threshold
    )

    return(res_dt)
}


FitSaChillComp2 <- function(surv_data, ifplot = FALSE) {
    # Survival model fit
    afit1 <- glm(event ~ acc_forc * doy,
        family = "binomial",
        data = surv_data$pp_dt_for_mod
    )

    afit2 <- glm(event ~ acc_forc * doy + acc_chill,
        family = "binomial",
        data = surv_data$pp_dt_for_mod
    )

    null_rmse <- RMSE(surv_data$pl_dt$SOS, mean(surv_data$pl_dt$SOS))
    afit1_pred <- SurvPredict(afit1, surv_data$pp_dt_for_pred)
    afit1_rmse <- RMSE(afit1_pred$SOS, afit1_pred$pred_SOS)
    afit2_pred <- SurvPredict(afit2, surv_data$pp_dt_for_pred)
    afit2_rmse <- RMSE(surv_data$pl_dt$SOS, afit2_pred$pred_SOS)

    if (ifplot) {
        null_1 <- PlotSurvGoodFit(afit1_pred)
        null_2 <- PlotSurvGoodFit(afit2_pred)
    }

    dev_diff <- anova(afit1, afit2)$Deviance[2]

    afit1_aicc <- AICc(afit1_pred$SOS, afit1_pred$pred_SOS, length(coef(afit1)))
    afit2_aicc <- AICc(afit2_pred$SOS, afit2_pred$pred_SOS, length(coef(afit2)))

    res_dt <- data.table(
        siteID = surv_data$siteID,
        NoChill_RMSE = afit1_rmse,
        Chill_rmse = afit2_rmse,
        NoChill_AICc = afit1_aicc,
        Chill_AICc = afit2_aicc,
        Null_RMSE = null_rmse,
        Diff_dev = dev_diff
    )

    return(res_dt)
}



#' Calculate chilling and forcing for person-period data format
#' 
#' @param pp_dt The person-period data.table.
#' @param cfun Chilling calculation function.
#' @param t0,Tb The start counting day and base temperature for forcing.
#' @param tc The start counting day for chilling.
#' @param ... The additional parameters passed to the `cfun`.
#' @return A person-period data.table with chilling and forcing values.
#' @export
CalCF <- function(pp_dt, cfun = NULL, t0, Tb, tc, ...) {
    # Calculate chilling
    pp_dt[, ":="(chill = 0, acc_chill = 0)]
    pp_dt[doy > tc,
        chill := .SD[, do.call(cfun, list(Tmean, ...))],
        by = .(id)
    ]
    # Accumulated chilling
    pp_dt[doy > tc,
        acc_chill := .SD[, cumsum(chill)],
        by = .(id)
    ]
    # pp_dt[doy <= tc, ":=" (chill = NA, acc_chill = NA)]

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
    # pp_dt[doy <= t0, ":="(forc = NA, acc_forc = NA)]

    return(pp_dt)
}




#' Format survival data using `ChillFix()` to calculate chilling
#' 
#' @description
#' Estimate temperature seasonality based on each site and calculate chilling
#' and forcing separately for each site. This type of chilling calculation may
#' not well represent all cases as there are places where the fixed temperature
#' thresholds rarely work.
#' 
#' @param dt A person-period data.table.
#' @param tc_thresh Threshold used to extract chilling.
#' @param t0_thresh Threshold used to extract forcing.
#' @param method Default "a", will be passed to the `ChillFix()` function.
#' @return Survival data object, ready for fitting survival models.
#' @export
FmSurvDataChillFix <- function(dt, tc_thresh, t0_thresh, method = "a") {
    site_dt <- by(dt, dt$siteID, function(id_dt) {
        # Retrieve temperature seasonality
        site_season <- TempSeason(id_dt,
            tc_thresh,
            t0_thresh
        )
        # Person-period full
        pp_dt_full <- ToPersonPeriod(id_dt)
        # Person-level table format
        pl_dt <- ToPersonLevel(pp_dt_full)
        
        # Calculate chilling and forcing
        Tb <- max(0, site_season$spring_min + 
            site_season$spring_amp * t0_thresh
        )
        pp_dt_full <- CalCF(
            pp_dt_full, cfun = ChillFix, 
            t0 = site_season$t0forc,
            Tb = Tb,
            tc = site_season$t0chill,
            method = method
        )
        
        # Person-period format for survival analysis
        pp_dt <- ToPersonPeriodForSurv(pp_dt_full)

        surv_data <- list(
            pp_dt_full = pp_dt_full,
            pl_dt = pl_dt,
            pp_dt = pp_dt
        )

        return(surv_data)
    })

    surv_data <- list(
        pp_dt_full = lapply(site_dt, "[[", "pp_dt_full") %>%
            do.call(rbind, .),
        pl_dt = lapply(site_dt, "[[", "pl_dt") %>%
            do.call(rbind, .),
        pp_dt = lapply(site_dt, "[[", "pp_dt") %>%
            do.call(rbind, .)
    )

    # Sort
    setorder(surv_data$pp_dt_full, Date)
    setorder(surv_data$pl_dt, Date)
    setorder(surv_data$pp_dt, Date)

    ids <- unique(surv_data$pp_dt$siteID)
    surv_data$siteID <- ifelse(length(ids) == 1, ids, paste(ids, sep = "_"))

    ts <- survfit(Surv(doy, event) ~ 1, data = surv_data$pl_dt)
    start_time <- ts$time[1]

    # Use the calculated start time to subset all data
    surv_data$pp_dt_for_pred <- surv_data$pp_dt_full[doy >= start_time]
    surv_data$pp_dt_for_mod <- surv_data$pp_dt[doy >= start_time]

    return(surv_data)
}


#' Format survival data using `Trichill()` function
#' 
#' @description
#' Estimate temperature seasonality based on each site and calculate chilling
#' and forcing separately for each site. The `Topt`, `Tmin`, `Tmax` will all
#' be estimated instead of fixed.
#' 
#' @param dt A person-period data.table.
#' @param tc_thresh Threshold used to extract chilling.
#' @param tc_min_thresh,tc_max_thresh Minimum and Maximum chilling temperature 
#'   relative to autumn temperature amplitude.
#' @param t0_thresh Threshold used to extract forcing.
#' @return Survival data object, ready for fitting survival models.
#' @export
FmSurvDataTriChill <- function(dt, tc_thresh, tc_min_thresh, tc_max_thresh, 
    t0_thresh
) {
    site_dt <- by(dt, dt$siteID, function(id_dt) {
        # Retrieve temperature seasonality
        site_season <- TempSeason(id_dt,
            chill_start_threshold = tc_thresh,
            forc_start_threshold = t0_thresh
        )

        # Person-period full
        pp_dt_full <- ToPersonPeriod(id_dt)
        # Person-level table format
        pl_dt <- ToPersonLevel(pp_dt_full)
        # Calculate chilling and forcing
        Tmin <- site_season$autumn_amp * tc_min_thresh
        Tmax <- site_season$autumn_amp * tc_max_thresh
        Topt <- (Tmax + Tmin) / 2
        pp_dt_full <- CalCF(
            pp_dt_full, cfun = TriChill, 
            t0 = site_season$t0forc,
            Tb = site_season$spring_amp * t0_thresh,
            tc = site_season$t0chill,
            Tmin = Tmin,
            Tmax = Tmax,
            Topt = Topt
        )
        # Person-period format for survival analysis
        pp_dt <- ToPersonPeriodForSurv(pp_dt_full)

        surv_data <- list(
            pp_dt_full = pp_dt_full,
            pl_dt = pl_dt,
            pp_dt = pp_dt
        )

        return(surv_data)
    })

    surv_data <- list(
        pp_dt_full = lapply(site_dt, "[[", "pp_dt_full") %>%
            do.call(rbind, .),
        pl_dt = lapply(site_dt, "[[", "pl_dt") %>%
            do.call(rbind, .),
        pp_dt = lapply(site_dt, "[[", "pp_dt") %>%
            do.call(rbind, .)
    )

    # Sort 
    setorder(surv_data$pp_dt_full, Date)
    setorder(surv_data$pl_dt, Date)
    setorder(surv_data$pp_dt, Date)

    ids <- unique(surv_data$pp_dt$siteID)
    surv_data$siteID <- ifelse(length(ids) == 1, ids, paste(ids, collapse = "&"))

    ts <- survfit(Surv(doy, event) ~ 1, data = surv_data$pl_dt)
    start_time <- ts$time[1]

    # Use the calculated start time to subset all data
    surv_data$pp_dt_for_pred <- surv_data$pp_dt_full[doy >= start_time]
    surv_data$pp_dt_for_mod <- surv_data$pp_dt[doy >= start_time]

    return(surv_data)
}
