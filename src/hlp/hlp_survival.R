#'******************************************************************************
#' Description: Base functions for survival analysis.
#'******************************************************************************


#' Complementary log-log function
#'
#' @param x Numeric vector.
#' @return clog-log transformed vector.
cloglog <- function(x) {
    cloglog_val <- log(-log(1 - x))
    return(cloglog_val)
}

#' Inverse clog-log function
#'
#' @param cloglog_val Numeric vector of values in clog-log scale.
#' @return Transformed vector.
icloglog <- function(cloglog_val) {
    prob <- 1 - 1 / exp(exp(cloglog_val))
    return(prob)
}

#' Logistic function
#'
#' @param x Numeric vector.
#' @return logit transformed vector.
logit <- function(x) {
    return(log(x / (1 - x)))
}

#' Inverse logistic function
#'
#' @param x Logitic values.
#' @return Transformed vector.
ilogit <- function(x) {
    return(1 / (1 + exp(-x)))
}


#' Convert pheno data list format to data table format
#'
#' @param data_li Spring phenology model data format list.
#' @return SOS data in data table format.
PhenoModelDataToDT <- function(data_li) {
    dt <- lapply(seq_along(data_li$site), function(idx) {
        siteID <- data_li$site[idx]
        PhenoYear <- data_li$year[idx]
        Date <- as_date(paste0(PhenoYear, "-01-01")) + data_li$doy
        Tmean <- data_li$Ti[, idx]
        Tmax <- data_li$T_max[, idx]
        Tmin <- data_li$T_min[, idx]
        SOS <- data_li$transition_dates[idx]

        sy_dt <- data.table(siteID, PhenoYear, SOS, Date, Tmean, Tmax, Tmin)

        return(sy_dt)
    }) %>%
        do.call(rbind, .)

    return(dt)
}

#' Convert SOS data.table format to person period format
#'
#' @param dt SOS data in data.table format.
#' @return Person period data format.
ToPersonPeriod <- function(dt) {
    dt$Date <- as_date(dt$Date)
    doy <- ifelse(year(dt$Date) < dt$PhenoYear,
        dt$Date - as.Date(paste0(dt$PhenoYear, "-01-01")),
        dt$Date - as.Date(paste0(dt$PhenoYear, "-01-01")) + 1
    )
    dt$doy <- doy
    dt$event <- ifelse(dt$SOS <= dt$doy, 1, 0)

    # Make an `id` column
    dt$id <- paste(dt$siteID, dt$PhenoYear, sep = "_")
    setcolorder(dt, "id")

    return(dt)
}

#' Convert full SOS data in person-period format to person-level format
#'
#' @param pp_dt_full SOS full data in person-period format.
#' @return Person-level period data.
ToPersonLevel <- function(pp_dt_full) {
    pl_dt <- pp_dt_full[,
        .SD[round(SOS) == doy, .(Date, SOS, doy, event)],
        by = .(id)
    ]
    pl_dt$event <- 1
    pl_dt$SOS <- round(pl_dt$SOS)

    return(pl_dt)
}

#' Retrieve person-period data for survival analysis, i.e., rows after event
#' happening will be removed
#'
#' @param pp_dt_full Full person-period data format.
#' @return Person-period data format for survival analysis.
ToPersonPeriodForSurv <- function(pp_dt_full) {
    # Person-period format for survival analysis
    pp_dt <- by(pp_dt_full, pp_dt_full$id, function(site) {
        sos <- site[event > 0][1]
        pre <- site[event < 1]

        res_dt <- rbind(pre, sos)

        return(res_dt)
    }) %>%
        do.call(rbind, .)

    return(pp_dt)
}


#' Survival model prediction
#'
#' @param mod The survival model.
#' @param pp_dt_for_pred Full person-period data table.
#' @return Prediction data.table
SurvPredict <- function(mod, pp_dt_for_pred) {
    # Predict raw hazard
    pp_dt_for_pred$pred_h <- predict(mod, 
        newdata = pp_dt_for_pred, 
        type = "response"
    )

    pred_dt <- by(pp_dt_for_pred, pp_dt_for_pred$id, function(id_dt) {
        n <- nrow(id_dt)
        h <- id_dt$pred_h
        pmf <- h
        for (t in 2:n) {
            pmf[t] <- h[t] * prod(1 - h[1:(t - 1)], na.rm = TRUE)
        }
        pmf[n] <- 1 - sum(pmf[1:(n - 1)], na.rm = TRUE)

        pmf[is.na(pmf)] <- 0

        # plot(pmf)

        cdf <- cumsum(pmf)
        # plot(cdf)

        # Compute summaries of the survival distribution
        q50 <- which.min(abs(cdf - 0.5))
        q05 <- max(which(cdf < 0.05), 0)
        q95 <- min(which(cdf > 0.95))
        mn <- sum((1:n) * pmf)

        # Make an ugly plot
        # plot(pmf)
        # abline(v = mn, col = 2)
        # abline(v = q50, col = 3)
        # abline(v = q05, col = 4)
        # abline(v = q95, col = 4)
        # legend("topright",
        #     c("Mean", "Median", "90% interval"),
        #     lty = 1,
        #     col = 2:4, bty = "n"
        # )

        return(data.table(
            siteID = id_dt$siteID[1],
            PhenoYear = id_dt$PhenoYear[1],
            SOS = id_dt$SOS[1],
            pred_SOS = id_dt[q50, doy],
            pred_lwr = ifelse(q05 == 0, NA, id_dt[q05, doy]),
            pred_upr = id_dt[q95, doy]
        ))
    })

    pred_dt <- do.call(rbind, pred_dt)

    return(pred_dt)
}


#' Prepare survival analysis data
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
    start_time <- ts$time[1]

    # Use the calculated start time to subset all data
    pp_dt_for_pred <- pp_dt_full[doy >= start_time - 30]
    pp_dt_for_mod <- pp_dt[doy >= start_time]

    return(list(
        pp_dt_full = pp_dt_full,
        pl_dt = pl_dt,
        pp_dt_for_mod = pp_dt_for_mod,
        pp_dt_for_pred = pp_dt_for_pred,
        start_time = ts$time[1]
    ))
}


#' Plot survival goodness-of-fit
#' 
#' @param pred_dt The survival model predicted data.table using `SurvPredict`.
#' @param ifcenter Indicates if the plot will be centered to anomalies 
#' @param rg Range of x and y axes.
#'   by siteID.
#' @return NULL
PlotSurvGoodFit <- function(pred_dt, ifcenter = FALSE, rg = NULL, ...) {
    require(lmodel2)

    x <- pred_dt$pred_SOS
    y <- pred_dt$SOS

    if (ifcenter == TRUE) {
        x <- pred_dt[, scale(pred_SOS, scale = FALSE), by = siteID]$V1
        y <- pred_dt[, scale(SOS, scale = FALSE), by = siteID]$V1
    }

    if (is.null(rg)) {
        # rg <- range(x, y, pred_dt$pred_lwr, na.rm = TRUE)
        rg <- range(x, y, na.rm = TRUE)
    }
    plot(NA, xlim = rg, ylim = rg, mgp = c(1.5, 0.5, 0), 
        xlab = "Est.", ylab = "Obs.",
        bty = "L",
        main = "Survival Model",
        ...
    )

    # Error bars
    # arrows(
    #     x0 = pred_dt$pred_lwr, x1 = pred_dt$pred_upr,
    #     y0 = pred_dt$SOS, y1 = pred_dt$SOS,
    #     angle = 90, code = 3, length = 0.02,
    #     col = "grey"
    # )
    # Data points
    points(x, y, pch = 16, cex = 0.8, col = adjustcolor("grey", 0.9), ...)

    # draw a 1:1 line
    abline(a = 0, b = 1, lty = 2)

    rmse <- NA

    fit <- tryCatch({
        suppressMessages(
            lmodel2(y ~ x, data.frame(x = x, y = y), "interval", "interval")
        )
    }, error = function(x) return(NULL))

    if (!is.null(fit) & any(is.na(coef(fit))) == FALSE) {
        fit_res <- data.table(fit$regression.results)
        inter_slp <- fit_res[Method == "RMA", .(Intercept, Slope)]
        abline(a = inter_slp$Intercept, b = inter_slp$Slope, col = "red")
        # # get equation
        # cf <- round(coef(fit), 2)
        # eq <- paste(
        #     "y=",
        #     ifelse(sign(cf[2]) == 1, "", "-"),
        #     format(abs(cf[2]), nsmall = 2),
        #     "x",
        #     ifelse(sign(cf[1]) == 1, "+", "-"),
        #     format(abs(cf[1]), nsmall = 2),
        #     sep = ""
        # )
        # # equation <- bquote(.(eq) ~ R^2 == .(r_sqr))

        # get stats
        r_sqr <- round(fit$rsquare, 2)
        r_sqr_str <- bquote(R^2 == .(r_sqr))

        rmse <- round(sqrt(mean((y - x)^2, na.rm = TRUE)), 2)
        rmse_str <- bquote(RMSE == .(rmse))

        p_val <- fit$P.param
        # f <- summary(fit)$fstatistic
        # p_val <- pf(f[1], f[2], f[3], lower.tail = FALSE)
        p_val_str <- paste(
            "p-value: ",
            formatC(p_val, format = "e", digits = 2)
        )

        # rmse <- round(sqrt(mean((y - x)^2, , na.rm = TRUE)), 2)
        # rmse_str <- bquote(RMSE == .(rmse))

        # legend("bottomright", legend = eq, bty = "n", y.intersp = 7)
        legend("bottomright", legend = r_sqr_str, bty = "n", y.intersp = 5)
        legend("bottomright", legend = rmse_str, bty = "n", y.intersp = 3)
        legend("bottomright", legend = p_val_str, bty = "n")
    }
}


PrepareSurvDataPerSite <- function(dt, 
    chill_start_threshold = 0.5, forc_start_threshold = 0.05
) {
    site_dt <- by(dt, dt$siteID, function(id_dt) {
        # Retrieve temperature seasonality
        site_season <- TempSeason(id_dt,
            chill_start_threshold = chill_start_threshold,
            forc_start_threshold = forc_start_threshold
        )
        # Prepare survival data format
        # Cannot directly use the `PrepareSurvData()` function b/c it removes
        #   records before the first SOS date
        # Person-period full
        pp_dt_full <- ToPersonPeriod(id_dt)
        # Person-level table format
        pl_dt <- ToPersonLevel(pp_dt_full)
        # Calculate chilling and forcing
        pp_dt_full <- do.call(CalChillForc,
            args = list(pp_dt_full,
                Tc = site_season$autumn_amp * chill_start_threshold,
                tc = site_season$t0chill,
                Tb = site_season$spring_amp * forc_start_threshold,
                t0 = site_season$t0forc
            )
        )
        # pp_dt_full <- CalTriChilForc(
        #     pp_dt_full, 
        #     Topt = site_season$autumn_amp * chill_start_threshold,
        #     Tmin = site_season$autumn_amp * 0.1,
        #     Tmax = site_season$autumn_amp * 0.8
        # )
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


CalTriChilForc <- function(pp_dt, 
    Topt, Tmin, Tmax,
    Tb, t0
) {
    # Calculate chilling
    # Let's use chill days for now
    pp_dt[, ":=" (chill = 0, acc_chill = 0)]
    pp_dt[doy > tc,
        chill := .SD[, TriChill(Tmean, Topt, Tmin, Tmax)],
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

    return(pp_dt)
}

#' Calculate photoperiod (i.e., daylength) variables including daily daylength
#' and accumulated daylength.
#' 
#' @param pp_dt Person-period data table.
#' @param doy_colname DOY column name, default is "doy".
#' @param lat_colname Latitude column name.
#'
#' @return Person-period data.table with photoperiod variables as columns.
CalDaylength <- function(pp_dt, doy_colname = "doy", lat_colname) {
    require(geosphere)

    if (length(lat_colname) != 1 && 
        lat_colname %in% colnames(pp_dt) == FALSE
    ) {
        stop("No latitude column found!")
    }
    
    if (length(doy_colname) != 1 && 
        doy_colname %in% colnames(pp_dt) == FALSE
    ) {
        stop("No DOY column found!")
    }
    
    pp_dt[, daylen := daylength(get(lat_colname), get(doy_colname))]
    pp_dt[get(doy_colname) > 0, 
        acc_daylen := .SD[, cumsum(daylen)],
        by = .(id)
    ]

    return(pp_dt)
}


#' Fit a survival model to a point time series
#' 
#' @param pt_dt The point time series w/ continuous temperatures.
#' @param param The `tc`, `t0`, and `mtd` parameters used to calculate chilling.
#' @param equ The survial model formula.
#'
#' @return A list containing formatted survival data, model fit, and prediction.
FitSaFromTS <- function(pt_dt, param, equ) {
    surv_data <- FmSurvDataChillFix(pt_dt, param$tc, param$t0, param$mtd)
    afit <- glm(as.formula(equ),
        family = "binomial",
        data = surv_data$pp_dt_for_mod
    )
    afit_pred <- SurvPredict(afit, surv_data$pp_dt_for_pred)
    
    return(list(
        surv_data = surv_data,
        fit = afit,
        pred = afit_pred
    ))
}
