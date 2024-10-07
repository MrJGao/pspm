#'******************************************************************************
#' Description: This file contains useful functions and variables for the entire 
#' project.
#'******************************************************************************
library(data.table)
library(magrittr)

LoadHelpers <- function() {
   tmp <- list.files("src/hlp", full.names = TRUE)
   invisible(lapply(tmp, source))
}


# ~ API ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse a file containing some credentials or private information 
# (e.g., AppEEARS account) that I don't want them to be public
ParsePrivate <- function(attr_name = NULL) {
    doc <- jsonlite::fromJSON(file.path(gdir, "private.json"))
    if(!is.null(attr_name)) {
        return(doc[[tolower(attr_name)]])
    }
    return(doc)
}


# ~ Helper functions ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Convert Kelvin to Celsius
#' 
#' @param k Temperature in Kelvin.
#' @return Temperature in Celsius.
K2C <- function(k) {
    return(k - 273.15)
}

#' Two way `setdiff`
SetDiff <- function(a, b) {
    diff <- setdiff(union(a, b), intersect(a, b))

    return(diff)
}

# Custom AIC function which accepts loess regressions
# Copied from `phenor` package
AICc <- function(measured, predicted, k) {
    # calculate number of observations
    n <- length(measured)

    # calculatue residual sum of squares
    RSS <- sum((measured - predicted)^2)

    # AIC
    AIC <- 2 * k + n * log(RSS / n)

    # AICc
    AICc <- AIC + (2 * k * (k + 1)) / (n - k - 1)

    # return both AIC
    return(list(
        "AIC" = AIC,
        "AICc" = AICc
    ))
}


# read in NLCD legend
LoadNlcdLegend <- function() {
    nlcd_legend <- read.csv(file.path(gdir, "Data", "ARD", "nlcd_legend.csv"))
    return(nlcd_legend)
}


#' Format data table to a list object that is more convenient for fitting 
#' sping phenology models
#' @param dt A data table containing `siteID`, `PhenoYear`, `SOS`, `Tmean`, 
#' `Tmax`, and `Tmin`
#' @return A list object.
#' @export
FormatDataForPhenoModel <- function(dt, doy = c((-110):(-1), 1:255)) {
    # get pheno obs
    pheno_dt <- unique(dt[, .(siteID, PhenoYear, SOS)])

    # get daily mean temperature matrix
    temp_mean_mat <- apply(pheno_dt, 1, function(sy) {
        temp <- dt[siteID == sy[["siteID"]] & PhenoYear == sy[["PhenoYear"]], ]
        if (nrow(temp) != 365) {
            return(rep(NA, 365))
        }
        return(temp$Tmean)
    })

    # get daily min temperature matrix
    temp_min_mat <- apply(pheno_dt, 1, function(sy) {
        temp <- dt[siteID == sy[["siteID"]] & PhenoYear == sy[["PhenoYear"]], ]
        if (nrow(temp) != 365) {
            return(rep(NA, 365))
        }
        return(temp$Tmin)
    })

    # get daily max temperature matrix
    temp_max_mat <- apply(pheno_dt, 1, function(sy) {
        temp <- dt[siteID == sy[["siteID"]] & PhenoYear == sy[["PhenoYear"]], ]
        if (nrow(temp) != 365) {
            return(rep(NA, 365))
        }
        return(temp$Tmax)
    })

    # The result list
    li <- list()
    li$site <- pheno_dt$siteID
    li$year <- pheno_dt$PhenoYear
    li$transition_dates <- pheno_dt$SOS
    li$doy <- doy
    li$Ti <- temp_mean_mat
    li$T_min <- temp_min_mat
    li$T_max <- temp_max_mat

    # Drop site-years with NA temperature
    na_idx <- apply(li$Ti, 2, function(ti) {
        return(any(is.na(ti)))
    })

    if (any(na_idx)) {
        warning("Some site years have NA in temperature, drop them!")

        li$site <- li$site[!na_idx]
        li$year <- li$year[!na_idx]
        li$transition_dates <- li$transition_dates[!na_idx]
        li$Ti <- li$Ti[, !na_idx]
        li$T_min <- li$T_min[, !na_idx]
        li$T_max <- li$T_max[, !na_idx]

    }
    
    return(li)
}


# Convert the formatted list back to data.table
Conv2Dt <- function(data_li) {
    # data_li <- readRDS(file.path(hpc_local_dir, "Pipeline", "flux_li.Rds"))
    # data_dt <- flux_dt

    data_dt <- data.table()
    
    data_dt <- lapply(1:length(data_li$site), function(i) {
        site <- data_li$site[i]
        year <- data_li$year[i]
        sos <- data_li$transition_dates[i]
        Date <- as_date(paste0(year, "-01-01")) + data_li$doy
        Tmean <- data_li$Ti[, i]
        Tmin <- data_li$T_min[, i]
        Tmax <- data_li$T_max[, i]

        return(data.table(
            siteID = rep(site, 365),
            Date = Date,
            Tmax = Tmax,
            Tmin = Tmin,
            Tmean = Tmean,
            PhenoYear = year,
            SOS = sos
        ))
    })
    
    data_dt <- do.call(rbind, data_dt)

    return(data_dt)
}


#' Get the `idx` specified site years data from the whole `data_li`.
#' 
#' @param data_li The formatted data list.
#' @param idx The interested site years index.
#' @return Result data list with only the interested site years.
GetDataLiIndex <- function(data_li, idx) {
    li <- list()
    li$site <- data_li$site[idx]
    li$year <- data_li$year[idx]
    li$transition_dates <- data_li$transition_dates[idx]
    li$doy <- data_li$doy
    # Keep the dimensions of below variables
    li$Ti <- array(data_li$Ti[, idx], dim = c(365, length(idx)))
    li$T_min <- array(data_li$T_min[, idx], dim = c(365, length(idx)))
    li$T_max <- array(data_li$T_max[, idx], dim = c(365, length(idx)))

    return(li)
}


#' Get the p-value of the linear regression model.
#'
#' @param mod_lm Linear regression model object.
#' @return Numeric. p-value.
GetLmPval <- function(mod_lm) {
    f <- summary(mod_lm)$fstatistic
    pval <- pf(f[1], f[2], f[3], lower.tail = FALSE)

    return(pval)
}

#' Compute RMSE value for two vectors
#'
#' @param x,y Vector x and y.
#' @return RMSE value.
RMSE <- function(x, y, na_rm = TRUE) {
    rmse <- sqrt(mean((x - y)^2, na.rm = na_rm))
    return(rmse)
}



#' Convert day of year to date.
#'
#' @param doy Day of year.
#' @param year Year.
#' @return Date format variable
#' @export
DoyToDate <- function(doy, year) {
    date_str <- as.Date(doy, origin = paste0(year, "-01-01"))
    return(date_str)
}

# ~ Survival analysis related functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# #' Calculate AICc and AIC for comparing models
# #' 
# #' @description
# #' The default mode computes AICc as it is better for small sample sizes and 
# #' converges to AIC when the sample size is large.
# #' 
# #' @param obs Observed data.
# #' @param est Model estimated data.
# #' @param k Number of parameters in the model
# #' @param type Default is `AICc`, it can be changed to `AIC`.
# #' @return Commputed AICc (or AIC) value.
# #' @export
# AICc <- function(obs, est, k, type = "AICc") {
#     n <- length(obs)
#     RSS <- sum((obs - est)^2)

#     aic <- 2 * k + n * log(RSS / n)
#     aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)

#     if (tolower(type) == "aic") {
#         return(aic)
#     }

#     return(aicc)
# }

# ~ Base visualization functions ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Compare 2 vectors with a linear equation legend
CompareXY_J <- function(x, y = NULL, ...) {
    if (is.null(y)) {
        tmp <- x[, 1]
        y <- unlist(x[, 2])
        x <- unlist(tmp)
    }
    range <- range(x, y, na.rm = TRUE)
    plot(x, y, pch = 16, xlim = range, ylim = range, mgp = c(1.5, 0.5, 0), ...)
    # draw a 1:1 line
    abline(a = 0, b = 1, lty = 2)

    fit <- lm(y ~ x)
    if (any(is.na(coef(fit))) == FALSE) {
        abline(fit, col = "red")
        # get equation
        cf <- round(coef(fit), 2)
        eq <- paste(
            "y=", 
            ifelse(sign(cf[2]) == 1, "", "-"), 
            format(abs(cf[2]), nsmall = 2), 
            "x",
            ifelse(sign(cf[1]) == 1, "+", "-"), 
            format(abs(cf[1]), nsmall = 2),
            sep = ""
        )
        # equation <- bquote(.(eq) ~ R^2 == .(r_sqr))

        # get stats
        r_sqr <- round(summary(fit)$r.squared, 2)
        r_sqr_str <- bquote(R^2 == .(r_sqr))

        f <- summary(fit)$fstatistic
        p_val <- pf(f[1], f[2], f[3], lower.tail = FALSE)
        p_val_str <- paste(
            "p-value: ", 
            formatC(p_val, format = "e", digits = 2)
        )

        rmse <- round(sqrt(mean((y - x)^2, , na.rm = TRUE)), 2)
        rmse_str <- bquote(RMSE == .(rmse))

        legend("bottomright", legend = eq, bty = "n", y.intersp = 7)
        legend("bottomright", legend = r_sqr_str, bty = "n", y.intersp = 5)
        legend("bottomright", legend = rmse_str, bty = "n", y.intersp = 3)
        legend("bottomright", legend = p_val_str, bty = "n")
    }
}

#' Plot stretched raster image with a nice color and legend
#' @example PlotStretchLegend(r, viridis(200))
PlotStretchLegend <- function(r, colorRamp, digits = 0, labSize = 1, ...) {
    pal <- colorRampPalette(colorRamp)
    qs <- quantile(r, c(0, 0.02, 0.98, 1))
    if (qs[1] == qs[2]) {
        qs[1] <- qs[1] - 1e-6
    }
    if (qs[3] == qs[4]) {
        qs[4] <- qs[4] + 1e-6
    }
    
    r_breaks <- c(qs[1], seq(qs[2], qs[3], len = 255), qs[4])
    plot(r, col = pal(length(r_breaks) - 1), breaks = r_breaks, axes = F, box = F, legend = F, ...)
    # add a reasonable legend
    legend_at <- round(seq(r_breaks[2], r_breaks[length(r_breaks) - 1], len = 7), digits)
    legend_labels <- c(
        paste("<", legend_at[1]),
        as.character(legend_at[2:(length(legend_at) - 1)]),
        paste(">", legend_at[length(legend_at)])
    )
    plot(raster(matrix(legend_at)),
        legend.only = T, col = pal(length(r_breaks) - 1),
        axis.args = list(at = legend_at, labels = legend_labels, cex.axis = labSize)
    )
}

#' Plot the density and the 95% / 50% CI of a variable
PlotDensity <- function(x) {
    plot(density(x))
    abline(v = quantile(x, c(0.025, 0.975)), lty = 1)
    abline(v = quantile(x, c(0.25, 0.75)), lty = 2)
}


#' Add the linear regression line and CI as well as stats to current plot. Only
#' support simple linear regression with a single intercept and a single slope.
#' 
#' @param mod_lm The `lm` object.
#' @param text_pos Text position. Default "bottomright".
#' @param text_col Text color. Default black.
#' @param line_col Line color. Default black.
#' @param ci_col CI color, Default transparent gray.
#' @param ifroundp Logical, indicates whether the p-value should be rounded to 
#' 2 digits.
#' @return Add linear regression elements to current plot.
LmLines <- function(mod_lm, text_pos = "bottomright", text_col = "black", 
    line_col = "black", ci_col = adjustcolor("grey", alpha.f = 0.5),
    ifroundp = FALSE
) {
    stopifnot("The input model is not a linear regression object" = 
        class(mod_lm) == "lm"
    )

    # Draw line
    abline(mod_lm, col = line_col, lwd = 1.5)
    
    # Compute and draw CI
    x_range <- grconvertX(c(0, 1), "npc")
    x_vec <- seq(x_range[1], x_range[2], len = 100)
    newdf <- data.frame(x = x_vec)
    colnames(newdf) <- all.vars(mod_lm$call$formula)[2]
    y_pred <- predict(mod_lm, 
        newdata = newdf, 
        interval = "confidence"
    )
    polygon(c(x_vec, rev(x_vec)), c(y_pred[, "lwr"], rev(y_pred[, "upr"])),
        col = ci_col, border = NA
    )

    # Get regression stats

    # Equation
    cf <- round(coef(mod_lm), 3)
    eq <- paste(
        "y=",
        ifelse(sign(cf[2]) == 1, "", "-"),
        format(abs(cf[2]), nsmall = 2),
        "x",
        ifelse(sign(cf[1]) == 1, "+", "-"),
        format(abs(cf[1]), nsmall = 2),
        sep = ""
    )

    # R2
    r_sqr <- round(summary(mod_lm)$r.squared, 2)
    r_sqr_str <- paste("R2 =", r_sqr)

    # p-value
    f <- summary(mod_lm)$fstatistic
    p_val <- pf(f[1], f[2], f[3], lower.tail = FALSE)
    if (ifroundp) {
        if (p_val < 0.05) {
            p_val_str <- "p-value < 0.05"
        } else {
            p_val_str <- paste("p-value", round(p_val, 2))
        }
    } else {
        p_val_str <- paste(
            "p-value: ", 
            formatC(p_val, format = "e", digits = 2)
        )
    }

    legend(text_pos,
        text.col = text_col,
        legend = c(eq, r_sqr_str, p_val_str),
        bty = "n", 
        x.intersp = 0
    )
}

