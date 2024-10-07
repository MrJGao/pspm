library(R6)

Visualization <- R6Class("Visualization",
    private = list(
        
    ),

    public = list(

        # Color palatte from Jenna
        jenna_pal = c(
            "#383aaa", 
            "#3A86FF", 
            "#8338EC", 
            "#FF006E", 
            "#FB5607", 
            "#FFBE0B"
        ), 

        #' Compare 2 vectors with a linear equation legend
        CompareXY_J = function(x, y = NULL, ...) {
            if (is.null(y)) {
                tmp <- x[, 1]
                y <- unlist(x[, 2])
                x <- unlist(tmp)
            }
            range <- range(x, y, na.rm = TRUE)
            plot(x, y, 
                pch = 16, xlim = range, ylim = range, 
                mgp = c(1.5, 0.5, 0), 
                ...
            )
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
                legend("bottomright", legend = r_sqr_str, 
                    bty = "n", y.intersp = 5
                )
                legend("bottomright", legend = rmse_str, 
                    bty = "n", y.intersp = 3
                )
                legend("bottomright", legend = p_val_str, bty = "n")
            }
        },

        #' Plot stretched raster image with a nice color and legend
        #' @example PlotStretchLegend(r, viridis(200))
        PlotStretchLegend = function(r, colorRamp, 
            digits = 0, labSize = 1, ...
        ) {
            pal <- colorRampPalette(colorRamp)
            qs <- quantile(r, c(0, 0.02, 0.98, 1))
            if (qs[1] == qs[2]) {
                qs[1] <- qs[1] - 1e-6
            }
            if (qs[3] == qs[4]) {
                qs[4] <- qs[4] + 1e-6
            }
            
            r_breaks <- c(qs[1], seq(qs[2], qs[3], len = 255), qs[4])
            plot(r, col = pal(length(r_breaks) - 1), 
                breaks = r_breaks, axes = F, box = F, legend = F, 
                ...
            )
            # add a reasonable legend
            legend_at <- round(seq(r_breaks[2], r_breaks[length(r_breaks) - 1], 
                len = 7), digits
            )
            legend_labels <- c(
                paste("<", legend_at[1]),
                as.character(legend_at[2:(length(legend_at) - 1)]),
                paste(">", legend_at[length(legend_at)])
            )
            plot(raster(matrix(legend_at)),
                legend.only = T, col = pal(length(r_breaks) - 1),
                axis.args = list(at = legend_at, labels = legend_labels, 
                    cex.axis = labSize
                )
            )
        },

        #' Plot the density and the 95% / 50% CI of a variable
        PlotDensity = function(x) {
            plot(density(x))
            abline(v = quantile(x, c(0.025, 0.975)), lty = 1)
            abline(v = quantile(x, c(0.25, 0.75)), lty = 2)
        },

        # Dark theme par setting
        DarkPar = function() {
            par(bg = NA, fg = "white", col.axis = "white", col.lab = "white",
                col.main = "white", col.sub = "white")
        },

        #' Do smoothScatter with a linear equation legend
        #' Usage is the same with basic smoothScatter.
        #' The background color can be specified to fit black background.
        smoothScatter_J = function(x, y, nbin = 300, nrpoints = 0, 
            bgCol = rgb(1, 1, 1, 0), ifroundp = FALSE, ...
        ) {
            require(RColorBrewer)
            pal <- colorRampPalette(
                c(bgCol, "#f9f9fd", "#a3a6e1", rev(hcl.colors(11, "Spectral"))), 
                alpha = TRUE
            )
            smoothScatter(x, y, colramp = pal, nbin = nbin, 
                nrpoints = nrpoints,
                ...
            )
            fit <- lm(y ~ x)
            abline(fit, col = "blue")

            # draw a 1:1 line
            abline(a = 0, b = 1, lty = 2)

            cf <- round(coef(fit), 2)
            eq <- paste("y=", ifelse(sign(cf[2]) == 1, "", "-"), abs(cf[2]), 
                "x",
                ifelse(sign(cf[1]) == 1, "+", "-"), abs(cf[1]), ", ",
                sep = ""
            )
            r_sqr <- round(summary(fit)$r.squared, 2)
            equation <- bquote(.(eq) ~ R^2 == .(r_sqr))

            f <- summary(fit)$fstatistic
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

            legend("bottomright", legend = equation, bty = "n", y.intersp = 3)
            legend("bottomright", legend = c("", p_val_str), bty = "n")
        },


        #' Create a residual plot for phenology models.
        #' @param x Model fitted values.
        #' @param resid Model residuals.
        ResidualPlot = function(x, resid, ...) {
            plot(x, resid, xlab = "Fitted", ylab = "Residuals", ...)
            lo <- loess(resid ~ x)
            newx <- seq(range(x)[1], range(x)[2], by = 1)
            lines(newx, predict(lo, newx), col = "red")
            abline(h = 0, lty = 2)
        },

        CalSiteAnomaly = function(model_fit, data_li) {
            anomaly_dt <- data.table()
            for (site in unique(data_li$site)) {
                # site <- unique(data_li$site)[1]
                site_dt <- data_li$transition_dates[data_li$site == site]
                site_anomaly <- scale(site_dt, scale = FALSE)
                site_TT_pred <- model_fit$TT$fitted_doy[data_li$site == site]
                site_PA_pred <- model_fit$PA$fitted_doy[data_li$site == site]
                site_SQ_pred <- model_fit$SQ$fitted_doy[data_li$site == site]
                site_AT_pred <- model_fit$AT$fitted_doy[data_li$site == site]
                site_UN_pred <- model_fit$UN$fitted_doy[data_li$site == site]
                site_TT_anomaly <- scale(site_TT_pred, scale = FALSE)
                site_PA_anomaly <- scale(site_PA_pred, scale = FALSE)
                site_SQ_anomaly <- scale(site_SQ_pred, scale = FALSE)
                site_AT_anomaly <- scale(site_AT_pred, scale = FALSE)
                site_UN_anomaly <- scale(site_UN_pred, scale = FALSE)

                anomaly_dt <- rbind(anomaly_dt, data.table(
                    site = site, obs_anom = site_anomaly,
                    TT_anom = site_TT_anomaly, PA_anom = site_PA_anomaly,
                    SQ_anom = site_SQ_anomaly, AT_anom = site_AT_anomaly,
                    UN_anom = site_UN_anomaly
                ))
            }

            return(anomaly_dt)
        },


        ScatterXY = function(x, y, ptcol = "grey", range = c(-30, 30), ...) {
            
            plot(x, y, pch = 21, xlim = range, ylim = range, 
                mgp = c(1.5, 0.5, 0), 
                bg = Transparent(ptcol, 1), lwd = 0.02, tck = 0.02, 
                cex.axis = 1.1,
                ...
            )
            # draw a 1:1 line
            abline(a = 0, b = 1, lty = 2)

            fit <- lm(y ~ x)
            if (any(is.na(coef(fit))) == FALSE) {
                abline(fit, lwd = 2)
                # get equation
                cf <- round(coef(fit), 2)
                eq <- paste("y=", ifelse(sign(cf[2]) == 1, "", "-"), 
                    abs(cf[2]), "x",
                    ifelse(sign(cf[1]) == 1, "+", "-"), abs(cf[1]),
                    sep = ""
                )
                # equation <- bquote(.(eq) ~ R^2 == .(r_sqr))

                # get stats
                r_sqr <- round(summary(fit)$r.squared, 2)
                r_sqr_str <- bquote(R^2 == .(r_sqr))

                f <- summary(fit)$fstatistic
                p_val <- pf(f[1], f[2], f[3], lower.tail = FALSE)
                p_val_str <- paste("p-value: ", 
                    formatC(p_val, format = "e", digits = 2)
                )

                rmse <- round(sqrt(mean((y - x)^2, na.rm = TRUE)), 2)
                rmse_str <- bquote(RMSE == .(rmse))

                if (is.na(r_sqr)) {
                    rmse_str <- "RMSE=NaN"
                }

                # legend("bottomright", legend = eq, bty = "n", y.intersp = 7)
                legend("bottomright", legend = r_sqr_str, bty = "n", 
                    cex = 1.2, y.intersp = 3)
                legend("bottomright", legend = rmse_str, bty = "n", 
                    cex = 1.2, y.intersp = 1)
            }
        },

        ScatterType2XY = function(x, y, ptcol = "grey", range = c(-30, 30), 
            ...
        ) {
            require(lmodel2)
            plot(x, y,
                pch = 21, xlim = range, ylim = range, mgp = c(1.5, 0.5, 0),
                bg = Transparent(ptcol, 1), lwd = 0.02, tck = 0.02, 
                cex.axis = 1.1,
                ...
            )
            # draw a 1:1 line
            abline(a = 0, b = 1, lty = 2)

            fit <- tryCatch({
                suppressMessages(
                    lmodel2(y ~ x, data.frame(x = x, y = y), 
                        "interval", "interval"
                    )
                )
            }, error = function(x) return(NULL))
            
            if (!is.null(fit) & any(is.na(coef(fit))) == FALSE) {
                fit_res <- data.table(fit$regression.results)
                inter_slp <- fit_res[Method == "RMA", .(Intercept, Slope)]
                abline(a = inter_slp$Intercept, b = inter_slp$Slope, lwd = 2)

                # get stats
                r_sqr <- round(fit$rsquare, 2)
                r_sqr_str <- bquote(R^2 == .(r_sqr))

                rmse <- round(sqrt(mean((y - x)^2, na.rm = TRUE)), 2)
                rmse_str <- bquote(RMSE == .(rmse))

                if (is.na(r_sqr)) {
                    rmse_str <- "RMSE=NaN"
                }

                legend("bottomright",
                    legend = r_sqr_str, bty = "n",
                    cex = 1.2, y.intersp = 3
                )
                legend("bottomright",
                    legend = rmse_str, bty = "n",
                    cex = 1.2, y.intersp = 1
                )
            }
        },


        Reformat2Dt = function(cv_li) {
            dt <- lapply(cv_li, function(fold) {
                dt <- data.table(
                    site = fold$site,
                    obs = fold$obs,
                    TT = fold$TT, PA = fold$PA, SQ = fold$SQ,
                    AT = fold$AT, UN = fold$UN
                )
            })
            dt <- do.call(rbind, dt)
            names(dt) <- c("site", "obs", "TT", "PA", "SQ", "AT", "UN")

            dt[dt == 9999] <- NA
            return(dt)
        }

    ),

    active = list(
        
    )
)
