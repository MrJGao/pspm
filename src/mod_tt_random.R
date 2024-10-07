#'******************************************************************************
#' Description: Apply generated random temperature to the Thermal Time Model.
#'******************************************************************************

# bsub < src/mod_hpc_tt_random.csh

library(data.table)
library(lubridate)

source("src/base.R")
source("src/vis/vis_base.R")
source("src/mod_compare_base.R")
source("src/hlp/hlp_model_cost_functions.R")
source("src/hlp/hlp_spring_pheno_models.R")
source("src/hlp/hlp_fit_spm.R")
source("src/vis/vis_spring_pheno_model_diag_figs.R")

library(Rmpi)
library(parallel)
library(snow)



# in:
hb_hf_pheno_dt_file <- file.path("data", "hb_hf_pheno_temp.csv")
# hb_hf_pheno_dt <- fread(file.path(gdir, "Data/ARD", "hb_hf_pheno_temp.csv"))

# out: TT_random_simulation_fit.csv
tt_random_sim_fit_file <- file.path(hpc_dir, 
    "Pipeline", 
    "TT_random_simulation_fit.csv"
)
# out: TT_random_simulation_cv.csv
tt_random_sim_cv_file <- file.path(hpc_dir, 
    "Pipeline", 
    "TT_random_simulation_cv.csv"
)
# out: LIN_random_simulation_fit.csv
lin_random_sim_fit_file <- file.path(hpc_dir,
    "Pipeline",
    "LIN_random_simulation_fit.csv"
)
# out: LIN_random_simulation_cv.csv
lin_random_sim_cv_file <- file.path(hpc_dir,
    "Pipeline",
    "LIN_random_simulation_cv.csv"
)

# ~ Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hb_hf_pheno_dt <- fread(hb_hf_pheno_dt_file)

# Use site "1B" as it contains 27 site years of observations.
site_dt <- hb_hf_pheno_dt[siteID == "1B"]


CompareXY_J2 <- function(x, y, ...) {
    range <- range(x, y)
    # plot(x, y,
    #     pch = 16, xlim = range, ylim = range, col = "grey",
    #     mgp = c(1.5, 0.5, 0), ...
    # )
    # draw a 1:1 line
    # abline(a = 0, b = 1, lty = 2, lwd = 1)

    fit <- lm(y ~ x)
    if (any(is.na(coef(fit))) == FALSE) {
        # abline(fit, col = "white", lwd = 2)
        # abline(fit, col = "black", lwd = 2)

        # Slope must be positive, negative slope means model doesn't work
        if (coef(fit)[2] < 0) {
            return(list(rmse = NA, r_sqr = NA, p_val = NA))
        }

        # get equation
        cf <- round(coef(fit), 2)
        eq <- paste("y=", ifelse(sign(cf[2]) == 1, "", "-"), abs(cf[2]), "x",
            ifelse(sign(cf[1]) == 1, "+", "-"), abs(cf[1]),
            sep = ""
        )
        # equation <- bquote(.(eq) ~ R^2 == .(r_sqr))

        # get stats
        r_sqr <- round(summary(fit)$r.squared, 2)
        r_sqr_str <- bquote(R^2 == .(r_sqr))

        f <- summary(fit)$fstatistic
        p_val <- pf(f[1], f[2], f[3], lower.tail = FALSE)
        p_val_str <- paste("p-value: ", formatC(p_val, format = "e", digits = 2))

        rmse <- round(sqrt(mean((y - x)^2)), 2)
        rmse_str <- bquote(RMSE == .(rmse))

        # legend("bottomright", legend = eq, bty = "n", y.intersp = 7)
        # legend("bottomright", legend = r_sqr_str, bty = "n", y.intersp = 5)
        # legend("bottomright", legend = rmse_str, bty = "n", y.intersp = 3)
        # legend("bottomright", legend = p_val_str, bty = "n")

        return(list(rmse = rmse, r_sqr = r_sqr, p_val = p_val))
    }
}


WhatAreWeModeling <- function(dt) {
    samp_li <- FormatDataForPhenoModel(dt[, .(
        siteID = siteID, Tmax = Tmean,
        Tmin = Tmean, Tmean = Tmean_new, PhenoYear, SOS
    )])


    mod_TT <- FitTheModel(
        model = ThermalTimeModel,
        cost_fun = CostRMSE,
        data = samp_li,
        control = list(maxit = 5000)
    )

    res <- CompareXY_J2(
        mod_TT$fitted_doy,
        samp_li$transition_dates,
        xlab = "Est.", ylab = "Obs."
    )

    return(res)
}


# Leave one year out cv
CvRandom <- function(dt) {
    samp_li <- FormatDataForPhenoModel(dt[, .(
        siteID = siteID, Tmax = Tmean,
        Tmin = Tmean, Tmean = Tmean_new, PhenoYear, SOS
    )])

    # Iteratively sample one site-year out, then cross validation
    pred <- lapply(samp_li$year, function(yr) {
        test_idx <- which(samp_li$year == yr)
        test_li <- GetDataLiIndex(samp_li, test_idx)

        train_idx <- which(samp_li$year != yr)
        trian_li <- GetDataLiIndex(samp_li, train_idx)

        # Train
        mod_TT <- FitTheModel(
            model = ThermalTimeModel,
            cost_fun = CostRMSE,
            data = trian_li,
            control = list(maxit = 3000)
        )
        # Test
        pred <- do.call(
            ThermalTimeModel,
            list(par = mod_TT$optim$par, data = test_li)
        )
        return(pred)
    })
    pred <- do.call(c, pred)

    res <- CompareXY_J2(
        pred,
        samp_li$transition_dates,
        xlab = "Est.", ylab = "Obs."
    )

    return(res)
}


# Linear regression
LinFit <- function(dt) {
    month_mean <- unique(
        dt[month(Date) == 5, 
            .(Tmean = mean(Tmean_new), SOS = SOS),
            by = .(siteID, PhenoYear)
        ]
    )

    lm_fit <- lm(SOS ~ Tmean, data = month_mean)

    res <- CompareXY_J2(
        fitted(lm_fit),
        month_mean$SOS,
        xlab = "Est.", ylab = "Obs."
    )

    return(res)
}

# Linear regression cross-validation
LinCV <- function(dt) {
    month_mean <- unique(
        dt[month(Date) == 5,
            .(Tmean = mean(Tmean_new), SOS = SOS),
            by = .(siteID, PhenoYear)
        ]
    )

    pred <- lapply(month_mean$PhenoYear, function(yr) {
        test_dt <- month_mean[PhenoYear == yr]
        train_dt <- month_mean[PhenoYear != yr]

        lm_fit <- lm(SOS ~ Tmean, data = train_dt)
        pred <- predict(lm_fit, newdata = test_dt)

        return(pred)
    })
    pred <- do.call(c, pred)

    res <- CompareXY_J2(
        pred,
        month_mean$SOS,
        xlab = "Est.", ylab = "Obs."
    )

    return(res)
}


# ~ Cluster config ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make cluster
cl <- makeCluster(mpi.universe.size() - 1, type = "MPI")
calls <- clusterCall(cl, function() {
    suppressWarnings({
        library(data.table)
        library(lubridate)

        source("src/base.R")
        source("src/vis/vis_base.R")
        source("src/mod_compare_base.R")
        source("src/hlp/hlp_model_cost_functions.R")
        source("src/hlp/hlp_spring_pheno_models.R")
        source("src/hlp/hlp_fit_spm.R")
        source("src/vis/vis_spring_pheno_model_diag_figs.R")
    })
})
clusterExport(cl, c(
    "site_dt", "CompareXY_J2", "WhatAreWeModeling", "CvRandom", 
    "LinFit", "LinCV"
))



# ~ Goodness-of-fit ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(8964)
fit_dt <- clusterApply(cl, x = 1:1000, function(x) {
    site_dt$Tmean_new <- rnorm(length(site_dt$Tmean), 4, 10)
    res <- WhatAreWeModeling(site_dt)
    return(as.data.frame(res))
})
fit_dt <- do.call(rbind, fit_dt)

# Export
fwrite(fit_dt, file = tt_random_sim_fit_file)


lin_fit_dt <- clusterApply(cl, x = 1:1000, function(x) {
    site_dt$Tmean_new <- rnorm(length(site_dt$Tmean), 4, 10)
    res <- LinFit(site_dt)
    return(as.data.frame(res))
})
lin_fit_dt <- do.call(rbind, lin_fit_dt)

# Export
fwrite(lin_fit_dt, file = lin_random_sim_fit_file)



# ~ Leave-one-year-out-cv ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(8964)
cv_dt <- clusterApply(cl, x = 1:1000, function(x) {
    site_dt$Tmean_new <- rnorm(length(site_dt$Tmean), 4, 10)
    res <- CvRandom(site_dt)
    return(as.data.frame(res))
})
cv_dt <- do.call(rbind, cv_dt)

# Export
fwrite(cv_dt, file = tt_random_sim_cv_file)


lin_cv_dt <- clusterApply(cl, x = 1:1000, function(x) {
    site_dt$Tmean_new <- rnorm(length(site_dt$Tmean), 4, 10)
    res <- LinCV(site_dt)
    return(as.data.frame(res))
})
lin_cv_dt <- do.call(rbind, lin_cv_dt)

# Export
fwrite(lin_cv_dt, file = lin_random_sim_cv_file)




# shut down cluster
snow::stopCluster(cl)
Rmpi::mpi.quit()