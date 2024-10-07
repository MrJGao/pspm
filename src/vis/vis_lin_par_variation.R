#'******************************************************************************
#' Description: Show the parameter variation in the LIN model, using the same
#' procedure as the Parrallel model figure in `vis_hpc_simulate_pheno.R`.
#'******************************************************************************
rm(list = ls())

library(data.table)
library(lubridate)

source("src/base.R")
source("src/vis/vis_base.R")
source("src/hlp/hlp_fit_spm.R")
source("src/hlp/hlp_spring_pheno_models.R")
source("src/vis/vis_spring_pheno_model_diag_figs.R")
source("src/mod_compare_base.R")




# ~ Read in some temperature data observed from the HB ground ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hb_dt <- fread(file.path(gdir, 
    "Data/ARD", 
    "hb_grd_temperature_pheno.csv"
))

# This site was selected b/c it has 27 site years
site_dt <- hb_dt[siteID == "1B"]

site_li <- FormatDataForPhenoModel(site_dt)



# ~ Simulation ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assume that the LIN model is the true model that works in nature.

# Parameter values
t1 <- 80 # Preseason start day
t2 <- 130 # Preseason end day

par_true <- c(t1, t2)
# Simulate SOS data
sim_sos <- LIN(par = par_true, data = site_li)

# Now, fit the similated sos data by the models
sim_li <- site_li
sim_li$transition_dates <- sim_sos


set.seed(8964)
est_par <- lapply(1:20, function(i) {
    # Simulate data
    sim_li <- site_li
    # Slightly change temperatures to add some randomness
    sim_li$Ti <- site_li$Ti + rnorm(length(site_li$Ti), 0, 0.01)
    # Simulate SOS
    sim_sos <- LIN(par = par_true, data = sim_li)
    sim_li$transition_dates <- sim_sos

    # Fit the model using the simulated data
    mod_LIN <- FitTheModel(
        model = LIN,
        cost_fun = CostRMSE,
        data = sim_li
    )
    pars <- mod_LIN$optim$par
    rmse <- mod_LIN$optim$value
    return(c(rmse, pars))
})
est_par <- do.call(rbind, est_par)
est_par <- rbind(c(0, par_true), est_par)
colnames(est_par) <- c("RMSE", "t1", "t2")



# fig: LIN parameter stability
{
    png(file.path("Output/Assets", "LIN_par_variation.png"),
        res = 150, width = 1000, height = 500
    )

    par(mfrow = c(1, ncol(est_par)), mgp = c(1.5, 0.5, 0))
    par(mar = c(4, 0, 2, 1), oma = c(0, 4, 0, 0), cex.lab = 1.2)
    null <- lapply(seq_along(est_par[1, ]), function(i) {
        xrg <- switch(colnames(est_par)[i],
            "RMSE" = c(-1, 1),
            "t1" = c(par_true[1] - 5, par_true[1] + 5),
            "t2" = c(par_true[2] - 5, par_true[2] + 5),
        )
        plot(NA,
            xlim = xrg,
            ylim = range(seq_along(est_par[, i])),
            bty = "n", yaxt = "n", ylab = "", xlab = "",
            xpd = NA
        )
        mtext(
            side = 1,
            text = switch(colnames(est_par)[i],
                "RMSE" = "RMSE",
                "t1" = expression(t[1]),
                "t2" = expression(t[2]),
            ),
            line = 2.5
        )
        grid(ny = nrow(est_par) + 1)
        if (i == 1) {
            axis(
                side = 2,
                at = 1:nrow(est_par),
                labels = c("Truth", 1:(nrow(est_par) - 1)),
                las = 1
            )
            mtext(side = 2, text = "Iteration", line = 1.5)
        }
        lines(est_par[, i], seq_along(est_par[, i]),
            type = "o",
            pch = 16
        )
        # Mark the true parameter value
        points(est_par[1, i], 1, pch = 16, col = mycolor[1], cex = 1.5)
    })

    dev.off()
}

