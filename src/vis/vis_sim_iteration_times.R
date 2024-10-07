#'******************************************************************************
#' Description: Visualize the effect of number of maximum iteration times on the 
#' model fitting and parameter estimation process.
#'******************************************************************************
library(data.table)
library(plotrix)

source("src/base.R")

# in:
# True paramters
par_true_file <- file.path(hpc_local_dir, "Pipeline", "PA_sim_par_true.csv")
# Simulation using the same data
sim_iter_same_file <- file.path(hpc_local_dir, 
    "Pipeline", 
    "PA_sim_iter_same_data.csv"
)
# Simulation using the noisy data
sim_iter_noisy_file <- file.path(hpc_local_dir, 
    "Pipeline", 
    "PA_sim_iter_noisy_data.csv"
)


par_true_dt <- fread(par_true_file)
sim_iter_same_dt <- fread(sim_iter_same_file)
sim_iter_noisy_dt <- fread(sim_iter_noisy_file)



# fig: Same data, different seeds, to test when the model converges
{
    png(file.path("Output/Assets", "PA_sim_iteration.png"),
        width = 1200, height = 1000, res = 150
    )

    par(mfcol = c(10, 2), mgp = c(2, 0.5, 0), 
        mar = c(0.5, 4, 0, 1), 
        oma = c(4, 4, 4, 1)
    )
    # RMSE
    plot(sim_iter_same_dt[, .(Maxit, RMSE)], type = "l", 
        xaxt = "n", xlab = "",
        ylim = c(0, 2)
    )
    mtext(side = 3, "Using measured temperature", line = 1)

    # Parameters
    null <- lapply(
        colnames(sim_iter_same_dt)[4:ncol(sim_iter_same_dt)], 
        function(par_name) {
            true_val <- par_true_dt[name == par_name, par_true]
            est_val <- sim_iter_same_dt[, get(par_name)]
            yrange <- range(c(true_val, est_val))
            
            plot(sim_iter_same_dt[, .(Maxit, get(par_name))], 
                ylim = yrange,
                type = "l", 
                xlab = "", xaxt = "n", 
                ylab = ""
            )
            mtext(side = 2, par_name, line = 2, cex = 0.7)
            abline(lty = 2, h = true_val)
        }
    )
    axis(side = 1, at = c(1, 5e4, 1e5, 5e5, 1e6), las = 2)
    mtext(side = 1, "Max number of iteration", line = 3)


    # Noisy data, different seeds, to test when the model converges
    # RMSE
    plot(x = sim_iter_noisy_dt$Maxit, y = sim_iter_noisy_dt$RMSE,
        type = "l",
        ylim = c(0, 2),
        xaxt = "n", xlab = "",
        ylab = "RMSE"
    )
    mtext(side = 3, "Using temperature with random noise", line = 1)

    # Parameters
    null <- lapply(
        colnames(sim_iter_noisy_dt)[4:ncol(sim_iter_noisy_dt)],
        function(par_name) {
            true_val <- par_true_dt[name == par_name, par_true]
            est_val <- sim_iter_noisy_dt[, get(par_name)]
            yrange <- range(c(true_val, est_val))

            plot(sim_iter_noisy_dt[, .(Maxit, get(par_name))],
                ylim = yrange,
                type = "l",
                xlab = "", xaxt = "n",
                ylab = ""
            )
            mtext(side = 2, par_name, line = 2, cex = 0.7)
            abline(lty = 2, h = true_val)
        }
    )
    axis(side = 1, at = c(1, 5e4, 1e5, 5e5, 1e6), las = 2)
    mtext(side = 1, "Max number of iteration", line = 3)


    dev.off()
}
