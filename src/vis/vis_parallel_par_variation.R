#'******************************************************************************
#' Description: Visualize Parallel model parameter variation result generated
#' by `mod_hpc_parallel_par_variation.R`.
#'******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")

library(data.table)
library(lubridate)


pa_par_var <- readRDS(file.path(hpc_local_dir, 
    "Pipeline", 
    "parallel_par_variation.Rds"
))

par_true <- pa_par_var$par_true
est_par <- pa_par_var$est_par


doy <- c((-110):(-1), 1:255)

# Convert index back to DOY
est_par[, 2] <- doy[round(est_par[, 2])]
est_par[, 3] <- doy[round(est_par[, 3])]


# # Bootstrapping results
# params <- readRDS(file.path(hpc_local_dir, "Pipeline", "mod_par_cor_1e5.Rds"))
# pa_boot_var <- params$PA


# fig: PA parameter unstability
{
    png(file.path("out", "PA_par_variation_1e5.png"),
        res = 150, width = 1500, height = 600
    )

    # Boxplots
    first_row <- matrix(1:ncol(est_par), nrow = 1)
    second_row <- matrix(rep(1:ncol(est_par) + ncol(est_par), 10), 
        nrow = 10, byrow = TRUE
    )
    layout(rbind(first_row, second_row))
    par(mgp = c(1.5, 0.5, 0), mar = c(0, 0, 2, 1), oma = c(0, 4, 0, 0), 
        bty = "n"
    )
    null <- lapply(seq_along(est_par[1, ]), function(i) {
        boxplot(est_par[, i], 
            xaxt = "n", yaxt = "n", xlab = "", ylab = "",
            horizontal = TRUE,
            width = 0.01, outline = FALSE
        )
    })

    # par(mfrow = c(1, ncol(est_par)), mgp = c(1.5, 0.5, 0))
    par(mar = c(4, 0, 0, 1), cex.lab = 1.5)
    null <- lapply(seq_along(est_par[1, ]), function(i) {
        plot(NA,
            xlim = range(est_par[, i], na.rm = TRUE),
            ylim = range(seq_along(est_par[, i])),
            bty = "n", yaxt = "n", ylab = "", xlab = "",
            xpd = NA,
            cex.axis = 1.2
        )
        mtext(
            side = 1,
            text = switch(colnames(est_par)[i],
                "RMSE" = "RMSE",
                "t0" = expression(t[0]),
                "t0chl" = expression(t[0[chil]]),
                "T_base" = expression(T[base]),
                "T_opt" = expression(T[opt]),
                "T_min" = expression(T[min]),
                "T_max" = expression(T[max]),
                "C_min" = expression(C[min]),
                "F_crit" = expression(F[crit]),
                "C_req" = expression(C[req])
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
            type = "b",
            pch = 16
        )
        # Mark the true parameter value
        points(est_par[1, i], 1, pch = 16, col = "red", cex = 1.5)
    })

    dev.off()
}
