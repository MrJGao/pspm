# ******************************************************************************
# Based on AICc's equation, plot the \Delta AICc space constructed by variations
# of RMSE and number of model parameters. The figure will show that how much
# RMSE ratio we should get to support the selection of a more complex model.
# ******************************************************************************



rm(list=ls())


rmse_ratio <- seq(0.1, 1.15, len = 100)
# RSS ratio between true model and alternative model
rss_ratio <- rmse_ratio^2

# Number of observations
nobs <- 12:length(rss_ratio)


calDeltaAICc <- function(k1, k2 = 3) {
    # Number of parameters for the alternative model
    # k1 <- 9 # for PA and UN model
    # k1 <- 8 # for SQ model
    # k1 <- 5 # for AT model
    # # Number of parameters for the true model
    # k2 <- 3

    # Delta AICc matrix
    aicc_mat <- matrix(NA, nrow = length(rss_ratio), ncol = length(nobs))
    for (i in seq_along(rss_ratio)) {
        i_rss <- rss_ratio[i]
        for (j in seq_along(nobs)) {
            j_n <- nobs[j]

            delta_aicc <- 2 * (k1 - k2) + j_n * log(i_rss) +
                (2 * k1 * (k1 + 1)) / (j_n - k1 - 1) -
                (2 * k2 * (k2 + 1)) / (j_n - k2 - 1)

            aicc_mat[i, j] <- delta_aicc
        }
    }
    return(aicc_mat)
}


# levelplot(aicc_mat, 
#     col.regions = hcl.colors(256, "rdbu"), 
#     xlab = "RMSE ratio", ylab = "Number of observations",
#     row.values = rmse_ratio, column.values = nobs,
#     aspect = "fill",
#     at = at_vec,
#     panel = function(..., at, region) {
#         panel.levelplot(..., at = at, region = TRUE)
#         panel.contourplot(..., at = -2, region = FALSE)
#     }
# )


# {
#     png("Output/AICc_sig_contour.png", width = 2500, height = 1000, res = 300)

#     par(mfrow = c(1, 3), oma = c(3, 0, 0, 0), mar = c(5, 4, 3, 2))
#     k1_vec <- c(9, 8, 5)
#     for (i in 1:length(k1_vec)) {
#         aicc_mat <- calDeltaAICc(k1_vec[i], k2 = 3)
#         qt <- range(aicc_mat)
#         at_vec <- c(
#             seq(qt[1], -1e-5, len = 127),
#             0,
#             seq(1e-5, qt[2], len = 127)
#         )

#         # Plot the image
#         image(
#             x = rmse_ratio, y = nobs, aicc_mat, breaks = at_vec,
#             col = hcl.colors(254, "rdbu"),
#             xlab = "RMSE ratio", ylab = "Number of observations",
#             mgp = c(1.5, 0.5, 0),
#             cex.lab = 1.3
#         )
#         contour(
#             x = rmse_ratio, y = nobs, aicc_mat, levels = -2, labcex = 1,
#             drawlabels = FALSE,
#             add = TRUE
#         )
#     }

#     fields::set.panel()
#     par(mar = c(0, 0, 0, 0), cex = 0.7)
#     fields::imagePlot(
#         aicc_mat,
#         breaks = at_vec,
#         col = hcl.colors(254, "rdbu"),
#         legend.cex = 0.1,
#         legend.line = 1,
#         legend.only = TRUE,
#         horizontal = TRUE,
#         smallplot = c(0.1, 0.9, 0, 0.03)
#     )
#     text(grconvertX(c(0.02, 0.34, 0.67), "ndc"), grconvertY(rep(0.89, 3), "ndc"),
#         labels = letters[1:3], xpd = NA,
#         cex = 1.5
#     )
#     dev.off()
# }

{
    # png("out/AICc_sig_contour.png", width = 1500, height = 1300, res = 300)
    
    png("out/AICc_sig_contour_dark.png", width = 1500, height = 1300, res = 300)
    par(bg = NA, fg = "white", col.axis = "white", col.lab = "white", 
        col.main = "white", col.sub = "white")

    par(mar = c(3, 3, 1, 1), cex = 1.2)
    plot(NA,
        xlim = range(rmse_ratio), ylim = range(nobs), mgp = c(1.5, 0.5, 0),
        xlab = "RMSE ratio", ylab = "Number of observations"
    )
    k1_vec <- c(9, 8, 5)
    for (i in 1:length(k1_vec)) {
        aicc_mat <- calDeltaAICc(k1_vec[i], k2 = 3)
        qt <- range(aicc_mat)
        at_vec <- c(
            seq(qt[1], -1e-5, len = 127),
            0,
            seq(1e-5, qt[2], len = 127)
        )

        contour(
            x = rmse_ratio, y = nobs, aicc_mat, levels = -2, labcex = 1,
            drawlabels = FALSE,
            add = TRUE, lty = i
        )
    }
    text(0.6, 50, label = expression(~Delta ~ AICc < -2))
    text(1, 20, label = expression(~Delta ~ AICc > -2))
    legend("topleft",
        lty = seq_along(k1_vec), bty = "n",
        legend = c("k=9", "k=8", "k=5")
    )

    dev.off()
}
