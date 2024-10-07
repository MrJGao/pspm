#'******************************************************************************
#' Description: Visualize TT random simulation results.
#'******************************************************************************
library(data.table)
library(lubridate)

source("src/base.R")
source("src/vis/vis_base.R")



# ~ Load-data ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Goodness-of-fit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tt_sim_fit <- fread(file.path(
    hpc_local_dir,
    "Pipeline",
    "TT_random_simulation_fit.csv"
))

lin_sim_fit <- fread(file.path(
    hpc_local_dir,
    "Pipeline",
    "LIN_random_simulation_fit.csv"
))

# Remove NA rows
tt_sim_fit <- tt_sim_fit[!is.na(p_val)]
lin_sim_fit <- lin_sim_fit[!is.na(p_val)]


# ~ CV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tt_sim_cv <- fread(file.path(
    hpc_local_dir,
    "Pipeline",
    "TT_random_simulation_cv.csv"
))

lin_sim_cv <- fread(file.path(
    hpc_local_dir,
    "Pipeline",
    "LIN_random_simulation_cv.csv"
))

# Remove NA rows
tt_sim_cv <- tt_sim_cv[!is.na(p_val), ]
lin_sim_cv <- lin_sim_cv[!is.na(p_val), ]



# ~ Fig ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png(file.path("Output/Assets", "TT_random.png"),
    width = 1200, height = 900, 
    res = 150
)

par(mfcol = c(3, 2), mgp = c(1.5, 1, 0), mar = c(3, 3, 2, 1), cex.axis = 1.2)
# ~ Goodness-of-fit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RMSE
boxplot(tt_sim_fit$rmse, lin_sim_fit$rmse, 
    horizontal = TRUE, axst = "n", 
    main = "RMSE"
)
axis(side = 2, at = 1:2, labels = c("TT", "LIN"), las = 1)

# R^2
boxplot(tt_sim_fit$r_sqr, lin_sim_fit$r_sqr,
    horizontal = TRUE, axst = "n", 
    main = expression(R^2)
)
axis(side = 2, at = 1:2, labels = c("TT", "LIN"), las = 1)

# P-value
boxplot(tt_sim_fit$p_val, lin_sim_fit$p_val,
    horizontal = TRUE, axst = "n", 
    main = "p-value"
)
axis(side = 2, at = 1:2, labels = c("TT", "LIN"), las = 1)
abline(v = 0.05, lty = 2)
text(0.05, 2.5, xpd = NA, "0.05", pos = 4)
tt_pct <- round(nrow(tt_sim_fit[p_val < 0.05]) / nrow(tt_sim_fit) * 100, 1)
lin_pct <- round(nrow(lin_sim_fit[p_val < 0.05]) / nrow(lin_sim_fit) * 100, 1)
text(0.25, 0.8, paste0("% < 0.05: ", tt_pct), pos = 1)
text(0.4, 2, paste0("% < 0.05: ", lin_pct), pos = 1)


# ~ CV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RMSE
boxplot(tt_sim_cv$rmse, lin_sim_cv$rmse,
    horizontal = TRUE, axst = "n", 
    main = "RMSE"
)
axis(side = 2, at = 1:2, labels = c("TT", "LIN"), las = 1)

boxplot(tt_sim_cv$r_sqr, lin_sim_cv$r_sqr, 
    horizontal = TRUE, axst = "n", 
    main = expression(R^2)
)
axis(side = 2, at = 1:2, labels = c("TT", "LIN"), las = 1)

boxplot(tt_sim_cv$p_val, lin_sim_cv$p_val,
    horizontal = TRUE, axst = "n", 
    main = "p-value"
)
axis(side = 2, at = 1:2, labels = c("TT", "LIN"), las = 1)
abline(v = 0.05, lty = 2)
text(0.05, 2.5, xpd = NA, "0.05", pos = 4)
tt_pct <- round(nrow(tt_sim_cv[p_val < 0.05]) / nrow(tt_sim_cv) * 100, 1)
lin_pct <- round(nrow(lin_sim_cv[p_val < 0.05]) / nrow(lin_sim_cv) * 100, 1)
text(0.32, 1, paste0("% < 0.05: ", tt_pct), pos = 1)
text(0.4, 2, paste0("% < 0.05: ", lin_pct), pos = 1)

dev.off()