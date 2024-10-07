# ******************************************************************************
# Combine the results and make the figure
# 
# Author: Xiaojie Gao
# Date: 2024-08-18
# ******************************************************************************
library(data.table)
library(ggplot2)



PlotLoess <- function(x, y) {
    mod_loess <- loess(y ~ x)
    newx <- seq(grconvertX(0, "npc"), grconvertX(1, "npc"), len = 100)
    lines(newx, predict(mod_loess, newx), col = "blue", lwd = 2)
}

files <- list.files("pipe/fake_chill_forc3", ".Rds$", full.names = TRUE)
res_li <- lapply(files, readRDS)
res_dt <- do.call(rbind, res_li)


png("out/fake_chill_forc3.png", width = 1200, height = 450, res = 150)
par(mfrow = c(1, 3), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 2, 1))

plot(res_dt[, .(tc, sig_pct)], pch = 16, ylab = "Sig.%", xlab = "t_chill")
PlotLoess(res_dt$tc, res_dt$sig_pct)

plot(res_dt[, .(Tb, sig_pct)], pch = 16, ylab = "Sig.%", xlab = "Tb")
PlotLoess(res_dt$Tb, res_dt$sig_pct)

plot(res_dt[, .(t0, sig_pct)], pch = 16, ylab = "Sig.%", xlab = "t0")
PlotLoess(res_dt$t0, res_dt$sig_pct)

dev.off()

