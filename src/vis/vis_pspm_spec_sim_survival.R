#'******************************************************************************
#' Description: Visualize the result of survival models fitting to simulations.
#'******************************************************************************

dev_diff <- readRDS(file.path(gdir, "Pipeline", "pspm_spec_sim_survival.Rds"))

png(file.path("out", "sim_mix_survival.png"), 
    width = 1000, height = 800, res = 150
)
# png(file.path("out", "sim_mix_survival_dark.png"), 
#     width = 2000, height = 1200, res = 300
# )
# par(bg = NA, fg = "white", col.axis = "white", col.lab = "white", 
#     col.main = "white", col.sub = "white")


par(bty = "L", mar = c(3, 9, 2, 1), mgp = c(1.5, 0.5, 0), tck = -0.02)

boxplot(dif ~ case,
    data = dev_diff,
    ylab = "", xlab = "Deviance difference",
    las = 1,
    horizontal = TRUE,
    names = c(
        "Equally mix", "Chilling mix",
        "70%TT+30%PA", "50%TT+50%PA", "30%TT+70%PA",
        "70%TT+30%SQ", "50%TT+50%SQ", "30%TT+70%SQ",
        "70%TT+30%AT", "50%TT+50%AT", "30%TT+70%AT",
        "70%TT+30%UN", "50%TT+50%UN", "30%TT+70%UN"
    ),
    outline = FALSE,
    tck = 0.01
)
mtext(side = 2, text = "Cases", line = 7)
abline(v = 5, lty = 2)
text(x = 5, y = 14, 
    labels = substitute(italic("*significance")), 
    adj = c(-0.05, -1.5),
    xpd = NA
)


dev.off()