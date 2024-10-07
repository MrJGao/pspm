#'******************************************************************************
#' Description: Show that T_base and F_crit in the process-based spring models 
#' (including the thermal time model, the parallel model, and the sequential 
#' model) are highly correlated.
#'******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")

# Bootstrapping result
params <- readRDS(file.path("pipe", "mod_par_cor_1e5.Rds"))

theme_color <- mcd_dm_color

PlotCorr <- function(par_dt, title, dark = FALSE) {
    if (dark == FALSE) {
        coloramp <- colorRampPalette(
            # c("white", RColorBrewer::brewer.pal(8, "OrRd"))
            c(
                "#ffffff", "#eff5f5", "#deebeb",
                "#cee1e1", "#bed7d8", "#adcece",
                "#9dc4c4", "#8dbabb", "#7cb1b2",
                "#6ba7a8", "#599d9f", "#469496"
            )
        )
    } else {
        coloramp <- colorRampPalette(viridis::magma(256))
    }
    smoothScatter(par_dt$T_base, par_dt$F_crit,
        xlab = bquote(T["base"]),
        ylab = bquote(F["crit"]), pch = 16,
        colramp = coloramp,
        nrpoints = 0,
        cex = 0.2, cex.lab = 1.5, las = 1
    )
    corr <- round(cor(par_dt$T_base, par_dt$F_crit), 2)
    text(grconvertX(0.65, "npc"), grconvertY(0.7, "npc"),
        label = paste("Correlation =", corr))
    title(title, line = 0.5)
}

# fig: corr T_base vs F_crit and w vs k in the UN model
png("out/corr_tf_wk.png", 
    res = 300, width = 2000, height = 800
)
ifdark <- FALSE

# png("Output/Assets/corr_Tbase_Fcrit_dark.png", 
#     res = 300, width = 1450, height = 4000
# )
# par(bg = NA, fg = "white", col.axis = "white", col.lab = "white",
#     col.main = "white", col.sub = "white")
# ifdark <- TRUE
# par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))


par(mgp = c(2, 0.5, 0), cex = 1.5)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))


# TT
PlotCorr(params$TT, "TT", dark = ifdark)
# PA
PlotCorr(params$PA, "PA", dark = ifdark)
# SQ
PlotCorr(params$SQ, "SQ", dark = ifdark)

dev.off()


