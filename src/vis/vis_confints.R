#'******************************************************************************
#' Description: Visualize confidence intervals for spring phenology models.
#'******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")


library(vioplot)
library(data.table)


# Example data
hb_hf_li <- readRDS("pipe/hb_hf_li.Rds")
# Bootstrapping results
params <- readRDS("pipe/mod_par_cor_1e5.Rds")

theme_color <- mcd_dm_color

par_TT <- params$TT
par_PA <- params$PA
par_SQ <- params$SQ
par_AT <- params$AT
par_UN <- params$UN

DOYs <- hb_hf_li$doy

# param = par_TT$T_base
# var_name = "T[base](~degree~C)"
PlotConfint <- function(param, var_name) {
    conf_int <- quantile(param, c(0.025, 0.975))
    var <- param[between(param, lower = conf_int[1], upper = conf_int[2])]

    # vioplot(var, las = 1, xaxt = "n", col = "grey90")
    vioplot(var, las = 1, xaxt = "n", col = theme_color)
    title(str2expression(var_name), line = 1)
}


png(file.path("out", "par_confint_1e5.png"), 
    width = 2000, height = 1800, res = 200)

# png(file.path("out", "par_confint_1e5_dark.png"), 
#     width = 2000, height = 1800, res = 200)
# par(bg = NA, fg = "white", col.axis = "white", col.lab = "white", 
#     col.main = "white", col.sub = "white")

par(mgp = c(1.5, 0.5, 0), mar = c(1, 2, 2, 1.5), oma = c(0, 3, 0, 0), 
    cex.main = 1.5, cex.axis = 1.2)

layout(matrix(c(
    1,  0,  2,  0,  0,  0,  0,  0, 3,  # TT
    4,  5,  6,  7,  8,  9, 10, 11, 12, # PA
    13, 14, 15, 16, 17, 18, 0, 19, 20, # SQ
    21,  0, 22, 23, 24, 25, 0,  0,  0, # AT
    26, 27, 28, 29, 30, 31, 32, 33, 34 # UN
), nrow = 5, byrow = TRUE))

# TT
PlotConfint(DOYs[par_TT$t0], "t[0](DOY)")
PlotConfint(par_TT$T_base, "T[base](~degree~C)")
PlotConfint(par_TT$F_crit, "F[crit]")

# PA
PlotConfint(DOYs[par_PA$t0], "t[0](DOY)")
PlotConfint(DOYs[par_PA$t0_chil], "t[0-chill](DOY)")
PlotConfint(par_PA$T_base, "T[base](~degree~C)")
PlotConfint(par_PA$T_opt, "T[opt](~degree~C)")
PlotConfint(par_PA$T_min, "T[min](~degree~C)")
PlotConfint(par_PA$T_max, "T[max](~degree~C)")
PlotConfint(par_PA$C_min, "C[min]")
PlotConfint(par_PA$C_req, "C[req]")
PlotConfint(par_PA$F_crit, "F[crit]")

# SQ
PlotConfint(DOYs[par_SQ$t0], "t[0](DOY)")
PlotConfint(DOYs[par_SQ$t0_chil], "t[0-chill](DOY)")
PlotConfint(par_SQ$T_base, "T[base](~degree~C)")
PlotConfint(par_SQ$T_opt, "T[opt](~degree~C)")
PlotConfint(par_SQ$T_min, "T[min](~degree~C)")
PlotConfint(par_SQ$T_max, "T[max](~degree~C)")
PlotConfint(par_SQ$C_req, "C[req]")
PlotConfint(par_SQ$F_crit, "F[crit]")

# AT
PlotConfint(DOYs[par_AT$t0], "t[0](DOY)")
PlotConfint(par_AT$T_base, "T[base](~degree~C)")
PlotConfint(par_AT$a, "a")
PlotConfint(par_AT$b, "b")
PlotConfint(par_AT$c, "c")

# UN
PlotConfint(DOYs[par_UN$tc], "t[c](DOY)")
PlotConfint(par_UN$a_c, "a[c]")
PlotConfint(par_UN$b_c, "b[c]")
PlotConfint(par_UN$c_c, "c[c]")
PlotConfint(par_UN$b_f, "b[f]")
PlotConfint(par_UN$c_f, "c[f]")
PlotConfint(par_UN$w, "w")
PlotConfint(par_UN$C_req, "C[req]")
PlotConfint(par_UN$k, "k")


# Annotations
row_height <- 1 / 5
row_ctr <- (0:4) * row_height + row_height / 2
row_bottom <- (0:4) * row_height
row_top <- (1:5) * row_height

rect(xleft = grconvertX(rep(0.005, 5), "ndc"), 
    ybottom = grconvertY(row_bottom + 0.01, "ndc"),
    xright = grconvertX(rep(0.035, 5), "ndc"),
    ytop = grconvertY(row_top - 0.02, "ndc"), 
    xpd = NA, col = theme_color, border = 0)
text(x = grconvertX(0.02, "ndc"), y = grconvertY(row_ctr, "ndc"),
    labels = rev(c("TT", "PA", "SQ", "AT", "UN")), xpd = NA, cex = 1.5)

segments(x0 = grconvertX(rep(0, 4), "ndc"), 
    y0 = grconvertY(row_bottom[2:5], "ndc"), 
    x1 = grconvertX(rep(1, 4), "ndc"), 
    xpd = NA, lwd = 1, lty = 2
)


dev.off()

