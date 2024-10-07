rm(list = ls())

library(data.table)
library(lubridate)

source("src/base.R")
source("src/vis/vis_base.R")
source("src/hlp/hlp_model_cost_functions.R")
source("src/hlp/hlp_spring_pheno_models.R")
source("src/hlp/hlp_fit_spm.R")
source("src/vis/vis_spring_pheno_model_diag_figs.R")





Forcing <- function(t0) {
    # rf <- ifelse(temp$Tmean - t_base > 0, 1, 0)
    rf <- temp$Tmean - t_base
    rf[rf < 0] <- 0
    rf[1:t0] <- 0
    rf <- cumsum(rf)
    return(rf)
}

hf <- flux_dt[siteID == "US-Ha1"]

# If I fix the base temperature to 0
t_base <- 0
temp <- hf[PhenoYear == 2001]
sos <- unique(temp$SOS)
t0 <- seq(12, 365, by = 5)


# fig: TT problem 
png("Output/Fig/TT_problem.png", res = 300, width = 2200, height = 2000)

# png("Output/Fig/TT_problem_dark.png", res = 300, width = 2200, height = 2000)
# par(bg = NA, fg = "white", col.axis = "white", col.lab = "white",
#     col.main = "white", col.sub = "white", cex.lab = 1.5)

par(mfrow = c(2, 1), lwd = 3, cex = 1.5)
# Plot possibilities
par(mar = c(0, 3, 1, 1), mgp = c(1.5, 0.5, 0))
plot(flux_li$doy, temp$Tmean, type = "l", lwd = 1, col = "grey", 
    ylab = "Temperature", xaxt = "n")

# T_base
T_base <- 0
abline(h = T_base, lty = 2)
text(-100, T_base, label = bquote(T[base]), pos = 3)

# SOS
abline(v = sos, lty = 2, col = "seagreen", lwd = 2)


par(new = TRUE)
plot(NA, xlim = range(flux_li$doy), ylim = c(0, 3000), type = "l", 
    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# plot(NA, xlim = c(150, 200), ylim = c(0, 3000), type = "l")
# colors <- viridis::viridis(length(t0))
greys <- rev(RColorBrewer::brewer.pal(8, "Greys"))
colors <- colorRampPalette(c(greys[1], color_frc, greys[1]))(length(t0))
f_crit <- NULL
for (i in seq_along(t0)) {
    frc <- Forcing(t0[i])
    lines(flux_li$doy, frc, type = "l", col = colors[i])
    # abline(h = frc[sos + 110])
    f_crit <- c(f_crit, frc[sos + 110])
}

# If I estimate F_crit for this site year
par(mar = c(3, 3, 0, 1))
# Convert t0 to real date
t0_date <- as.Date("2001-01-01") + flux_li$doy[t0]
plot(t0_date, f_crit, type = "l", 
    xlab = bquote(t["0"]), ylab = bquote(F["crit"]))
# SOS
sos_date <- as.Date("2001-01-01") + sos
abline(v = sos_date, lty = 2, col = "seagreen", lwd = 2)
text(sos_date, grconvertY(0.9, "npc"),
    label = "Start-of-season",
    col = "seagreen", xpd = NA, pos = 4
)


dev.off()