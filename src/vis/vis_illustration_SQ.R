#'******************************************************************************
#' Description: Sequential model illustration.
#'******************************************************************************
rm(list=ls())

source("src/base.R")
source("src/vis/vis_base.R")
source("src/hlp/hlp_spring_pheno_models.R")

library(data.table)



# in: Use flux tower temperature as an example.
temp_dt <- fread(
    file.path(gdir, "Data/ARD/flux_measured_daily_temperature.csv")
)
temp_dt$Date <- as.Date(temp_dt$Date)
temp_dt <- temp_dt[siteID == "US-Ha1" &
    between(Date, as.Date("2001-09-01"), as.Date("2002-08-31")), ]



# Fake parameters
# temp_dt <- data.table(Date = -110:254, temp = rnorm(365, 0, 20))
sos <- as.Date("2002-06-10")
t0_chill <- as.Date("2001-11-01")
t0 <- as.Date("2002-1-1")
T_base <- 0
T_opt <- 2
T_min <- -2
T_max <- 5


# Chilling
Rc <- triangular_temperature_response(
    T = temp_dt[between(Date, t0_chill, t0), Temp], 
    T_opt = T_opt, T_min = T_min, T_max = T_max
)
Sc <- cumsum(Rc)

# Forcing
temp_dt[between(Date, as.Date("2001-09-01"), t0), tt := 0]
temp_dt[between(Date, t0, as.Date("2002-08-31")), 
    tt := ifelse(Temp > 0, Temp, 0)]
temp_dt[, frc := cumsum(tt)]
temp_dt[, tt := NULL]


# fig: SQ illustration
png("Output/Fig/SQ_illustration.png", width = 2560, height = 1440, res = 300)

# png("Output/Fig/SQ_illustration_dark.png", 
#     width = 2560, height = 1440, res = 300
# )
# par(bg = NA, fg = "white", col.axis = "white", col.lab = "white", 
#     col.main = "white", col.sub = "white")

par(mar = c(3, 3, 1, 1), mgp = c(0.5, 0.5, 0), lwd = 3, cex = 1.5)

# Plot temperature
plot(temp_dt[, .(Date, Temp)], type = "l", 
    xlab = "Date", ylab = "",
    col = "grey", bty = "L", lwd = 1,
    xaxt = "n", yaxt = "n")

# T_base
abline(h = T_base, lty = 3)
text(t0 - 100, T_base, label = bquote(T[base]), pos = 3)

# Chiling state
par(new = TRUE)
plot(NA, xlim = range(temp_dt$Date), ylim = 2 * range(Sc), 
    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(t0_chill:t0, Sc, lwd = 3, col = color_chi)
segments(t0_chill, grconvertY(0, "npc"), t0_chill, grconvertY(0.05, "npc"))
text(t0_chill, grconvertY(0.05, "npc"), label = bquote(t["0_chill"]), pos = 2)

segments(t0, max(Sc), t0 + 20, max(Sc), col = color_chi, lwd = 2, lty = 2)
text(t0 + 10, max(Sc), label = bquote(C["req"]), pos = 3, col = color_chi)

segments(t0, max(Sc), t0, grconvertY(0, "npc"), lty = 2)
text(t0, grconvertY(0.05, "npc"), label = bquote(t["0"]), pos = 2)


# Forcing state
par(new = TRUE)
plot(NA, xlim = range(temp_dt$Date), ylim = range(temp_dt$frc),
    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(temp_dt[Date > t0, .(Date, frc)],
    type = "l", lwd = 3, col = color_frc)
segments(sos, temp_dt[Date == sos, frc], grconvertX(1, "npc"), 
    temp_dt[Date == sos, frc], 
    lty = 3, lwd = 3, col = color_frc)
text(sos + 30, temp_dt[Date == sos, frc], pos = 3, labels = bquote(F[crit]), 
    col = color_frc)
# SOS
segments(sos, temp_dt[Date == sos, frc], sos, grconvertY(0, "npc"),
    lwd = 3, lty = 2, col = "seagreen"
)
text(sos, grconvertY(0, "npc", "user"),
    label = "Start-of-season",
    col = "seagreen", xpd = NA, pos = 1
)



# Legend
legend("topleft",
    lty = c(1, 1, 1),
    col = c("grey", color_chi, color_frc),
    lwd = c(1, 3, 3),
    legend = c("Temperature", "Chilling", "Forcing"),
    bty = "n", cex = 0.8
)
title("Sequential model", line = -1)

dev.off()


