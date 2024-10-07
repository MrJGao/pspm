#'******************************************************************************
#' Description: Parallel model illustration.
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
T_opt <- -2
T_min <- -4
T_max <- 10
C_ini <- 0.01



# Chilling
Rc <- triangular_temperature_response(T = temp_dt[Date > t0_chill, Temp], 
    T_opt = T_opt, T_min = T_min, T_max = T_max)
temp_dt[, Sc := 0]
temp_dt[Date > t0_chill, Sc := cumsum(Rc)]

C_req <- 80

temp_dt[, k := ifelse(Sc < C_req, C_ini + Sc * (1 - C_ini) / C_req, 1)]




# Forcing
temp_dt[between(Date, as.Date("2001-09-01"), t0), tt := 0]
temp_dt[, Rf := 0]
temp_dt[between(Date, t0, as.Date("2002-08-31")), 
    Rf := ifelse(Temp > 0, Temp, 0)]
temp_dt[, frc := cumsum(k * Rf)]
temp_dt[, frc_no_chil := cumsum(Rf)]


# fig: PA illustration
png("Output/Fig/PA_illustration.png", width = 2560, height = 1440, res = 300)

# png("Output/Fig/PA_illustration_dark.png", 
#     width = 2560, height = 1440, res = 300
# )
# par(bg = NA, fg = "white", col.axis = "white", col.lab = "white", 
#     col.main = "white", col.sub = "white", lwd = 3, cex = 1.5)

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
plot(NA, xlim = range(temp_dt$Date), ylim = range(temp_dt$Sc) * 1.1, 
    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(temp_dt[Date > t0_chill, .(Date, Sc)], lwd = 2, col = color_chi)
segments(t0_chill, grconvertY(0, "npc"), t0_chill, grconvertY(0.05, "npc"))
text(t0_chill, grconvertY(0.05, "npc"), label = bquote(t["0_chill"]), pos = 2)

segments(temp_dt[Sc >= C_req, Date][1], C_req, temp_dt[Sc >= C_req, Date][30], 
    C_req, col = color_chi, lwd = 2, lty = 2)
text(temp_dt[Sc >= C_req, Date][15], C_req, label = bquote(C["req"]), 
    pos = 1, col = color_chi)

# Forcing state
par(new = TRUE)
plot(NA, xlim = range(temp_dt$Date), ylim = range(temp_dt$frc),
    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# Forcing w/o chilling
lines(temp_dt[Date > t0, .(Date, frc_no_chil)],
    type = "l", lwd = 2, col = Transparent(color_frc, 0.3))
# Forcing w/ chilling
lines(temp_dt[Date > t0, .(Date, k * frc)], lwd = 3, col = color_frc)
segments(sos, temp_dt[Date == sos, k * frc], grconvertX(1, "npc"), 
    temp_dt[Date == sos, k * frc], 
    lty = 3, lwd = 3, col = color_frc)
text(sos + 30, temp_dt[Date == sos, k * frc], pos = 3, labels = bquote(F[crit]), 
    col = color_frc)
segments(t0, grconvertY(0.05, "npc"), t0, grconvertY(0, "npc"))
text(t0, grconvertY(0.05, "npc"), label = bquote(t["0"]), pos = 2)



# SOS
segments(sos, temp_dt[Date == sos, k * frc], sos, grconvertY(0, "npc"),
    lwd = 4, lty = 2, col = "seagreen"
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
title("Parallel model", line = -1)

dev.off()


