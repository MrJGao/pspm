#'******************************************************************************
#' Description: Alternating model illustration.
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
a <- 100
b <- 400
c <- -0.2


# Chilling
temp_dt[, Rc := 0]
temp_dt[Date > t0_chill, Rc := ifelse(Temp < T_base, 1, 0)]
temp_dt[, Sc := 0]
temp_dt[Date > t0_chill, Sc := cumsum(Rc)]


# Forcing
temp_dt[, Rf := 0]
temp_dt[
    between(Date, t0, as.Date("2002-08-31")), 
    Rf := ifelse(Temp > 0, Temp, 0)
]
temp_dt[, frc := cumsum(Rf)]

F_crit <- a + b * exp(c * temp_dt$Sc)





# old_par <- par()
# par(old_par)

# fig: AT illustration
png("Output/Fig/AT_illustration.png", width = 2560, height = 1440, res = 300)

# png("Output/Fig/AT_illustration_dark.png", 
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
plot(NA, xlim = range(temp_dt$Date), ylim = range(temp_dt$Sc) * 1.1, 
    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(temp_dt[Date > t0_chill, .(Date, Sc)], lwd = 3, col = color_chi)
segments(t0_chill, grconvertY(0, "npc"), t0_chill, grconvertY(0.05, "npc"))
text(t0_chill, grconvertY(0.05, "npc"), label = bquote(t["0_chill"]), pos = 2)


# Forcing state
par(new = TRUE)
plot(NA, xlim = range(temp_dt$Date), ylim = range(temp_dt$frc),
    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# Forcing w/ chilling
lines(temp_dt[Date > t0, .(Date, frc)], lwd = 3, col = color_frc)
segments(sos, temp_dt[Date == sos, frc], grconvertX(1, "npc"), 
    temp_dt[Date == sos, frc], 
    lty = 3, lwd = 3, col = color_frc)
text(sos + 30, temp_dt[Date == sos, frc], pos = 3, labels = bquote(F[crit]), 
    col = color_frc)
segments(t0, grconvertY(0.05, "npc"), t0, grconvertY(0, "npc"))
text(t0, grconvertY(0.05, "npc"), label = bquote(t["0"]), pos = 2)

# Initial F_crit
segments(sos + 30, temp_dt[Date == sos + 30, frc], grconvertX(1, "npc"), 
    temp_dt[Date == sos + 30, frc], 
    lty = 3, lwd = 3, col = Transparent(color_frc, 0.5))
text(sos + 30, temp_dt[Date == sos + 30, frc], adj = c(-1, -0.5),
    labels = bquote(F["crit"]), col = Transparent(color_frc, 0.5))



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
title("Alternating model", line = -1)


# Subplot
par(new = TRUE, fig = c(0.15, 0.35, 0.45, 0.75), mar = c(1, 2, 0.5, 0.5))
rect(
    grconvertX(par("fig")[1], "ndc"), 
    grconvertY(par("fig")[3], "ndc"),
    grconvertX(par("fig")[2], "ndc"), 
    grconvertY(par("fig")[4], "ndc"),
    lwd = 0,
    col = Transparent("white", 0.8), xpd = NA
)

# For dark background
# rect(
#     grconvertX(par("fig")[1], "ndc"), 
#     grconvertY(par("fig")[3], "ndc"),
#     grconvertX(par("fig")[2], "ndc"), 
#     grconvertY(par("fig")[4], "ndc"),
#     lwd = 0,
#     col = Transparent("black", 0.8), xpd = NA
# )
par(new = TRUE, fig = c(0.15, 0.35, 0.45, 0.75), mar = c(1, 1.2, 0.5, 0.5))
plot(NA,
    xlim = range(temp_dt$Sc), ylim = range(F_crit),
    xlab = "", ylab = "",
    bty = "L", xaxt = "n", yaxt = "n"
)
lines(temp_dt$Sc, F_crit, col = color_frc)
mtext(side = 1, text = "Chill days", col = color_chi)
mtext(side = 2, text = bquote(F["crit"]), adj = 1.1, col = color_frc, las = 1)

dev.off()


