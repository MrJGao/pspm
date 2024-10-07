#'******************************************************************************
#' Description: Thermal time model machanism illustration.
#'******************************************************************************
rm(list=ls())

source("src/base.R")
source("src/vis/vis_base.R")
source("src/hlp/hlp_spring_pheno_models.R")

# library(ggplot2)
# library(egg)
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
t0 <- as.Date("2002-1-1")
T_base <- 0

temp_dt[between(Date, as.Date("2001-09-01"), t0), tt := 0]
temp_dt[between(Date, t0, as.Date("2002-08-31")), 
    tt := ifelse(Temp > 0, Temp, 0)]
temp_dt[, frc := cumsum(tt)]
temp_dt[, tt := NULL]

# Base temperature figure
# p_base <- ggplot(temp_dt) +
#     geom_line(aes(x = Date, y = frc)) +
#     xlab("Date") +
#     ylab(expression("Temperature ("~degree~C~")")) +
#     theme_classic() +
#     theme(
#         axis.ticks = element_blank(),
#         axis.text = element_blank()
#     )
# p_base

# p_base + 
#     # t_0
#     geom_segment(aes(x = t0, y = -500, xend = t0, yend = -600),
#         lineend = "round", linejoin = "round") +
#     geom_text(aes(x = t0, y = -500), 
#         label = "t[0]", vjust = -0.5, parse = TRUE
#     ) +
#     # F_crit
#     geom_segment(aes(x = sos, y = temp_dt[sos, frc], 
#         xend = 254, yend = temp_dt[sos, frc]), linetype = "dashed") +
#     geom_text(aes(x = sos + 50, y = temp_dt[sos, frc]), 
#         label = "F[crit]", parse = TRUE, vjust = -0.5) +
#     # SOS
#     geom_segment(aes(x = sos, y = -600, xend = sos, yend = 0), 
#         linetype = "dashed", size = 1.5, color = "seagreen") +
#     geom_text(aes(x = sos, y = -300), label = "start of season", 
#         hjust = 1.1, color = "seagreen")




# fig: TT illustration
png("Output/Fig/TT_illustration.png", width = 2560, height = 1440, res = 300)

# png("Output/Fig/TT_illustration_dark.png", 
#     width = 2560, height = 1440, res = 300
# )
# par(bg = NA, fg = "white", col.axis = "white", col.lab = "white", 
#     col.main = "white", col.sub = "white")

par(mar = c(3, 3, 1, 1), mgp = c(0.5, 0.5, 0), lwd = 3, cex = 1.5)
plot(temp_dt[, .(Date, Temp)], type = "l", 
    xlab = "Date", ylab = "",
    col = "#dfdddd", bty = "L",
    xaxt = "n", yaxt = "n", lwd = 1)
# t0
segments(t0, -5, t0, grconvertY(0, "npc"))
text(t0, -5, label = bquote(t[0]), pos = 2)
# T_base
abline(h = T_base, lty = 3)
text(t0 - 100, T_base, label = bquote(T[base]), pos = 3)

# Forcing
par(new = TRUE)
plot(temp_dt[, .(Date, frc)], type = "l", lwd = 3, col = color_frc,
    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# SOS
segments(sos, temp_dt[Date == sos, frc], sos, grconvertY(0, "npc"), 
    lwd = 3, lty = 2, col = "seagreen")
# legend(sos, 400, legend = "start-of-season", text.col = "seagreen", 
    # box.col = "white", box.lwd = 0, bg = "white", xjust = 0.54, cex = 0.7)
text(sos, grconvertY(0, "npc", "user"), label = "Start-of-season", 
    col = "seagreen", xpd = NA, pos = 1)

# Forcing critical value
segments(sos, temp_dt[Date == sos, frc], grconvertX(1, "npc"), 
    temp_dt[Date == sos, frc], 
    lty = 3, lwd = 3, col = color_frc)
text(sos + 30, temp_dt[Date == sos, frc], pos = 3, labels = bquote(F["crit"]), 
    col = color_frc)

# Legend
legend("topleft", 
    lty = c(1, 1), 
    col = c("grey", color_frc), 
    lwd = c(1, 3),
    legend = c("Temperature", "Forcing"), 
    bty = "n", cex = 0.8)
title("Thermal Time model", line = -1)

dev.off()
