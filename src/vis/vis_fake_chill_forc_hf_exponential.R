library(data.table)
library(sp)
library(lubridate)

source("src/base.R")
source("src/vis/vis_base.R")

dt <- fread("data/flux_mcd12q2_daymet_dt.csv")

T_base <- 2
t_c <- "-01-01"
t_f <- "-02-01"

site_dt <- dt[siteID == "US-Ha1"]
frc_chil_dt <- by(site_dt, list(site_dt$siteID, site_dt$PhenoYear), function(sy) {
    yr <- as.character(sy[1, PhenoYear])
    sos_date <- as.Date(paste0(yr, "-01-01")) + sy[1, SOS]
    # fake sos
    sos_date <- as.Date(paste0(yr, "-01-01")) + 80
    frc <- sy[between(Date, as.Date(paste0(yr, t_f)), sos_date), Tmean]
    frc <- ifelse(frc > T_base, frc, 0)
    frc <- sum(frc)

    prev_yr <- as.numeric(yr)
    chil <- sy[between(Date, as.Date(paste0(prev_yr, t_c)), sos_date), Tmean]
    chil <- ifelse(chil <= T_base, 1, 0)
    chil <- sum(chil)

    return(data.table(frc = frc, chil = chil))
})
frc_chil_dt <- do.call(rbind, frc_chil_dt)

# fig: Chiling-forcing negative exponential
png("out/chil_frc_exp.png",
    width = 1600, height = 1600, res = 300
)
par(cex.lab = 1.5, mgp = c(1.5, 0.5, 0))

plot(frc_chil_dt[, .(chil, frc)],
    pch = 16,
    xlab = "Chill days", ylab = "Forcing"
)
text(grconvertX(0.5, "npc"), grconvertY(0.9, "npc"),
    labels = bquote(T[base] == .(T_base) ~ degree ~ C), pos = 4, cex = 1.5
)
text(grconvertX(0.5, "npc"), grconvertY(0.8, "npc"),
    labels = bquote(t["0_chil"] == "Jan 1st"), pos = 4, cex = 1.5
)
text(grconvertX(0.5, "npc"), grconvertY(0.7, "npc"),
    labels = bquote(t["0_frc"] == "Feb 1st"), pos = 4, cex = 1.5
)
newx <- seq(range(frc_chil_dt$chil)[1], range(frc_chil_dt$chil)[2], length.out = 100)
exp_fit <- lm(log(frc) ~ chil, data = frc_chil_dt)
pred <- exp(predict(exp_fit, newdata = data.frame(chil = newx)))
lines(newx, pred, col = "red")

dev.off()