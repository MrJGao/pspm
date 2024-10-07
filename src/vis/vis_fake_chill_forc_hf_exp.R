# ******************************************************************************
# Fake negative exponential chilling-forcing relationship at Harvard Forest
# ******************************************************************************
library(data.table)
library(sp)
library(lubridate)

source("src/base.R")
source("src/vis/vis_base.R")

hf_dt <- fread("data/hf_grd_temperature_pheno.csv")
# hf_dt <- fread("data/flux_mcd12q2_daymet_dt.csv")[siteID == "US-Ha1"]

# hf_dt <- fread("data/raw/hf300-05-daily-m.csv")
# hf_dt <- hf_dt[between(date, as.Date("1990-01-01"), as.Date("2018-12-31"))]

# hf_dt[, SOS := 140]
# hf_dt[, PhenoYear := year(date) + 1]
# hf_dt[, ":=" (siteID = "HF", Tmean = airt, Tmax = airtmax, Tmin = airtmin)]

# hf_li <- FormatDataForPhenoModel(hf_dt)
# hf_dt <- Conv2Dt(hf_li)

# hf_dt[, Temp := Tmean]


DoFit <- function(T_base, t_c, t_f) {
    site_dt <- hf_dt
    frc_chil_dt <- by(site_dt, site_dt$PhenoYear, function(sy) {
        yr <- as.character(sy[1, PhenoYear])
        # sos_date <- as.Date(paste0(yr, "-01-01")) + sy[1, SOS]

        # fake sos, can be changed to other values
        sos_date <- as.Date(paste0(yr, "-01-01")) + 140
        
        frc <- sy[between(Date, as.Date(paste0(yr, t_f)), sos_date), Temp]
        frc <- ifelse(frc > T_base, frc, 0)
        frc <- sum(frc)

        prev_yr <- ifelse (as.numeric(substr(t_c, 2, 3)) > 8, 
            as.numeric(yr) - 1, 
            as.numeric(yr)
        )
        chil <- sy[between(Date, as.Date(paste0(prev_yr, t_c)), sos_date), Temp]
        chil <- ifelse(chil <= T_base, 1, 0)
        chil <- sum(chil)

        return(data.table(frc = frc, chil = chil))
    })
    frc_chil_dt <- do.call(rbind, frc_chil_dt)

    exp_fit <- lm(log(frc) ~ chil, data = frc_chil_dt)
    if (is.na(coef(exp_fit)[2]) || coef(exp_fit)[2] > -0.025) {
        return(NULL)
    }

    par(cex.lab = 1.5, mgp = c(1.5, 0.5, 0))
    plot(frc_chil_dt[, .(chil, frc)], pch = 16,
        xlab = "Chill days", ylab = "Forcing"
    )
    text(grconvertX(0.5, "npc"), grconvertY(0.9, "npc"), 
        labels = bquote(T[base] == .(T_base)~degree~C), pos = 4, cex = 1.5
    )
    text(grconvertX(0.5, "npc"), grconvertY(0.8, "npc"), 
        # labels = bquote(t["0_chil"] == "Jan 1st"), 
        labels = bquote(t["0_chil"] == .(t_c)), 
        pos = 4, cex = 1.5
    )
    text(grconvertX(0.5, "npc"), grconvertY(0.7, "npc"), 
        # labels = bquote(t["0_frc"] == "Feb 1st"), 
        labels = bquote(t["0_frc"] == .(t_f)), 
        pos = 4, cex = 1.5
    )
    # newx <- seq(
    #     range(frc_chil_dt$chil)[1], 
    #     range(frc_chil_dt$chil)[2], 
    #     length.out = 100
    # )
    newx <- seq(
        0, 
        300, 
        by = 1
    )
    exp_fit <- lm(log(frc) ~ chil, data = frc_chil_dt)
    pred <- exp(predict(exp_fit, newdata = data.frame(chil = newx)))
    lines(newx, pred, col = "red")
}

# T_base <- 2
# t_c <- "-01-01"
# t_f <- "-02-01"


# The following can be used to see this is not a unique case

pdf("out/zzz_hf.pdf", width = 5, height = 5) 
for (i in 1:1000) {
    # The variables are T_base, t_c, t_f, 
    T_base <- round(runif(1, 0, 10))
    tf_date <- sample(seq(1, 90, by = 5), 1) %>%
        as.Date(origin = "2000-01-01")
    tc_date <- tf_date - sample(seq(1, 60, by = 5), 1)

    t_f <- substr(tf_date, 5, 10)
    t_c <- substr(tc_date, 5, 10)

    fake_sos_date <- 130

    DoFit(T_base, t_c, t_f)
}

dev.off()


# # fig: Chiling-forcing negative exponential
# png("out/chil_frc_exp_hf.png", 
#     width = 1600, height = 1600, res = 300
# )
# # png("Output/Assets/chil_frc_exp_dark.png", 
# #     width = 1600, height = 1600, res = 300
# # )
# # par(bg = NA, fg = "white", col.axis = "white", col.lab = "white", 
# #     col.main = "white", col.sub = "white", lwd = 3)

# DoFit(T_base, t_c, t_f)

# dev.off()


