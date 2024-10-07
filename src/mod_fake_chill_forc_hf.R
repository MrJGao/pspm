# ******************************************************************************
# Fake chilling-forcing relationship simulation by using HF+ground temperature.
# ******************************************************************************
source("src/base.R")
LoadHelpers()
library(data.table)
library(magrittr)
library(lubridate)


hf_dt <- fread("data/hf_grd_temperature_pheno.csv")

sim_dt <- lapply(1:1000, function(i) {
    # The variables are T_base, t_c, t_f, 
    T_base <- runif(1, 0, 10)
    tf_date <- sample(seq(1, 90, by = 5), 1) %>%
        as.Date(origin = "2000-01-01")
    tc_date <- tf_date - sample(seq(1, 90, by = 5), 1)

    t_f <- substr(tf_date, 6, 10)
    t_c <- substr(tc_date, 6, 10)

    # Could try different ranges
    # fake_sos_date <- runif(1, 60, 200)
    fake_sos_date <- runif(1, 120, 160)

    frc_chil_dt <- by(hf_dt, list(hf_dt$siteID, hf_dt$PhenoYear), function(sy) {
        yr <- as.character(sy[1, PhenoYear])
        # sos_date <- as.Date(paste0(yr, "-01-01")) + sy[1, SOS]
        sos_date <- as.Date(paste0(yr, "-01-01")) + fake_sos_date
        frc <- sy[
            between(Date, as.Date(paste0(yr, "-", t_f)), sos_date), 
            Temp
        ]
        frc <- ifelse(frc > T_base, frc, 0)
        frc <- sum(frc)

        prev_yr <- as.numeric(yr) - 1
        chil <- sy[
            between(Date, as.Date(paste0(prev_yr, "-", t_c)), sos_date), 
            Temp
        ]
        chil <- ifelse(chil <= T_base, 1, 0)
        chil <- sum(chil)

        return(data.table(frc = frc, chil = chil))
    })
    frc_chil_dt <- do.call(rbind, frc_chil_dt)

    return(frc_chil_dt)
})


saveRDS(sim_dt, "pipe/hf_chil_forc_sim.Rds")






