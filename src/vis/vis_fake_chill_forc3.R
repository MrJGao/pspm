# ******************************************************************************
# See how chilling and forcing calculations would affect significance proportion
# of the chilling-forcing relationship
# 
# Author: Xiaojie Gao
# Date: 2024-08-16
# ******************************************************************************
source("src/base.R")
LoadHelpers()
library(data.table)
library(magrittr)
library(lubridate)


args <- commandArgs(trailingOnly = TRUE)
idx <- as.numeric(args[1])


# Calculate chilling
# As in Fu et al 2015, chilling was calculated as chill days during which the
#   daily mean temperature is between 0 and 5, start counting from Nov. 1st and
#   end by the average SOS date.
ChillForc <- function(pp_dt, tc = -30, Tb = 5, t0 = 30) {
    # Calculate chilling
    # Let's use chill days for now
    pp_dt[, ":="(chill = 0, acc_chill = 0)]
    pp_dt[doy > tc,
        chill := .SD[, ChillFix(Tmean, "d")],
        by = .(id)
    ]
    # Accumulated chilling
    pp_dt[doy > tc,
        acc_chill := .SD[, cumsum(chill)],
        by = .(id)
    ]
    # pp_dt[doy <= tc, ":=" (chill = NA, acc_chill = NA)]

    # Calculate forcing
    pp_dt[, ":="(forc = 0, acc_forc = 0)]
    pp_dt[doy > t0,
        forc := .SD[
            ,
            ifelse(Tmean - Tb > 0, Tmean - Tb, 0)
        ],
        by = .(id)
    ]
    # Accumulated forcing
    pp_dt[doy > t0,
        acc_forc := .SD[, cumsum(forc)],
        by = .(id)
    ]
    # pp_dt[doy <= t0, ":="(forc = NA, acc_forc = NA)]

    return(pp_dt)
}


FitLmPerSite <- function(st) {
    fit <- lm(acc_forc ~ acc_chill, data = st)
    f <- summary(fit)$fstatistic
    if (is.null(f) | nrow(st) < 40) {
        # browser() # For debug
        return(data.table(r2 = 9999, pval = 9999, inter = 9999, slp = 9999))
    }
    p_val <- pf(f[1], f[2], f[3], lower.tail = FALSE)
    r_sqr <- round(summary(fit)$r.squared, 2)

    inter <- coef(fit)[1]
    slp <- coef(fit)[2]

    return(data.table(r2 = r_sqr, pval = p_val, inter = inter, slp = slp))
}


outdir <- "pipe/fake_chill_forc3"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# sim_spec_dt <- readRDS("pipe/sim_spec_aesculus_dt_est_fcrit.Rds")
# sim_spec_dt <- readRDS("pipe/sim_spec_aesculus_dt_est_fcrit.Rds")

sim_spec_dt <- readRDS("pipe/sim_spec_aesculus_dt_by_at_1.Rds")

sim_spec_dt <- ToPersonPeriod(sim_spec_dt)


set.seed(idx)
tc <- sample(seq(-60, -1), 1)
Tb <- sample(seq(0, 10), 1)
t0 <- sample(seq(0, 30), 1)
sim_spec_dt <- ChillForc(sim_spec_dt, tc = tc, Tb = Tb, t0 = t0)


sim_spec_dt[, avg_sos := mean(SOS), by = siteID]
# For each site, only account chilling before the average SOS date
chill_forc_dt <- by(sim_spec_dt, sim_spec_dt$id, function(st) {
    st_dt <- st[
        which.min(abs(doy - avg_sos)),
        .(siteID, PhenoYear, acc_chill, acc_forc, SOS)
    ]
    return(st_dt)
})
chill_forc_dt <- do.call(rbind, chill_forc_dt)

fit_dt <- chill_forc_dt[, FitLmPerSite(.SD), by = siteID]
fit_dt <- fit_dt[pval != 9999, ]

# Significant fit percentage
sig_fit_dt <- fit_dt[pval < 0.05 & slp < 0, ]
sig_pct <- round(nrow(sig_fit_dt) / nrow(fit_dt) * 100, 1)

# rm(sim_spec_dt)

res <- data.table(
    tc = tc, Tb = Tb, t0 = t0,
    sig_pct = sig_pct
)


saveRDS(res, file = file.path(outdir, paste0(idx, ".Rds")))

