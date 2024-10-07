# ******************************************************************************
# Vis simulations to show that recent results about SOS advancing and chilling
# decreasing can be found in simulated data w/o chilling as well. And, the
# negative chilling-forcing relationship may be an artifact.
# 
# Data is generated from `src/mod_fake_chill_forc_pep.R`
# ******************************************************************************
rm(list=ls())
source("src/base.R")
LoadHelpers()
library(data.table)
library(magrittr)
library(lubridate)



# sim_spec_dt <- readRDS("pipe/sim_spec_aesculus_dt.Rds")
# sim_spec_dt <- readRDS("pipe/sim_spec_aesculus_dt_est_fcrit.Rds")
sim_spec_dt <- readRDS("pipe/sim_spec_aesculus_dt_by_at_1.Rds")

periods_1 <- c(1980, 1994, 1999, 2013)
periods_2 <- c(1951, 1979, 1980, 1999, 2000, 2019)

# ~ SOS and Chilling trends ####
# ~ ----------------------------------------------------------------------------
pdf("out/sos_trends3.pdf", width = 8)
# png("out/sos_trends.png", width = 1300, height = 1200, res = 150)

par(mfrow = c(2, 2), mgp = c(1.5, 0.5, 0))

# Split method 1
sim_dt_p1 <- sim_spec_dt[PhenoYear > periods_1[1] & PhenoYear < periods_1[2]]
sim_dt_p2 <- sim_spec_dt[PhenoYear > periods_1[3] & PhenoYear < periods_1[4]]

cols <- hcl.colors(8, "dynamic")[c(1, 5)]
hist(sim_dt_p1$SOS, breaks = 50, 
    col = adjustcolor(cols[1], 0.5),
    freq = FALSE, 
    xlab = "SOS",
    xlim = c(50, 250),
    ylim = c(0, 0.05),
    border = adjustcolor("grey", 0.5),
    main = ""
)
hist(sim_dt_p2$SOS, breaks = 50, col = adjustcolor(cols[2], 0.5), 
    freq = FALSE,
    border = adjustcolor("grey", 0.5),
    add = TRUE
)
legend("topright", bty = "n", fill = cols, legend = c("1980-1994", "1999-2013"))


# Split method 2
sim_dt_p1 <- sim_spec_dt[PhenoYear > periods_2[1] & PhenoYear < periods_2[2]]
sim_dt_p2 <- sim_spec_dt[PhenoYear > periods_2[3] & PhenoYear < periods_2[4]]
sim_dt_p3 <- sim_spec_dt[PhenoYear > periods_2[5] & PhenoYear < periods_2[6]]

cols <- hcl.colors(8, "dynamic")[c(1, 5, 8)]
boxplot(sim_dt_p1$SOS, sim_dt_p2$SOS, sim_dt_p3$SOS,
    outline = FALSE,
    names = c("1951-1979", "1980-1999", "2000-2019"),
    ylab = "SOS", xlab = "Period"
)
# legend("topright", bty = "n", fill = cols, 
#     legend = c("1951-1979", "1980-1999", "2000-2019")
# )


# Calculate chilling
# As in Fu et al 2015, chilling was calculated as chill days during which the 
#   daily mean temperature is between 0 and 5, start counting from Nov. 1st and
#   end by the average SOS date.

ChillForc <- function(pp_dt, tc = -30, Tb = 5, t0 = 30) {
    # Calculate chilling
    # Let's use chill days for now
    pp_dt[, ":=" (chill = 0, acc_chill = 0)]
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

sim_spec_dt <- ToPersonPeriod(sim_spec_dt)
sim_spec_dt <- ChillForc(sim_spec_dt)

sim_spec_dt[, avg_sos := mean(SOS), by = siteID]
# For each site, only account chilling before the average SOS date
chill_forc_dt <- by(sim_spec_dt, sim_spec_dt$id, function(st) {
    st_dt <- st[which.min(abs(doy - avg_sos)), 
        .(siteID, PhenoYear, acc_chill, acc_forc, SOS)
    ]
    return(st_dt)
})
chill_forc_dt <- do.call(rbind, chill_forc_dt)





# Split method 1
cf_p1 <- chill_forc_dt[PhenoYear > periods_1[1] & PhenoYear < periods_1[2]]
cf_p2 <- chill_forc_dt[PhenoYear > periods_1[3] & PhenoYear < periods_1[4]]

cols <- hcl.colors(8, "dynamic")[c(1, 5)]
barplot(c(mean(cf_p1$acc_chill), mean(cf_p2$acc_chill)), 
    names = c("1980-1994", "1999-2013"), col = cols,
    ylim = c(0, 80)
)
err_1 <- mean(cf_p1$acc_chill) - sd(cf_p1$acc_chill) * c(-1, 1)
err_2 <- mean(cf_p2$acc_chill) - sd(cf_p2$acc_chill) * c(-1, 1)
segments(0.7, err_1[1], 0.7, err_1[2], lwd = 2)
segments(1.9, err_2[1], 1.9, err_2[2], lwd = 2)
legend("topright", bty = "n", fill = cols, legend = c("1980-1994", "1999-2013"))

# Split method 2
cf_p1 <- chill_forc_dt[PhenoYear > periods_2[1] & PhenoYear < periods_2[2]]
cf_p2 <- chill_forc_dt[PhenoYear > periods_2[3] & PhenoYear < periods_2[4]]
cf_p3 <- chill_forc_dt[PhenoYear > periods_2[5] & PhenoYear < periods_2[6]]

cols <- hcl.colors(8, "dynamic")[c(1, 5, 8)]
boxplot(cf_p1$acc_chill, cf_p2$acc_chill, cf_p3$acc_chill, 
    outline = FALSE, names = c("1951-1979", "1980-1999", "2000-2019"),
    ylab = "Chilling accumulation", xlab = "Period"
)
# legend("topright", bty = "n", fill = cols, 
#     legend = c("1951-1979", "1980-1999", "2000-2019")
# )

text(
    grconvertX(c(0.03, 0.53, 0.03, 0.53), "ndc"), 
    grconvertY(c(0.92, 0.92, 0.43, 0.43), "ndc"),
    label = letters[1:4],
    xpd = NA, cex = 1.5
)

dev.off()





# ~ Neg. chilling-forcing relationship ####
# ~ ----------------------------------------------------------------------------

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

fit_dt <- chill_forc_dt[, FitLmPerSite(.SD), by = siteID]
fit_dt <- fit_dt[pval != 9999, ]

# Significant fit percentage
sig_fit_dt <- fit_dt[pval < 0.05 & slp < 0,]
sig_pct <- round(nrow(sig_fit_dt) / nrow(fit_dt) * 100, 1)

# Average stats from significant fits
avg_inter <- mean(sig_fit_dt$inter, na.rm = TRUE)
avg_slp <- mean(sig_fit_dt$slp, na.rm = TRUE)
avg_r2 <- round(mean(sig_fit_dt$r2, na.rm = TRUE), 2)


pdf("out/fake_chill_forc3.pdf", width = 4, height = 4)
# png("out/fake_chill_forc.png", width = 600, height = 600, res = 150)

par(mfrow = c(1, 1), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 2, 2))
smoothScatter(chill_forc_dt[, .(acc_chill, acc_forc)], 
    ylim = c(20, 500), xlim = c(0, 120), 
    pch = 16, cex = 0.5, 
    colramp = colorRampPalette(
        c("white", rev(hcl.colors(4, "gnbu"))), 
        alpha = TRUE
    ),
    xlab = "Accumulated chilling", ylab = "Accumulated forcing",
    nrpoints = 0
)
# All slope
null <- apply(sig_fit_dt, 1, function(a_row) {
    abline(a = a_row[["inter"]], b = a_row[["slp"]], 
        lwd = 0.3, col = adjustcolor("grey", 0.4))
})
# The average slope
abline(a = avg_inter, b = avg_slp, lwd = 2)

# legend("topright", legend = bquote(Avg.R^2 == .(avg_r2)), 
#     bty = "n", y.intersp = 3)
legend("topright", legend = c(paste0("Sig.% = ", sig_pct)), bty = "n")

dev.off()


