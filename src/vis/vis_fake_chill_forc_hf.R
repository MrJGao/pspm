# ******************************************************************************
# Visualize fake chilling and forcing relationship using HF data.
# ******************************************************************************


sim_dt <- readRDS("pipe/hf_chil_forc_sim.Rds")
sim_dt <- lapply(1:length(sim_dt), function(i) {
    cur_dt <- sim_dt[[i]]
    cur_dt$siteID <- i
    return(cur_dt)
})

sim_dt <- do.call(rbind, sim_dt)


FitLmPerSite <- function(st) {
    fit <- lm(frc ~ chil, data = st)
    f <- summary(fit)$fstatistic
    if (is.null(f)) {
        # browser() # For debug
        return(data.table(r2 = 9999, pval = 9999, inter = 9999, slp = 9999))
    }
    p_val <- pf(f[1], f[2], f[3], lower.tail = FALSE)
    r_sqr <- round(summary(fit)$r.squared, 2)

    inter <- coef(fit)[1]
    slp <- coef(fit)[2]

    return(data.table(r2 = r_sqr, pval = p_val, inter = inter, slp = slp))
}

fit_dt <- sim_dt[, FitLmPerSite(.SD), by = siteID]
fit_dt <- fit_dt[pval != 9999, ]

# Significant fit percentage
sig_fit_dt <- fit_dt[pval < 0.05 & slp < 0,]
sig_pct <- round(nrow(sig_fit_dt) / nrow(fit_dt) * 100, 2)

# Average stats from significant fits
avg_inter <- mean(sig_fit_dt$inter, na.rm = TRUE)
avg_slp <- mean(sig_fit_dt$slp, na.rm = TRUE)
avg_r2 <- round(mean(sig_fit_dt$r2, na.rm = TRUE), 2)


pdf("out/fake_chill_forc_hf.pdf", width = 4, height = 4)
# png("out/fake_chill_forc.png", width = 600, height = 600, res = 150)

par(mfrow = c(1, 1), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 2, 2))
smoothScatter(sim_dt[frc > 100, .(chil, frc)], 
    ylim = c(100, 1100), xlim = c(30, 220), 
    pch = 16, cex = 0, 
    colramp = colorRampPalette(
        c("white", rev(hcl.colors(4, "gnbu"))), 
        # c("white"), 
        alpha = TRUE
    ),
    xlab = "Accumulated chilling", ylab = "Accumulated forcing",
    nrpoints = 0
)
# All slope
null <- apply(sig_fit_dt, 1, function(a_row) {
    abline(a = a_row[["inter"]], b = a_row[["slp"]], 
        lwd = 0.3, col = adjustcolor("grey", 0.3))
})
# The average slope
abline(a = avg_inter, b = avg_slp, lwd = 2)

legend("topright", legend = c(paste0("Sig.% = ", sig_pct)), bty = "n")

dev.off()







