#'******************************************************************************
#' Description: Visualize model fit for multiple data sources. The input is 
#' generated from `mod_fit_multi_data.R`
#'******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")



# png("out/model_fit_comparison.png", 
#     res = 150, width = 2000, height = 1300
# )

pdf("out/model_fit_comparison.pdf", width = 16, height = 9)

# png("Output/Fig/model_fit_comparison_dark.png", 
#     res = 300, width = 3600, height = 2600)
# par(bg = NA, fg = "white", col.axis = "white", col.lab = "white", 
#     col.main = "white", col.sub = "white")
# Rewrite base color palette


layout(matrix(1:28, nrow = 4, byrow = TRUE))
par(mgp = c(1.5, 0.5, 0), oma = c(3, 4, 0, 2))


# ~ MCD12Q2 + Daymet at flux sites ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Pooled ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_mcd_dm_flux_pooled <- readRDS(
    file.path("pipe/goodness-of-fit", "mcd_dm_flux_pool",
        "mod_mcd_dm_flux_pooled.Rds"
    )
)
flux_li <- readRDS(file.path("pipe", "flux_li.Rds"))

# Compute site anomalies
anomaly_dt <- CalSiteAnomaly(mod_mcd_dm_flux_pooled, flux_li)

# Linear regression
mod_mcd_dm_flux_lm_pooled <- readRDS(file.path("pipe", "goodness-of-fit",
    "mod_mcd_dm_flux_lm_pooled.Rds")
)
mod_mcd_dm_flux_lm_pooled$obs <- fitted(mod_mcd_dm_flux_lm_pooled$mod) +
    residuals(mod_mcd_dm_flux_lm_pooled$mod)


{ #fig: MCD12Q2 + Daymet at flux sites pooled model comparison
    par(mar = c(0, 0, 2, 0))
    ScatterType2XY(mod_mcd_dm_flux_lm_pooled$obs, 
        fitted(mod_mcd_dm_flux_lm_pooled$mod),
        ptcol = mcd_dm_color, 
        xlab = "", ylab = "", xaxt = "n",
        main = "LIN"
    )
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$TT_anom, 
        ptcol = mcd_dm_color, xlab = "", ylab = "", main = "TT", 
        xaxt = "n", yaxt = "n"
    )
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$PA_anom, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "PA", xaxt = "n", yaxt = "n")
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$SQ_anom, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "SQ", xaxt = "n", yaxt = "n")
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$AT_anom, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "AT", xaxt = "n", yaxt = "n")
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$UN_anom, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "UN", xaxt = "n", yaxt = "n")
}


# ~ Site-specific ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_mcd_dm_flux_site <- lapply(
    list.files(
        file.path("pipe/goodness-of-fit", "mcd_dm_flux_site"), 
        full.names = TRUE
    ),
    readRDS
)
names(mod_mcd_dm_flux_site) <- sapply(
    list.files(file.path("pipe/goodness-of-fit", "mcd_dm_flux_site")), 
    tools::file_path_sans_ext
)


# Linear regression
mod_mcd_dm_flux_lm_site <- readRDS(
    file.path("pipe", "goodness-of-fit", "mod_mcd_dm_flux_lm_site.Rds")
)


mod_stats_mcd_dm <- data.table()
for(site in names(mod_mcd_dm_flux_site)) {
    obs <- mod_mcd_dm_flux_site[[site]]$obs
    if (is.null(obs)) { next }

    rmse_LIN <- mod_mcd_dm_flux_lm_site[[site]]$rmse
    rmse_TT <- CalRmse(obs, mod_mcd_dm_flux_site[[site]]$TT$fitted_doy)
    rmse_PA <- CalRmse(obs, mod_mcd_dm_flux_site[[site]]$PA$fitted_doy)
    rmse_SQ <- CalRmse(obs, mod_mcd_dm_flux_site[[site]]$SQ$fitted_doy)
    rmse_AT <- CalRmse(obs, mod_mcd_dm_flux_site[[site]]$AT$fitted_doy)
    rmse_UN <- CalRmse(obs, mod_mcd_dm_flux_site[[site]]$UN$fitted_doy)

    mod_stats_mcd_dm <- rbind(mod_stats_mcd_dm, data.table(
        siteID = rep(site, 6), 
        rmse = c(rmse_LIN, rmse_TT, rmse_PA, rmse_SQ, rmse_AT, rmse_UN), 
        model = c("LIN", "TT", "PA", "SQ", "AT", "UN")
    ))
}
mod_stats_mcd_dm$model <- factor(mod_stats_mcd_dm$model, 
    levels = c("LIN", "TT", "PA", "SQ", "AT", "UN")
)


# fit <- lm(rmse ~ factor(model), data = mod_stats_mcd_dm)
# summary(fit)
# fig: MCD12Q2 + Daymet at flux sites site-specific model comparison
par(mar = c(0, 3, 2, 0))
boxplot(rmse ~ model,
    data = mod_stats_mcd_dm, outline = FALSE, col = mcd_dm_color,
    xlab = "", ylab = "RMSE", xaxt = "n"
)





# ~ ground pheno and ground temperature at HF & HB ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Pooled ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_hb_hf_pooled <- readRDS(file.path(
    "pipe/goodness-of-fit", "hf_hb_pool",
    "mod_hb_hf_pooled.Rds"
))
hb_hf_li <- readRDS("pipe/hb_hf_li.Rds")

anomaly_dt <- CalSiteAnomaly(mod_hb_hf_pooled, hb_hf_li)

# Linear regression
mod_hb_hf_lm_pooled <- readRDS(
    file.path("pipe/goodness-of-fit", "mod_hb_hf_lm_pooled.Rds")
)
mod_hb_hf_lm_pooled$obs <- fitted(mod_hb_hf_lm_pooled$mod) +
    residuals(mod_hb_hf_lm_pooled$mod)

{ # fig: ground pheno and temperature at HF & HB for pooled model comparison
    par(mar = c(0, 0, 0, 0))
    ScatterType2XY(mod_hb_hf_lm_pooled$obs, fitted(mod_hb_hf_lm_pooled$mod),
        ptcol = grd_pheno_temp_color, xlab = "", ylab = "", xaxt = "n")
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$TT_anom, 
        ptcol = grd_pheno_temp_color, 
        xlab = "", ylab = "", xaxt = "n", yaxt = "n"
    )
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$PA_anom, 
        ptcol = grd_pheno_temp_color, 
        xlab = "", ylab = "", xaxt = "n", yaxt = "n"
    )
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$SQ_anom, 
        ptcol = grd_pheno_temp_color, 
        xlab = "", ylab = "", xaxt = "n", yaxt = "n"
    )
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$AT_anom,
        ptcol = grd_pheno_temp_color, 
        xlab = "", ylab = "", xaxt = "n", yaxt = "n"
    )
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$UN_anom, 
        ptcol = grd_pheno_temp_color, 
        xlab = "", ylab = "", xaxt = "n", yaxt = "n"
    )
    # dev.off()
}


# ~ Site-specific ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_hb_hf_site <- readRDS("pipe/goodness-of-fit/hf_hb_site/mod_hb_hf_site.Rds")

# Linear regression
mod_hb_hf_lm_site <- readRDS(file.path(
    "pipe", "goodness-of-fit",
    "mod_hb_hf_lm_site.Rds"
))


mod_stat_hb_hf_site <- data.table()
for (site in names(mod_hb_hf_site)) {
    obs <- mod_hb_hf_site[[site]]$obs
    if (is.null(obs)) {
        next
    }
    rmse_LIN <- mod_hb_hf_lm_site[[site]]$rmse
    rmse_TT <- CalRmse(obs, mod_hb_hf_site[[site]]$TT$fitted_doy)
    rmse_PA <- CalRmse(obs, mod_hb_hf_site[[site]]$PA$fitted_doy)
    rmse_SQ <- CalRmse(obs, mod_hb_hf_site[[site]]$SQ$fitted_doy)
    rmse_AT <- CalRmse(obs, mod_hb_hf_site[[site]]$AT$fitted_doy)
    rmse_UN <- CalRmse(obs, mod_hb_hf_site[[site]]$UN$fitted_doy)
    mod_stat_hb_hf_site <- rbind(mod_stat_hb_hf_site, data.table(
        siteID = rep(site, 6), 
        rmse = c(rmse_LIN, rmse_TT, rmse_PA, rmse_SQ, rmse_AT, rmse_UN),
        model = c("LIN", "TT", "PA", "SQ", "AT", "UN")
    ))
}
mod_stat_hb_hf_site$model <- factor(mod_stat_hb_hf_site$model,
    levels = c("LIN", "TT", "PA", "SQ", "AT", "UN")
)

# fit <- lm(rmse ~ factor(model), data = mod_stat_hb_hf_site)
# summary(fit)
# fig: ground pheno and temperature at HF & HB for site-specific model comparison
par(mar = c(0, 3, 0, 0))
boxplot(rmse ~ model,
    data = mod_stat_hb_hf_site, outline = FALSE, col = grd_pheno_temp_color,
    xlab = "Models", ylab = "RMSE", gap.axis = 0.01
)



# ~ MCD12Q2 + Flux_temperature ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_mcd_flux_temp_pooled <- readRDS(
    file.path("pipe", "goodness-of-fit", "mcd_flux_temp_pool",
        "mod_mcd_flux_temp_pooled.Rds"
    )
)
flux_temp_li <- readRDS(file.path("pipe", "flux_temp_li.Rds"))

anomaly_dt <- CalSiteAnomaly(mod_mcd_flux_temp_pooled, flux_temp_li)

# Linear regression
mod_mcd_flux_temp_lm_pooled <- readRDS(file.path("pipe", "goodness-of-fit",
    "mod_mcd_flux_temp_lm_pooled.Rds"
))
mod_mcd_flux_temp_lm_pooled$obs <- fitted(mod_mcd_flux_temp_lm_pooled$mod) +
    residuals(mod_mcd_flux_temp_lm_pooled$mod)


{ #fig: MCD12Q2 + flux_temperature for pooled model comparison
    par(mar = c(3, 0, 0, 0))
    ScatterType2XY(mod_mcd_flux_temp_lm_pooled$obs, 
        fitted(mod_mcd_flux_temp_lm_pooled$mod),
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "")
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$TT_anom, 
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "", yaxt = "n")
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$PA_anom, 
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "", yaxt = "n")
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$SQ_anom, 
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "", yaxt = "n")
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$AT_anom,
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "", yaxt = "n")
    ScatterType2XY(anomaly_dt$obs_anom, anomaly_dt$UN_anom, 
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "", yaxt = "n")
    # dev.off()
}


# This empty plot is for occupying the space
plot(NA, xlim = c(0, 1), ylim = c(0, 1), 
    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")




# ~ PEP725 + E-OBS ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Site-specific ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Pep725SiteModels <- function(species_name) {
    mod_site <- readRDS(file.path("pipe", "goodness-of-fit",
        paste0("mod_", species_name, "_pep725_eobs.Rds")))
    mod_lm_site <- readRDS(file.path("pipe", "goodness-of-fit", 
        paste0("mod_pep725_lm_site_", species_name, ".Rds")))

    mod_stat_site <- data.table()
    for (site in names(mod_site)) {
        # extract rmse for linear regression
        lm_site <- strsplit(site, "_")[[1]][2]
        rmse_LIN <- mod_lm_site[[lm_site]]$rmse
        
        # extract rmse for others
        obs <- mod_site[[site]]$obs
        if (is.null(obs)) {
            next
        }
        rmse_TT <- CalRmse(obs, mod_site[[site]]$TT$fitted_doy)
        rmse_PA <- CalRmse(obs, mod_site[[site]]$PA$fitted_doy)
        rmse_SQ <- CalRmse(obs, mod_site[[site]]$SQ$fitted_doy)
        rmse_AT <- CalRmse(obs, mod_site[[site]]$AT$fitted_doy)
        rmse_UN <- CalRmse(obs, mod_site[[site]]$UN$fitted_doy)
        if (rmse_UN > 999) rmse_UN <- NA
        mod_stat_site <- rbind(mod_stat_site, data.table(
            siteID = rep(site, 6), 
            rmse = c(rmse_LIN, rmse_TT, rmse_PA, rmse_SQ, rmse_AT, rmse_UN),
            model = c("LIN", "TT", "PA", "SQ", "AT", "UN")
        ))
    }
    mod_stat_site$model <- factor(mod_stat_site$model,
        levels = c("LIN", "TT", "PA", "SQ", "AT", "UN")
    )

    return(mod_stat_site)
}


# ~ alnus
mod_stat_alnus_site <- Pep725SiteModels("alnus")

# fit <- lm(rmse ~ factor(model), data = mod_stat_alnus_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
par(mar = c(4, 0, 2, 0))
boxplot(rmse ~ model,
    data = mod_stat_alnus_site, outline = FALSE, col = pep_color,
    xlab = "", ylab = "RMSE", ylim = c(0, 25)
)
legend("topleft", legend = "Alnus", bty = "n")

# ~ aesculus
mod_stat_aesculus_site <- Pep725SiteModels("aesculus")

# fit <- lm(rmse ~ factor(model), data = mod_stat_aesculus_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
boxplot(rmse ~ model,
    data = mod_stat_aesculus_site, outline = FALSE, col = pep_color,
    xlab = "", yaxt = "n", ylab = "", ylim = c(0, 25)
)
legend("topleft", legend = "Aesculus", bty = "n")


# ~ betula
mod_stat_betula_site <- Pep725SiteModels("betula")

# fit <- lm(rmse ~ factor(model), data = mod_stat_betula_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
boxplot(rmse ~ model,
    data = mod_stat_betula_site, outline = FALSE, col = pep_color,
    xlab = "", ylab = "", yaxt = "n", ylim = c(0, 25)
)
legend("topleft", legend = "Betula", bty = "n")


# ~ fagus
mod_stat_fagus_site <- Pep725SiteModels("fagus")

# fit <- lm(rmse ~ factor(model), data = mod_stat_fagus_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
boxplot(rmse ~ model,
    data = mod_stat_fagus_site, outline = FALSE, col = pep_color,
    xlab = "", ylab = "", yaxt = "n", ylim = c(0, 25)
)
legend("topleft", legend = "Fagus", bty = "n")


# ~ fraxinus
mod_stat_fraxinus_site <- Pep725SiteModels("fraxinus")

# fit <- lm(rmse ~ factor(model), data = mod_stat_fraxinus_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
boxplot(rmse ~ model,
    data = mod_stat_fraxinus_site, outline = FALSE, col = pep_color,
    xlab = "", ylab = "", yaxt = "n", ylim = c(0, 25)
)
legend("topleft", legend = "Fraxinus", bty = "n")


# ~ quercus
mod_stat_quercus_site <- Pep725SiteModels("quercus")

# fit <- lm(rmse ~ factor(model), data = mod_stat_quercus_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
boxplot(rmse ~ model,
    data = mod_stat_quercus_site, outline = FALSE, col = pep_color,
    xlab = "", ylab = "", yaxt = "n", ylim = c(0, 25)
)
legend("topleft", legend = "Quercus", bty = "n")




# ~ Lable annotations ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
text(grconvertX(0.015, "ndc"), grconvertY(0.63, "ndc"),
    label = "Est. SOS Anomaly",
    srt = 90, xpd = NA, cex = 1.5
)
text(grconvertX(0.44, "ndc"), grconvertY(0.29, "ndc"),
    label = "Obs. SOS Anomaly",
    xpd = NA, cex = 1.5
)

text(grconvertX(0.015, "ndc"), grconvertY(0.18, "ndc"),
    label = "RMSE", srt = 90, xpd = NA, cex = 1.5
)
text(grconvertX(0.44, "ndc"), grconvertY(0.07, "ndc"),
    label = "Models", cex = 1.5, xpd = NA
)

legend(grconvertX(0.85, "ndc"), grconvertY(0.3, "ndc"),
    bty = "n",
    legend = c(
        "MCD12Q2 + Daymet",
        "HF & HB",
        "MCD12Q2 + \n Flux temperature",
        "PEP725 + E-OBS"
    ),
    pch = 21,
    ncol = 1,
    pt.bg = c(
        mcd_dm_color, 
        grd_pheno_temp_color, 
        mcd_grd_temp_color,
        pep_color
    ),
    y.intersp = 1.5,
    cex = 1.5, pt.cex = 2, 
    xpd = NA
)

# legend(grconvertX(0.07, "ndc"), grconvertY(0.05, "ndc"), bty = "n", 
#     legend = c(
#         "MCD12Q2 + Daymet", 
#         "Ground phenology + Ground temperature",
#         "MCD12Q2 + Ground temperature",
#         "PEP725 + E-OBS"
#     ), pch = 16, ncol = 4,
#     col = c(mcd_dm_color, grd_pheno_temp_color, mcd_grd_temp_color, 
#         pep_color), 
#     cex = 1.3, pt.cex = 1.5, xpd = NA)

# text(grconvertX(0.01, "ndc"), grconvertY(0.99, "ndc"),
#     label = "a", cex = 2, xpd = NA
# )
# text(grconvertX(0.01, "ndc"), grconvertY(0.23, "ndc"),
#     label = "b", cex = 2, xpd = NA
# )



dev.off()
