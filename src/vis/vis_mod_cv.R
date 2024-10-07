#'******************************************************************************
#' Description: Cross-validation result figure for considered process-based 
#' spring phenology models.
#'******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")

pdf("out/model_cv_comparison.pdf", width = 16, height = 9, pointsize = 13.5)

# png("out/model_cv_comparison.png", 
#     res = 150, width = 2100, height = 1300)

# # Rewrite base color palette
# png("Output/Assets/model_cv_comparison_dark.png",
#     res = 300, width = 4100, height = 2600)
# par(bg = NA, fg = "white", col.axis = "white", col.lab = "white",
#     col.main = "white", col.sub = "white")

layout(matrix(1:28, nrow = 4, byrow = TRUE))
par(mgp = c(1.5, 0.5, 0), oma = c(3, 4, 0, 2))


# ~ MCD12Q2 + Daymet at flux sites ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Pooled ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_mcd_dm_flux_pooled_cv <- lapply(
    list.files(file.path("pipe/mod_cv", "mcd_dm_flux_pool"), full.names = TRUE),
    readRDS
)
names(mod_mcd_dm_flux_pooled_cv) <- sapply(
    list.files(file.path("pipe", "mod_cv", "mcd_dm_flux_pool")), 
    tools::file_path_sans_ext
)

cv_li <- mod_mcd_dm_flux_pooled_cv
cv_dt <- Reformat2Dt(mod_mcd_dm_flux_pooled_cv)

# Conver to site anomalies
cv_dt[, obs_ano := scale(obs, scale = FALSE), by = .(site)]
cv_dt[, TT_ano := scale(TT, scale = FALSE), by = .(site)]
cv_dt[, PA_ano := scale(PA, scale = FALSE), by = .(site)]
cv_dt[, SQ_ano := scale(SQ, scale = FALSE), by = .(site)]
cv_dt[, AT_ano := scale(AT, scale = FALSE), by = .(site)]
cv_dt[, UN_ano := scale(UN, scale = FALSE), by = .(site)]


# Linear regression
lm_cv <- readRDS(file.path("pipe/mod_cv", "mod_mcd_dm_flux_lm_pooled_cv.Rds"))

{ #fig: MCD12Q2 + Daymet at flux sites pooled model comparison
    par(mar = c(0, 0, 2, 0))
    
    ScatterType2XY(lm_cv$mod$obs, lm_cv$mod$pred, 
        ptcol = mcd_dm_color, xlab = "", ylab = "", main = "LIN", 
        xaxt = "n")
    ScatterType2XY(cv_dt$obs_ano, cv_dt$TT_ano, 
        ptcol = mcd_dm_color, xlab = "", ylab = "", main = "TT", 
        xaxt = "n", yaxt = "n")
    ScatterType2XY(cv_dt$obs_ano, cv_dt$PA_ano, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "PA", xaxt = "n", yaxt = "n")
    ScatterType2XY(cv_dt$obs_ano, cv_dt$SQ_ano, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "SQ", xaxt = "n", yaxt = "n")
    ScatterType2XY(cv_dt$obs_ano, cv_dt$AT_ano, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "AT", xaxt = "n", yaxt = "n")
    
    ScatterType2XY(cv_dt[UN_ano != 0, obs_ano], cv_dt[UN_ano != 0, UN_ano], 
        ptcol = mcd_dm_color, xlab = "", ylab = "", main = "UN", 
        xaxt = "n", yaxt = "n"
    )
    points(cv_dt[UN_ano == 0, obs_ano], cv_dt[UN_ano == 0, UN_ano], 
        col = adjustcolor(1, alpha.f = 0.2), 
        xlab = "", ylab = "", main = "UN", 
        xaxt = "n", yaxt = "n"
    )
}


# ~ Site-specific ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_mcd_dm_flux_site_cv <- lapply(
    list.files(file.path("pipe/mod_cv", "mcd_dm_flux_site"), full.names = TRUE),
    readRDS
)
names(mod_mcd_dm_flux_site_cv) <- sapply(
    list.files(file.path("pipe", "mod_cv", "mcd_dm_flux_site")), 
    tools::file_path_sans_ext
)
lm_cv <- readRDS(file.path(
    "pipe/mod_cv",
    "mod_mcd_dm_flux_lm_site_cv.Rds"
))

mod_stats_mcd_dm <- data.table()
for (site in names(mod_mcd_dm_flux_site_cv)) {
    rmse_LIN <- lm_cv[[site]]$rmse
    
    obs <- mod_mcd_dm_flux_site_cv[[site]]$obs
    if (is.null(obs)) {
        next
    }
    rmse_TT <- CalRmse(obs, mod_mcd_dm_flux_site_cv[[site]]$pred_TT)
    rmse_PA <- CalRmse(obs, mod_mcd_dm_flux_site_cv[[site]]$pred_PA)
    rmse_SQ <- CalRmse(obs, mod_mcd_dm_flux_site_cv[[site]]$pred_SQ)
    rmse_AT <- CalRmse(obs, mod_mcd_dm_flux_site_cv[[site]]$pred_AT)
    rmse_UN <- CalRmse(obs, mod_mcd_dm_flux_site_cv[[site]]$pred_UN)

    mod_stats_mcd_dm <- rbind(mod_stats_mcd_dm, data.table(
        siteID = rep(site, 6),
        rmse = c(rmse_LIN, rmse_TT, rmse_PA, rmse_SQ, rmse_AT, rmse_UN),
        model = c("LIN", "TT", "PA", "SQ", "AT", "UN")
    ))
}
mod_stats_mcd_dm$model <- factor(mod_stats_mcd_dm$model,
    levels = c("LIN", "TT", "PA", "SQ", "AT", "UN")
)

par(mar = c(0, 3, 2, 0))
boxplot(rmse ~ model,
    data = mod_stats_mcd_dm[mod_stats_mcd_dm$rmse < 20], 
    outline = FALSE, col = mcd_dm_color,
    xlab = "", ylab = "RMSE", xaxt = "n", yaxt = "n"
)
axis(
    side = 2,
    tck = 0.02,
    at = axTicks(side = 2),
    cex.axis = 1.2
)





# ~ ground pheno and ground temperature at HF & HB ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Pooled ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_hb_hf_pooled_cv <- readRDS(file.path(
    "pipe", "mod_cv",
    "mod_hf_hb_pooled_cv.Rds"
))
cv_li <- mod_hb_hf_pooled_cv

cv_dt <- Reformat2Dt(mod_hb_hf_pooled_cv)

# Conver to site anomalies
cv_dt[, obs_ano := scale(obs, scale = FALSE), by = .(site)]
cv_dt[, TT_ano := scale(TT, scale = FALSE), by = .(site)]
cv_dt[, PA_ano := scale(PA, scale = FALSE), by = .(site)]
cv_dt[, SQ_ano := scale(SQ, scale = FALSE), by = .(site)]
cv_dt[, AT_ano := scale(AT, scale = FALSE), by = .(site)]
cv_dt[, UN_ano := scale(UN, scale = FALSE), by = .(site)]


# Linear regression
lm_cv <- readRDS(file.path(
    "pipe", "mod_cv",
    "mod_hb_hf_lm_pooled_cv.Rds"
))

{ # fig: MCD12Q2 + Daymet at flux sites pooled model comparison
    par(mar = c(0, 0, 0, 0))
    ScatterType2XY(lm_cv$mod$obs, lm_cv$mod$pred,
        ptcol = grd_pheno_temp_color, xlab = "", ylab = "",
        xaxt = "n"
    )
    ScatterType2XY(cv_dt$obs_ano, cv_dt$TT_ano,
        ptcol = grd_pheno_temp_color, xlab = "", ylab = "",
        xaxt = "n", yaxt = "n"
    )
    ScatterType2XY(cv_dt$obs_ano, cv_dt$PA_ano,
        ptcol = grd_pheno_temp_color,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n"
    )
    ScatterType2XY(cv_dt$obs_ano, cv_dt$SQ_ano,
        ptcol = grd_pheno_temp_color,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n"
    )
    ScatterType2XY(cv_dt$obs_ano, cv_dt$AT_ano,
        ptcol = grd_pheno_temp_color,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n"
    )
    ScatterType2XY(cv_dt[UN_ano != 0, obs_ano], cv_dt[UN_ano != 0, UN_ano],
        ptcol = grd_pheno_temp_color, xlab = "", ylab = "",
        xaxt = "n", yaxt = "n"
    )
    points(cv_dt[UN_ano == 0, obs_ano], cv_dt[UN_ano == 0, UN_ano],
        col = adjustcolor(1, alpha.f = 0.2),
        xlab = "", ylab = "",
        xaxt = "n", yaxt = "n"
    )
}



# ~ Site-specific ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summarize caches
mod_hf_hb_site_cv <- lapply(
    list.files(file.path("pipe", "mod_cv", "hf_hb_site"), full.names = TRUE), 
    function(x) {
        li <- readRDS(x)
    }
)
names(mod_hf_hb_site_cv) <- sapply(
    list.files(file.path("pipe", "mod_cv", "hf_hb_site")), 
    tools::file_path_sans_ext
)

# mod_hf_hb_site_cv <- readRDS(file.path(
#     "pipe", "mod_cv",
#     "mod_hf_hb_site_cv.Rds"
# ))
lm_cv <- readRDS(file.path(
    "pipe", "mod_cv",
    "mod_hb_hf_lm_site_cv.Rds"
))

mod_stats_hfhb <- data.table()
for (site in names(mod_hf_hb_site_cv)) {
    rmse_LIN <- lm_cv[[site]]$rmse

    obs <- mod_hf_hb_site_cv[[site]]$obs
    if (is.null(obs)) {
        next
    }
    
    rmse_TT <- CalRmse(obs, mod_hf_hb_site_cv[[site]]$pred_TT)
    rmse_PA <- CalRmse(obs, mod_hf_hb_site_cv[[site]]$pred_PA)
    rmse_SQ <- CalRmse(obs, mod_hf_hb_site_cv[[site]]$pred_SQ)
    rmse_AT <- CalRmse(obs, mod_hf_hb_site_cv[[site]]$pred_AT)
    rmse_UN <- CalRmse(obs, mod_hf_hb_site_cv[[site]]$pred_UN)

    mod_stats_hfhb <- rbind(mod_stats_hfhb, data.table(
        siteID = rep(site, 6),
        rmse = c(rmse_LIN, rmse_TT, rmse_PA, rmse_SQ, rmse_AT, rmse_UN),
        model = c("LIN", "TT", "PA", "SQ", "AT", "UN")
    ))
}
mod_stats_hfhb$model <- factor(mod_stats_hfhb$model,
    levels = c("LIN", "TT", "PA", "SQ", "AT", "UN")
)

par(mar = c(0, 3, 0, 0))
boxplot(rmse ~ model,
    data = mod_stats_hfhb[mod_stats_hfhb$rmse < 7], 
    outline = FALSE, col = grd_pheno_temp_color,
    ylab = "RMSE",  xaxt = "n", xlab = "", yaxt = "n"
)
axis(side = 1, 
    gap.axis = 0.01,
    tck = 0.02, 
    at = 1:6, 
    labels = levels(mod_stats_hfhb$model),
    cex.axis = 1.2
)
axis(side = 2, 
    tck = 0.02, 
    at = axTicks(side = 2),
    cex.axis = 1.2
)



# ~ MCD12Q2 + Flux_temperature ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Pooled ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_mcd_flux_temp_pooled_cv <- lapply(
    list.files(
        file.path("pipe/mod_cv", "mcd_flux_temp_pool"), 
        full.names = TRUE
    ),
    readRDS
)
names(mod_mcd_flux_temp_pooled_cv) <- sapply(
    list.files(file.path("pipe", "mod_cv", "mcd_flux_temp_pool")), 
    tools::file_path_sans_ext
)

cv_li <- mod_mcd_flux_temp_pooled_cv

cv_dt <- Reformat2Dt(mod_mcd_flux_temp_pooled_cv)

# Conver to site anomalies
cv_dt[, obs_ano := scale(obs, scale = FALSE), by = .(site)]
cv_dt[, TT_ano := scale(TT, scale = FALSE), by = .(site)]
cv_dt[, PA_ano := scale(PA, scale = FALSE), by = .(site)]
cv_dt[, SQ_ano := scale(SQ, scale = FALSE), by = .(site)]
cv_dt[, AT_ano := scale(AT, scale = FALSE), by = .(site)]
cv_dt[, UN_ano := scale(UN, scale = FALSE), by = .(site)]


# Linear regression
lm_cv <- readRDS(file.path(
    "pipe/mod_cv",
    "mod_mcd_flux_temp_lm_pooled_cv.Rds"
))

{ # fig: MCD12Q2 + Daymet at flux sites pooled model comparison
    par(mar = c(3, 0, 0, 0))
    ScatterType2XY(lm_cv$mod$obs, lm_cv$mod$pred,
        ptcol = mcd_grd_temp_color, xlab = "", ylab = ""
    )
    ScatterType2XY(cv_dt$obs_ano, cv_dt$TT_ano,
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "",
        yaxt = "n"
    )
    ScatterType2XY(cv_dt$obs_ano, cv_dt$PA_ano,
        ptcol = mcd_grd_temp_color,
        xlab = "", ylab = "", yaxt = "n"
    )
    ScatterType2XY(cv_dt$obs_ano, cv_dt$SQ_ano,
        ptcol = mcd_grd_temp_color,
        xlab = "", ylab = "", yaxt = "n"
    )
    ScatterType2XY(cv_dt$obs_ano, cv_dt$AT_ano,
        ptcol = mcd_grd_temp_color,
        xlab = "", ylab = "", yaxt = "n"
    )
    ScatterType2XY(cv_dt[UN_ano != 0, obs_ano], cv_dt[UN_ano != 0, UN_ano],
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "",
        yaxt = "n"
    )
    points(cv_dt[UN_ano == 0, obs_ano], cv_dt[UN_ano == 0, UN_ano],
        col = adjustcolor(1, alpha.f = 0.2),
        xlab = "", ylab = "",
        xaxt = "n", yaxt = "n"
    )
}



# This empty plot is for occupying the space
plot(NA,
    xlim = c(0, 1), ylim = c(0, 1),
    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = ""
)


# ~ PEP725 + E-OBS ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~ Site-specific ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Pep725SiteModels <- function(species_name) {
    caches <- list.files(
        file.path("pipe", "mod_cv", "pep"),
        pattern = species_name,
        full.names = TRUE
    )
    mod_site <- lapply(caches, readRDS)
    names(mod_site) <- lapply(caches, function(x) {
        strsplit(basename(x), "\\.")[[1]][1]
    })
    
    mod_lm_site <- readRDS(file.path(
        "pipe", "mod_cv",
        paste0("mod_pep725_lm_site_", species_name, "_cv.Rds")
    ))

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
        rmse_TT <- CalRmse(obs, mod_site[[site]]$pred_TT)
        rmse_PA <- CalRmse(obs, mod_site[[site]]$pred_PA)
        rmse_SQ <- CalRmse(obs, mod_site[[site]]$pred_SQ)
        rmse_AT <- CalRmse(obs, mod_site[[site]]$pred_AT)
        rmse_UN <- CalRmse(obs, mod_site[[site]]$pred_UN)
        if (rmse_UN > 999 | is.na(rmse_UN)) {
            rmse_UN <- NA
        }
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
    xlab = "", ylab = "RMSE", ylim = c(0, 30), xaxt = "n", yaxt = "n"
)
axis(side = 1,
    gap.axis = 0.01,
    tck = 0.02,
    at = 1:6,
    labels = levels(mod_stat_alnus_site$model),
    cex.axis = 1.2
)
axis(side = 2,
    tck = 0.02,
    at = axTicks(side = 2),
    cex.axis = 1.2,
    gap.axis = 0.5
)
legend("topleft", legend = "Alnus", bty = "n", 
    cex = 1.5, 
    x.intersp = 0, y.intersp = 0.1
)


# ~ aesculus
mod_stat_aesculus_site <- Pep725SiteModels("aesculus")

# fit <- lm(rmse ~ factor(model), data = mod_stat_aesculus_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
boxplot(rmse ~ model,
    data = mod_stat_aesculus_site, outline = FALSE, col = pep_color,
    xlab = "", yaxt = "n", ylab = "", ylim = c(0, 30), xaxt = "n"
)
axis(side = 1,
    gap.axis = 0.01,
    tck = 0.02,
    at = 1:6,
    labels = levels(mod_stat_aesculus_site$model),
    cex.axis = 1.2
)
legend("topleft", legend = "Aesculus", bty = "n", 
    cex = 1.5, 
    x.intersp = 0, y.intersp = 0.1
)


# ~ betula
mod_stat_betula_site <- Pep725SiteModels("betula")

# fit <- lm(rmse ~ factor(model), data = mod_stat_betula_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
boxplot(rmse ~ model,
    data = mod_stat_betula_site, outline = FALSE, col = pep_color,
    xlab = "", ylab = "", yaxt = "n", ylim = c(0, 25), xaxt = "n"
)
axis(side = 1,
    gap.axis = 0.01,
    tck = 0.02,
    at = 1:6,
    labels = levels(mod_stat_betula_site$model),
    cex.axis = 1.2
)
legend("topleft", legend = "Betula", bty = "n",
    cex = 1.5, 
    x.intersp = 0, y.intersp = 0.1
)


# ~ fagus
mod_stat_fagus_site <- Pep725SiteModels("fagus")

# fit <- lm(rmse ~ factor(model), data = mod_stat_fagus_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
boxplot(rmse ~ model,
    data = mod_stat_fagus_site, outline = FALSE, col = pep_color,
    xlab = "", ylab = "", yaxt = "n", ylim = c(0, 25), xaxt = "n"
)
axis(side = 1,
    gap.axis = 0.01,
    tck = 0.02,
    at = 1:6,
    labels = levels(mod_stat_fagus_site$model),
    cex.axis = 1.2
)
legend("topleft", legend = "Fagus", bty = "n",
    cex = 1.5, 
    x.intersp = 0, y.intersp = 0.1
)


# ~ fraxinus
mod_stat_fraxinus_site <- Pep725SiteModels("fraxinus")

# fit <- lm(rmse ~ factor(model), data = mod_stat_fraxinus_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
boxplot(rmse ~ model,
    data = mod_stat_fraxinus_site, outline = FALSE, col = pep_color,
    xlab = "", ylab = "", yaxt = "n", ylim = c(0, 25), xaxt = "n"
)
axis(side = 1,
    gap.axis = 0.01,
    tck = 0.02,
    at = 1:6,
    labels = levels(mod_stat_fraxinus_site$model),
    cex.axis = 1.2
)
legend("topleft", legend = "Fraxinus", bty = "n",
    cex = 1.5, 
    x.intersp = 0, y.intersp = 0.1
)


# ~ quercus
mod_stat_quercus_site <- Pep725SiteModels("quercus")

# fit <- lm(rmse ~ factor(model), data = mod_stat_quercus_site)
# summary(fit)

# fig: grd pheno and temperature at HF & HB for site-specific model comparison
boxplot(rmse ~ model,
    data = mod_stat_quercus_site, outline = FALSE, col = pep_color,
    xlab = "", ylab = "", yaxt = "n", ylim = c(0, 25), xaxt = "n"
)
axis(side = 1,
    gap.axis = 0.01,
    tck = 0.02,
    at = 1:6,
    labels = levels(mod_stat_quercus_site$model),
    cex.axis = 1.2
)
legend("topleft", legend = "Quercus", bty = "n",
    cex = 1.5, 
    x.intersp = 0, y.intersp = 0.1
)





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
    cex = 1.5, pt.cex = 2, xpd = NA
)

# legend(grconvertX(0.01, "ndc"), grconvertY(0.05, "ndc"), 
#     bty = "n", 
#     legend = c(
#         "MCD12Q2 + Daymet", 
#         "Ground phenology + Ground temperature",
#         "MCD12Q2 + Ground temperature",
#         "PEP725 + E-OBS"
#     ), 
#     pch = 16, 
#     ncol = 4,
#     col = c(mcd_dm_color, grd_pheno_temp_color, mcd_grd_temp_color, 
#         pep_color), 
#     x.intersp = 0.5,
#     cex = 1.8, pt.cex = 2, xpd = NA
# )

# text(grconvertX(0.01, "ndc"), grconvertY(0.99, "ndc"),
#     label = "a", cex = 2, xpd = NA
# )
# text(grconvertX(0.01, "ndc"), grconvertY(0.23, "ndc"),
#     label = "b", cex = 2, xpd = NA
# )



dev.off()
