#'******************************************************************************
#' Description: Compare spring phenology model fits by Daymet and flux 
#' temperature.
#'******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")


png("out/model_compare_mcd_dm_flux_temp.png", 
    res = 150, width = 1500, height = 700)

layout(matrix(1:12, nrow = 2, byrow = TRUE))
par(mgp = c(1.5, 0.5, 0), oma = c(3, 4, 0, 2))



# ~ MCD12Q2 + Daymet at flux sites that also have flux temperature ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_mcd_dm_flux_pooled <- readRDS(file.path(gdir, "Pipeline", 
    "mod_mcd_dm_flux-temp-sites_pooled.Rds"))
# Read the saved Rds
flux_dm_temp_li <- readRDS(file.path(gdir, "Pipeline", "flux_dm_temp_li.Rds"))

# Compute site anomalies
anomaly_dt <- CalSiteAnomaly(mod_mcd_dm_flux_pooled, flux_dm_temp_li)

# Linear regression
mod_mcd_dm_flux_lm_pooled <- readRDS(file.path(gdir, "Pipeline",
    "mod_mcd_dm_flux_temp_pooled.Rds"))
mod_mcd_dm_flux_lm_pooled$obs <- fitted(mod_mcd_dm_flux_lm_pooled$mod) +
    residuals(mod_mcd_dm_flux_lm_pooled$mod)


{ #fig: MCD12Q2 + Daymet at flux sites pooled model comparison
    par(mar = c(0, 0, 2, 0))
    ScatterXY(mod_mcd_dm_flux_lm_pooled$obs, 
        fitted(mod_mcd_dm_flux_lm_pooled$mod),
        ptcol = mcd_dm_color, xlab = "", ylab = "", main = "LIN", xaxt = "n")
    ScatterXY(anomaly_dt$obs_anom, anomaly_dt$TT_anom, 
        ptcol = mcd_dm_color, xlab = "", ylab = "", main = "TT", 
        xaxt = "n", yaxt = "n")
    ScatterXY(anomaly_dt$obs_anom, anomaly_dt$PA_anom, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "PA", xaxt = "n", yaxt = "n")
    ScatterXY(anomaly_dt$obs_anom, anomaly_dt$SQ_anom, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "SQ", xaxt = "n", yaxt = "n")
    ScatterXY(anomaly_dt$obs_anom, anomaly_dt$AT_anom, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "AT", xaxt = "n", yaxt = "n")
    ScatterXY(anomaly_dt$obs_anom, anomaly_dt$UN_anom, ptcol = mcd_dm_color, 
        xlab = "", ylab = "", main = "UN", xaxt = "n", yaxt = "n")
}


# ~ MCD12Q2 + Flux_temperature ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_mcd_flux_temp_pooled <- readRDS(file.path(gdir, "Pipeline", 
    "mod_mcd_flux_temp_pooled.Rds"))
flux_temp_li <- readRDS(file.path(gdir, "Pipeline", "flux_temp_li.Rds"))

anomaly_dt <- CalSiteAnomaly(mod_mcd_flux_temp_pooled, flux_temp_li)

# Linear regression
mod_mcd_flux_temp_lm_pooled <- readRDS(file.path(gdir, "Pipeline",
    "mod_mcd_flux_temp_lm_pooled.Rds"))
mod_mcd_flux_temp_lm_pooled$obs <- fitted(mod_mcd_flux_temp_lm_pooled$mod) +
    residuals(mod_mcd_flux_temp_lm_pooled$mod)


{ #fig: MCD12Q2 + flux_temperature for pooled model comparison
    par(mar = c(3, 0, 0, 0))
    ScatterXY(mod_mcd_flux_temp_lm_pooled$obs, 
        fitted(mod_mcd_flux_temp_lm_pooled$mod),
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "")
    ScatterXY(anomaly_dt$obs_anom, anomaly_dt$TT_anom, 
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "", yaxt = "n")
    ScatterXY(anomaly_dt$obs_anom, anomaly_dt$PA_anom, 
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "", yaxt = "n")
    ScatterXY(anomaly_dt$obs_anom, anomaly_dt$SQ_anom, 
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "", yaxt = "n")
    ScatterXY(anomaly_dt$obs_anom, anomaly_dt$AT_anom,
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "", yaxt = "n")
    ScatterXY(anomaly_dt$obs_anom, anomaly_dt$UN_anom, 
        ptcol = mcd_grd_temp_color, xlab = "", ylab = "", yaxt = "n")
    # dev.off()
}



# ~ Lable annotations ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
text(grconvertX(0.02, "ndc"), grconvertY(0.5, "ndc"),
    label = "Est. SOS Anomaly",
    srt = 90, xpd = NA, cex = 1.5
)
text(grconvertX(0.5, "ndc"), grconvertY(0.1, "ndc"),
    label = "Obs. SOS Anomaly",
    xpd = NA, cex = 1.5
)

legend(grconvertX(0.3, "ndc"), grconvertY(0.1, "ndc"),
    bty = "n",
    legend = c(
        "MCD12Q2 + Daymet",
        "MCD12Q2 + Ground temperature"
    ), pch = 16, ncol = 4,
    col = c(
        mcd_dm_color, mcd_grd_temp_color
    ),
    cex = 1.3, pt.cex = 1.5, xpd = NA
)


dev.off()

