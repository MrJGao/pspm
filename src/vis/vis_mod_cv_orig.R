# ******************************************************************************
# The Cross-validation figure with original SOS dates instead of anomalies
# ******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")


DrawRow <- function(pooled_cv, lm_cv_file, color, 
    ifmain = FALSE, ifxaxt = FALSE
) {
    # pooled_cv <- readRDS(pooled_cv_file)
    cv_li <- pooled_cv
    cv_dt <- Reformat2Dt(pooled_cv)

    # Linear regression
    lm_cv <- readRDS(lm_cv_file)

    rg <- range(
        lm_cv$pred,
        cv_dt$obs, cv_dt$TT, cv_dt$PA, cv_dt$SQ, cv_dt$AT, cv_dt[UN > 0]$UN,
        na.rm = TRUE
    )

    ScatterType2XY(lm_cv$mod$obs, lm_cv$mod$pred,
        ptcol = color, xlab = "", ylab = "",
        range = rg,
        xaxt = "n"
    )
    if (ifmain) {
        title("LIN")
    }
    if (ifxaxt) {
        axis(side = 1)
    }

    ScatterType2XY(cv_dt$obs, cv_dt$TT,
        ptcol = color, xlab = "", ylab = "",
        xaxt = "n", yaxt = "n",
        range = rg
    )
    if (ifmain) {
        title("TT")
    }
    if (ifxaxt) {
        axis(side = 1)
    }

    ScatterType2XY(cv_dt$obs, cv_dt$PA,
        ptcol = color,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n",
        range = rg
    )
    if (ifmain) {
        title("PA")
    }
    if (ifxaxt) {
        axis(side = 1)
    }

    ScatterType2XY(cv_dt$obs, cv_dt$SQ,
        ptcol = color,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n",
        range = rg
    )
    if (ifmain) {
        title("SQ")
    }
    if (ifxaxt) {
        axis(side = 1)
    }
    ScatterType2XY(cv_dt$obs, cv_dt$AT,
        ptcol = color,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n",
        range = rg
    )
    if (ifmain) {
        title("AT")
    }
    if (ifxaxt) {
        axis(side = 1)
    }

    ScatterType2XY(cv_dt[UN >= 0, obs], cv_dt[UN >= 0, UN],
        ptcol = color, xlab = "", ylab = "",
        xaxt = "n", yaxt = "n",
        range = rg
    )
    if (ifmain) {
        title("UN")
    }
    if (ifxaxt) {
        axis(side = 1)
    }
    points(cv_dt[UN < 0, obs], cv_dt[UN < 0, UN],
        col = adjustcolor(1, alpha.f = 0.2),
        xlab = "", ylab = "",
        xaxt = "n", yaxt = "n"
    )

}


png("out/model_cv_comparison_orig.png",
    res = 150, width = 1600, height = 1000
)
layout(matrix(1:18, nrow = 3, byrow = TRUE))
par(mgp = c(1.5, 0.5, 0), oma = c(6, 4, 0, 2))


# ~ MCD12Q2 + Daymet at flux sites ####
mod_mcd_dm_flux_pooled_cv <- lapply(
    list.files(file.path("pipe/mod_cv", "mcd_dm_flux_pool"), full.names = TRUE),
    readRDS
)
names(mod_mcd_dm_flux_pooled_cv) <- sapply(
    list.files(file.path("pipe", "mod_cv", "mcd_dm_flux_pool")), 
    tools::file_path_sans_ext
)

par(mar = c(0, 0, 2, 0))
DrawRow(
    pooled_cv = mod_mcd_dm_flux_pooled_cv,
    lm_cv_file = file.path(
        "pipe/mod_cv",
        "mod_mcd_dm_flux_lm_pooled_cv_orig.Rds"
    ),
    color = mcd_dm_color,
    ifmain = TRUE
)


# ~ ground pheno and ground temperature at HF & HB ####
mod_hb_hf_pooled_cv <- readRDS(file.path(
    "pipe", "mod_cv",
    "mod_hf_hb_pooled_cv.Rds"
))
par(mar = c(0, 0, 0, 0))
DrawRow(
    pooled_cv = mod_hb_hf_pooled_cv,
    lm_cv_file = file.path(
        "pipe/mod_cv",
        "mod_hb_hf_lm_pooled_cv_orig.Rds"
    ),
    color = grd_pheno_temp_color,
    ifmain = FALSE
)

# ~ MCD12Q2 + Flux_temperature ####
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

par(mar = c(0, 0, 0, 0))
DrawRow(
    pooled_cv = mod_mcd_flux_temp_pooled_cv,
    lm_cv_file = file.path(
        "pipe/mod_cv",
        "mod_mcd_flux_temp_lm_pooled_cv_orig.Rds"
    ),
    color = mcd_grd_temp_color,
    ifmain = FALSE,
    ifxaxt = TRUE
)


# ~ Label annotations ####
# ~ ----------------------------------------------------------------------------
text(grconvertX(0.02, "ndc"), grconvertY(0.5, "ndc"),
    label = "Est. SOS",
    srt = 90, xpd = NA, cex = 1.5
)
text(grconvertX(0.5, "ndc"), grconvertY(0.07, "ndc"),
    label = "Obs. SOS",
    xpd = NA, cex = 1.5
)

legend(grconvertX(0.5, "ndc"), grconvertY(0.07, "ndc"),
    bty = "n",
    legend = c(
        "MCD12Q2 + Daymet",
        "HF & HB",
        "MCD12Q2 + Flux temperature"
    ),
    pch = 21,
    ncol = 3,
    pt.bg = c(
        mcd_dm_color,
        grd_pheno_temp_color,
        mcd_grd_temp_color
    ),
    xjust = 0.5,
    cex = 1.5, pt.cex = 2,
    xpd = NA
)


dev.off()
