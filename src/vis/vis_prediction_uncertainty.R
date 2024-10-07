#'******************************************************************************
#' Description: Prediction uncertainties for the process-based spring phenology
#' models.
#'******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")


library(data.table)

# Example data
hb_hf_li <- readRDS(file.path(hpc_local_dir, "Pipeline", "hb_hf_li.Rds"))
# Bootstrapping results
params <- readRDS(file.path(hpc_local_dir, "Pipeline", "mod_par_cor_1e5.Rds"))


par_TT <- params$TT
par_PA <- params$PA
par_SQ <- params$SQ
par_AT <- params$AT
par_UN <- params$UN

DOYs <- hb_hf_li$doy

sites <- hb_hf_li$site %>%
    unique()



# ~ Calculate predtion data tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CalPredDt <- function(fun, par_dt) {
    pred_dt <- lapply(1:nrow(par_dt), function(i) {
            pred <- do.call(
                fun,
                list(par = unlist(par_dt[i, ]), data = hb_hf_li)
            )
            pred[pred == 9999] <- NA
            return(pred)
        }) %>%
            do.call(cbind, .) %>%
            # Calculate median and CI
            apply(1, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
            t() %>%
            data.table() %>%
            .[, PhenoYear := hb_hf_li$year] %>%
            .[, siteID := hb_hf_li$site] %>%
            setcolorder(c("siteID", "PhenoYear")) %>%
            set_colnames(c("siteID", "PhenoYear", "lwr", "med", "upr")) %>%
            # Calculate uncertainty range
            .[, unc_range := upr - lwr]
}

TT_pred_dt <- CalPredDt(ThermalTimeModel, par_TT)
PA_pred_dt <- CalPredDt(ParallelModel, par_PA)
SQ_pred_dt <- CalPredDt(SequentialModel, par_SQ)
AT_pred_dt <- CalPredDt(AlternatingModel, par_AT)
UN_pred_dt <- CalPredDt(UnifiedModel, par_UN)


# ~ Make the plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PlotPredCI <- function(site_name) {

    PlotModel <- function(pred_dt, color) {
        site_ha <- pred_dt[siteID == site_name]
        lines(site_ha[, .(PhenoYear, med)],
            type = "l",
            col = color[1],
            lwd = 2
        )
        polygon(c(site_ha$PhenoYear, rev(site_ha$PhenoYear)),
            c(site_ha$lwr, rev(site_ha$upr)),
            border = NA,
            col = adjustcolor(color[1], alpha.f = 0.1)
        )
    }

    PlotModel(TT_pred_dt, jenna_pal[1])
    PlotModel(PA_pred_dt, jenna_pal[2])
    PlotModel(SQ_pred_dt, jenna_pal[3])
    PlotModel(AT_pred_dt, jenna_pal[4])
    PlotModel(UN_pred_dt, jenna_pal[5])

    legend("bottomleft", 
        bty = "n", 
        lwd = 2,
        lty = rep(1, 5),
        col = jenna_pal,
        legend = c("TT", "PA", "SQ", "AT", "UN")
    )
}




# fig: model_pred_uncertainty.png
png(file.path("Output/Assets", "model_pred_uncertainty.png"), 
    width = 1500, height = 1500, res = 150
)

par(mfrow = c(3, 2), mar = c(4, 4, 3, 2), oma = c(1, 1, 0, 0), 
    xpd = NA,
    bty = "L"
)
for (i in seq_along(sites)) {
    plot(NA,
        xlim = range(hb_hf_li$year), ylim = c(120, 170),
        xlab = "Year", ylab = "Predicted SOS (DOY)",
        las = 1,
        cex.lab = 1.5,
        cex.axis = 1.5
    )
    PlotPredCI(sites[i])
    title(paste0("Site:", sites[i]), line = -1, adj = 0.9)
}

boxplot(NA, xlim = c(0.5, 5.5), ylim = c(0, 20), 
    xlab = "Site", 
    ylab = "Uncertainty range (days)",
    cex.lab = 1.5,
    cex.axis = 1.5
)
axis(side = 1, at = 1:5, labels = sites, cex.axis = 1.5)

for (i in seq_along(sites)) {
    PlotModelUncBox <- function(st, pred_dt, color, at) {
        boxplot(pred_dt[siteID == st, unc_range],
            add = TRUE,
            yaxt = "n",
            col = color,
            boxwex = 0.15,
            at = at
        )
    }
    
    at_idx <- seq(-0.25, 0.25, len = 5) + i

    PlotModelUncBox(sites[i], TT_pred_dt, jenna_pal[1], at_idx[1])
    PlotModelUncBox(sites[i], PA_pred_dt, jenna_pal[2], at_idx[2])
    PlotModelUncBox(sites[i], SQ_pred_dt, jenna_pal[3], at_idx[3])
    PlotModelUncBox(sites[i], AT_pred_dt, jenna_pal[4], at_idx[4])
    PlotModelUncBox(sites[i], UN_pred_dt, jenna_pal[5], at_idx[5])
}

dev.off()
