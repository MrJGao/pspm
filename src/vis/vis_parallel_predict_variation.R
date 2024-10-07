#'******************************************************************************
#' Description: Visualize Parallel model prediction variation result generated
#' by `mod_hpc_parallel_par_variation.R`.
#'******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")

library(data.table)
library(lubridate)
library(magrittr)


hb_dt <- fread(file.path(hpc_local_dir, 
    "Data/ARD", 
    "hb_grd_temperature_pheno.csv")
)
# unique(hb_dt[, .N / 365, by = .(siteID)])

# This site was selected b/c it has 27 site years
site_dt <- hb_dt[siteID == "1B"]

site_li <- FormatDataForPhenoModel(site_dt)


pa_var <- readRDS(file.path(
    hpc_local_dir,
    "Pipeline",
    "parallel_prediction_variation.Rds"
))


pred_dt <- lapply(pa_var, function(x) {
    x$fit
}) %>%
    do.call(cbind, .) %>%
    apply(1, quantile, c(0.025, 0.5, 0.975)) %>%
    t() %>%
    data.table() %>%
    .[, PhenoYear := site_li$year] %>%
    setcolorder("PhenoYear") %>%
    set_colnames(c("PhenoYear", "lwr", "med", "upr"))


par(bty = "L")
layout(matrix(c(1, 1, 1, 2), nrow = 1))
plot(pred_dt[, .(PhenoYear, med)],
    type = "l",
    lwd = 2,
    col = jenna_pal[2],
    xlab = "Year",
    ylab = "Predicted SOS (DOY)"
)
polygon(
    c(pred_dt$PhenoYear, rev(pred_dt$PhenoYear)),
    c(pred_dt$lwr, rev(pred_dt$upr)),
    col = adjustcolor(jenna_pal[2], alpha.f = 0.1),
    border = NA
)

boxplot(pred_dt[, upr - lwr], col = jenna_pal[2],
    ylab = "Uncertainty range (days)"
)

