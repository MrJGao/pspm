# ******************************************************************************
# This figure is to show that since photoperiod does not have interannual
# variability, it may be problematic to use photoperiod to model site-specific
# phenology.
# ******************************************************************************
source("src/base.R")
source("src/vis/vis_base.R")
library(geosphere)

# Let's grab HF data as an example
hfhb_dt <- fread(file.path(hpc_dir, "Data/ARD", "hb_hf_pheno_temp.csv"))
hf_dt <- hfhb_dt[siteID == "US-Ha1"]

format(hf_dt$Date[1], "%j")
# Calculate daylength
hf_dt[, doy := format(Date, "%j")]
hf_dt[, daylen := daylength(42.5378, as.numeric(doy)), by = .(PhenoYear)]


png(file.path("Output", "photo_sos.png"), 
    width = 1200, height = 700, res = 150
)

plot(hf_dt[PhenoYear == 1992, .(Date, daylen)], 
    type = "l",
    xlab = "Date", ylab = "Daylength"
)
uni_sos <- unique(hf_dt[doy == SOS, .(Date, PhenoYear, SOS)])
uni_sos[, fake_yr_sos := as.Date(paste0("1992", "-", substr(Date, 6, 10)))]
abline(v = uni_sos$fake_yr_sos, col = "green")

dev.off()





