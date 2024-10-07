# ******************************************************************************
# Fake chilling-forcing relationship simulation by using PEP725+E-OBS data
# ******************************************************************************
source("src/base.R")
LoadHelpers()
library(data.table)
library(magrittr)
library(lubridate)
library(parallel)


# ~ Simulate Fcrit ####
# ~ ----------------------------------------------------------------------------

dt_files <- list.files(file.path("data/raw", "pep725_pheno_temperature"), 
    pattern = ".csv$", full.names = TRUE
)

# Use one species as an example for now
spec_dt <- fread(dt_files[1])
# Number of site-years
spec_dt[, .(PEP_ID)] %>% unique() %>% nrow()

uniqueN(spec_dt$PEP_ID)

# Randomly sample 100 sites out
# sample_id <- sample(spec_dt$PEP_ID, 100)
# spec_dt <- spec_dt[PEP_ID %in% sample_id, ]

spec_li <- FormatDataForPhenoModel(
    spec_dt[, .(siteID = PEP_ID, Date = as_date(Date), Tmean, SOS, PhenoYear)]
)


# Prescribed parameters
t0 <- 30 + 110
T_base <- 5
F_crit <- 200

sim_sos <- do.call(ThermalTimeModel, 
    list(par = c(t0, T_base, F_crit), data = spec_li)
)

sim_spec_li <- copy(spec_li)
sim_spec_li$T_min <- sim_spec_li$Ti
sim_spec_li$T_max <- sim_spec_li$Ti
sim_spec_li$transition_dates <- sim_sos
# Remove invalidate site-years
sim_spec_li <- GetDataLiIndex(sim_spec_li, which(sim_sos != 9999))


# Conver back to data.table
sim_spec_dt <- Conv2Dt(sim_spec_li)

saveRDS(sim_spec_dt, "pipe/sim_spec_aesculus_dt.Rds")



# ~ Use TT model to fit PEP725 to get Fcrit ####
# ~ ----------------------------------------------------------------------------
# Now, fix t0 and Tbase
TTfix <- function(par, data, diagnose = FALSE) {
    # extract the parameter values from the
    # par argument for readability
    t0 <- 110 + 30
    T_base <- 5
    F_crit <- par

    # simple degree day sum setup
    Rf <- data$Ti - T_base
    Rf[Rf < 0] <- 0
    Rf[1:t0, ] <- 0

    # DOY of budburst criterium
    doy <- apply(Rf, 2, function(xt) {
        data$doy[which(cumsum(xt) >= F_crit)[1]]
    })

    doy <- StandardizeOutput(doy)

    if (diagnose == TRUE) {
        return(list(doy = doy, Rf = Rf))
    } else {
        return(doy)
    }
}



dt_files <- list.files(file.path("data/raw", "pep725_pheno_temperature"), 
    pattern = ".csv$", full.names = TRUE
)

# Use one species as an example for now
spec_dt <- fread(dt_files[1])
# Number of site-years
spec_dt[, .(PEP_ID)] %>% unique() %>% nrow()

sites <- unique(spec_dt$PEP_ID)

cl <- makeCluster(detectCores() - 1)
calls <- clusterCall(cl, function() {
    suppressWarnings({
        source("src/base.R")
        LoadHelpers()
        library(data.table)
        library(magrittr)
        library(lubridate)
        library(parallel)
    })
})
clusterExport(cl, c("spec_dt", "TTfix"))

sim_spec_dt2 <- clusterApplyLB(cl, sites, function(st) {
    # For one site
    site_dt <- spec_dt[PEP_ID == st,]

    site_li <- FormatDataForPhenoModel(
        site_dt[, .(siteID = PEP_ID, Date = as_date(Date), Tmean, SOS, PhenoYear)]
    )

    # Fit the model
    est_par <- GenSA::GenSA(
        par = NULL,
        model = TTfix,
        fn = CostRMSE,
        data = site_li,
        lower = c(0),
        upper = c(20000),
        control = NULL
    )
    est <- do.call(TTfix, list(par = est_par$par, data = site_li))

    site_li$transition_dates <- est

    site_dt <- Conv2Dt(site_li)

    return(site_dt)
})
sim_spec_dt2 <- do.call(rbind, sim_spec_dt2)


saveRDS(sim_spec_dt2, "pipe/sim_spec_aesculus_dt_est_fcrit.Rds")


