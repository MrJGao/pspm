
source("src/base.R")
LoadHelpers()
library(data.table)
library(magrittr)
library(lubridate)
library(parallel)
library(Rmpi)



#' Alternating model.
#' 
#' @param par Ordered parameter vector including: `t0`, `T_base`, `a`, `b`, `c`.
#' @param data Data containing spring phenophase transition dates and the
#' corresponding temperature records.
#' @param diagnose (logical) `TRUE` will return diagnostic information.
#' @return Fitted phenophase transition dates, in day of year (DOY) format.
#' @export
AlternatingModelFix <- function(par, data, diagnose = FALSE) {

    # extract the parameter values
    t0 <- 110 - 30
    T_base <- 5
    a <- 50
    b <- 5000
    c <- par

    # chilling
    Rc <- ifelse(data$Ti >= T_base, 0, 1)
    # Rc[Rc < 0] <- 1
    # Rc[Rc >= 0] <- 0
    Rc[1:t0, ] <- 0
    Sc <- apply(Rc, 2, cumsum)

    # forcing
    Rf <- data$Ti - T_base
    Rf[Rf <= 0] <- 0
    Rf[1:t0, ] <- 0
    Sf <- apply(Rf, 2, cumsum)

    Sfc <- Sf - (a + b * exp(c * Sc))

    doy <- apply(Sfc, 2, function(x) {
        data$doy[which(x > 0)[1]]
    })

    doy <- StandardizeOutput(doy)

    if (diagnose == TRUE) {
        return(list(doy = doy, Rc = Rc, Sf = Sf, Sfc = Sfc, Sc = Sc))
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
spec_dt[, .(PEP_ID)] %>%
    unique() %>%
    nrow()

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
clusterExport(cl, c("spec_dt", "AlternatingModelFix"))

sim_spec_at_dt <- clusterApplyLB(cl, sites, function(st) {
    # For one site
    site_dt <- spec_dt[PEP_ID == st, ]

    site_li <- FormatDataForPhenoModel(
        site_dt[, .(siteID = PEP_ID, Date = as_date(Date), Tmean, SOS, PhenoYear)]
    )

    # Fit the model
    est_par <- GenSA::GenSA(
        par = NULL,
        model = AlternatingModelFix,
        fn = CostRMSE,
        data = site_li,
        lower = c(-1),
        upper = c(0),
        control = NULL
    )
    est <- do.call(AlternatingModelFix, list(par = est_par$par, data = site_li))

    site_li$transition_dates <- est

    site_dt <- Conv2Dt(site_li)

    return(site_dt)
})
sim_spec_at_dt <- do.call(rbind, sim_spec_at_dt)


saveRDS(sim_spec_at_dt, "pipe/sim_spec_aesculus_dt_by_at.Rds")
