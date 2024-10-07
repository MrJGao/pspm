#'******************************************************************************
#' Description: Use survial models to test chilling in the simulated data.
#'******************************************************************************

rm(list=ls())

source("src/base.R")
source("src/hlp/hlp_chilling_calc.R")
source("src/mod_hpc_pspm_spec_sim_base.R")

library(survival)


# Calculate chilling and forcing from a person-period dataset
CalChillForc <- function(pp_dt) {
    # Calculate chilling
    Tc <- 0
    tc <- -20

    # Let's use chill days for now
    pp_dt[doy > tc,
        chill := .SD[, ChillFix(Tmean, "h")],
        by = .(id)
    ]
    # Accumulated chilling
    pp_dt[doy > tc,
        acc_chill := .SD[, cumsum(chill)],
        by = .(id)
    ]


    # Calculate forcing
    Tb <- 5
    t0 <- 1
    pp_dt[doy > t0,
        forc := .SD[
            ,
            ifelse(Tmean - Tb > 0, Tmean - Tb, 0)
        ],
        by = .(id)
    ]
    # Accumulated forcing
    pp_dt[doy > t0,
        acc_forc := .SD[, cumsum(forc)],
        by = .(id)
    ]

    return(pp_dt)
}



# Deviance difference for each case, run 100 iterations
dev_diff <- data.table(iter = numeric(), case = numeric(), dif = numeric())
for (j in 1:100) {
    for (i in 1:14) {
        data_li <- mix_li[[i]]
        # Add random noise
        data_li$transition_dates <- data_li$transition_dates +
            rnorm(length(data_li$transition_dates), 0, 1)
        data_li$transition_dates <- round(data_li$transition_dates)

        # Prepare data for survival analysis
        sos_dt <- PhenoModelDataToDT(data_li)
        # Person-period full
        pp_dt_full <- ToPersonPeriod(sos_dt)
        # Person-level table format
        pl_dt <- ToPersonLevel(pp_dt_full)
        # Calculate chilling and forcing
        pp_dt_full <- CalChillForc(pp_dt_full)
        # Person-period format for survival analysis
        pp_dt <- ToPersonPeriodForSurv(pp_dt_full)

        # Caculate survival object
        ts <- survfit(Surv(doy, event) ~ 1, data = pl_dt)

        # Survival model fit
        afit1 <- glm(event ~ acc_forc * doy,
            family = "binomial", data = pp_dt[doy >= ts$time[1], ]
        )
        afit2 <- glm(event ~ acc_forc * doy + acc_chill,
            family = "binomial", data = pp_dt[doy >= ts$time[1], ]
        )

        ano <- anova(afit1, afit2)
        dif <- ano$Deviance[2]

        dev_diff <- rbind(dev_diff, data.table(
            iter = j,
            case = i,
            dif = dif
        ))
    }
}

# out: Export the survival simulation result
saveRDS(dev_diff, 
    file.path(gdir, "Pipeline", "pspm_spec_sim_survival.Rds")
)
