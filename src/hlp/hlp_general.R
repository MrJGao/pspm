require(R6)

General <- R6Class("General",
    private = list(
        
    ),

    public = list(
        #' Convert Kelvin to Celsius
        #' 
        #' @param k Temperature in Kelvin.
        #' @return Temperature in Celsius.
        K2C = function(k) {
            return(k - 273.15)
        },
        
        #' Two way `setdiff`
        SetDiff = function(a, b) {
            diff <- setdiff(union(a, b), intersect(a, b))
            return(diff)
        },

        # Custom AIC function which accepts loess regressions
        # Copied from `phenor` package
        AICc = function(measured, predicted, k) {
            # calculate number of observations
            n <- length(measured)

            # calculatue residual sum of squares
            RSS <- sum((measured - predicted)^2)

            # AIC
            AIC <- 2 * k + n * log(RSS / n)

            # AICc
            AICc <- AIC + (2 * k * (k + 1)) / (n - k - 1)

            # return both AIC
            return(list(
                "AIC" = AIC,
                "AICc" = AICc
            ))
        },

        #' Format data table to a list object that is more convenient for
        #' fitting sping phenology models 
        #' 
        #' @param dt A data table containing `siteID`, `PhenoYear`, `SOS`, 
        #' `Tmean`, `Tmax`, and `Tmin` 
        #' @return A list object.
        FormatDataForPhenoModel = function(dt, doy = c((-110):(-1), 1:255)) {
            # get pheno obs
            pheno_dt <- unique(dt[, .(siteID, PhenoYear, SOS)])

            # get daily mean temperature matrix
            temp_mean_mat <- apply(pheno_dt, 1, function(sy) {
                temp <- dt[siteID == sy[["siteID"]] & 
                    PhenoYear == sy[["PhenoYear"]], 
                ]
                if (nrow(temp) != 365) {
                    return(rep(NA, 365))
                }
                return(temp$Tmean)
            })

            # get daily min temperature matrix
            temp_min_mat <- apply(pheno_dt, 1, function(sy) {
                temp <- dt[siteID == sy[["siteID"]] & 
                    PhenoYear == sy[["PhenoYear"]], 
                ]
                if (nrow(temp) != 365) {
                    return(rep(NA, 365))
                }
                return(temp$Tmin)
            })

            # get daily max temperature matrix
            temp_max_mat <- apply(pheno_dt, 1, function(sy) {
                temp <- dt[siteID == sy[["siteID"]] & 
                    PhenoYear == sy[["PhenoYear"]], 
                ]
                if (nrow(temp) != 365) {
                    return(rep(NA, 365))
                }
                return(temp$Tmax)
            })

            # The result list
            li <- list()
            li$site <- pheno_dt$siteID
            li$year <- pheno_dt$PhenoYear
            li$transition_dates <- pheno_dt$SOS
            li$doy <- doy
            li$Ti <- temp_mean_mat
            li$T_min <- temp_min_mat
            li$T_max <- temp_max_mat

            # Drop site-years with NA temperature
            na_idx <- apply(li$Ti, 2, function(ti) {
                return(any(is.na(ti)))
            })

            if (any(na_idx)) {
                warning("Some site years have NA in temperature, drop them!")

                li$site <- li$site[!na_idx]
                li$year <- li$year[!na_idx]
                li$transition_dates <- li$transition_dates[!na_idx]
                li$Ti <- li$Ti[, !na_idx]
                li$T_min <- li$T_min[, !na_idx]
                li$T_max <- li$T_max[, !na_idx]

            }
            
            return(li)
        },

        #' Get the p-value of the linear regression model.
        #'
        #' @param mod_lm Linear regression model object.
        #' @return Numeric. p-value.
        GetLmPval = function(mod_lm) {
            f <- summary(mod_lm)$fstatistic
            pval <- pf(f[1], f[2], f[3], lower.tail = FALSE)

            return(pval)
        },

        #' Calculate p-value.
        #' 
        #' @param true The true values vector.
        #' @param est The estimated values vector.
        #' @return p-value.
        CalPvalue = function(true, est) {
            fit <- lm(true ~ est)
            p_val <- self$GetLMPval(fit)
            return(p_val)
        },


        #' Calculate RMSE value.
        #' 
        #' @param true The true values vector.
        #' @param est The estimated values vector.
        #' @return RMSE.
        CalRmse = function(true, est) {
            rmse <- round(sqrt(mean((true - est)^2)), 2)
            return(rmse)
        },

        
        #' Convert day of year to date.
        #'
        #' @param doy Day of year.
        #' @param year Year.
        #' @return Date format variable
        #' @export
        DoyToDate = function(doy, year) {
            date_str <- as.Date(doy, origin = paste0(year, "-01-01"))
            return(date_str)
        }
    ),

    active = list(
        
    )
)
