#'******************************************************************************
#' Description: Spring phenology model functions. THIS IS NOT MY ORIGINAL!
#' I referred the `phenor` package when writing these functions. 
#' Actually, many functions were directly copied from the `phenor` package.
#'******************************************************************************

# Some common parameter values:
# `t0`: starting day of the heat sum calculation.
# `T_base`: base temperature.
# `F_crit`: threshold of state of forcing to budburst/leaf-out/flowering.
# `C_crit`: threshold of state of chilling to transit from rest to quiescence.


#' Standardize output day-of-year values.
#' @param doy DOY vector.
#' @return Standardized DOY vector.
#' @export
StandardizeOutput <- function(doy) {
    doy[is.na(doy)] <- 9999
    doy[is.infinite(doy)] <- 9999
    return(doy)
}



# ~ Temperature responses ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Growing degree day temperature response.
#' As defined in Basler et al 2016.
#' @param T a vector or matrix of temperatures
#' @param T_base Base temperature threshold
#' @return Temperature responses
#' @export
#' @example
#' /dontrun {
#' T <- -30:30
#' T_response <- GddTempRsp(T = T, T_base = 10)
#' plot(T, T_response)
#' }
GddTempRsp <- function(T, T_base) {
    T[T <= T_base] <- 0
    T[T > T_base] <- (T - T_base)[T > T_base]
    return(T)
}



#' Triangular temperature response function. (Copied from `phenor`)
#' As defined in: Basler et al. 2016 (AgForMet)
#'
#' @param T a vector or matrix of temperatures
#' @param T_opt optimal temperature
#' @param T_min minimum viable temperature
#' @param T_max maximum viable temperature
#' @keywords phenology, model, temperature response
#' @export
#' @examples
#'
#' T_response <- triangular_temperature_response(T = 0:45)
#' \dontrun{
#' plot(0:45, T_response, type = "l")
#' }
triangular_temperature_response <- function(T = -10:45, 
    T_opt = 10, T_min = 1, T_max = 15
) {
    # sanity checks
    if (T_opt >= T_max || T_opt <= T_min || T_max <= T_min) {
        T[] <- NA
        return(T)
    }

    # find locations of rising and falling
    # part of the triangular function
    loc_rising <- which(T < T_opt & T >= T_min)
    loc_falling <- which(T <= T_max & T >= T_opt)

    # fill this vector according to a triangular
    # ruleset function

    # set out of range values to 0
    T[T < T_min | T > T_max] <- 0

    # convert temperatures
    T[loc_rising] <- (T[loc_rising] - T_min) / (T_opt - T_min)
    T[loc_falling] <- 1 - (T[loc_falling] - T_opt) / (T_max - T_opt)

    # returns the converted temperature data
    return(T)
}


#' Bell shaped temperature response function.
#' As defined in: Basler et al 2016.
#' @param T A vector or matrix of temperatures.
#' @param a Bell shape parameter.
#' @param b Bell shape parameter.
#' @param c Bell shape parameter.
#' @return Temperature responses.
#' @export
#' @example
#' /dontrun {
#' T_response <- BellshpTempRsp(T = c(-20:30), 1, 1, 1)
#' plot(T_response)
#' }
BellShpTempRsp <- function(T, a, b, c) {
    T <- 1 / (1 + exp(a * (T - c)^2) + b * (T - c))
    return(T)
}


#' Sigmoid temperature response function.
#' As defined in: Basler et al 2016.
#' @param T A vector or matrix of temperatures.
#' @param b Sigmoid function parameter.
#' @param c Sigmoid function parameter.
#' @return Temperature responses.
#' @export
#' @example
#' /dontrun {
#' T_response <- SigmoidTempRsp(T = c(-20:30), 1, 1)
#' plot(T_response)
#' }
SigmoidTempRsp <- function(T, b, c) {
    T <- 1 / (1 + exp(-b * (T - c)))
    return(T)
}



# ~ Models ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Linear regression model with preseason mean temperature. The period of 
#' preseason is estimated in the model.
#' 
#' NOTE: Cannot use `do.call` to predict from this model as it depends on 
#' training dataset to get the slope and intercept of the linear regression. 
#' 
#' @param par Ordered parameter vector including `t1` and `t2`, representing 
#' start and end of the preseason period used to compute mean temperature.
#' @param data Data containing spring phenophase transition dates and the.
#' @return Fitted phenophase transition dates, in day of year (DOY) format.
#' @export
LIN <- function(par, data) {
    t1 <- round(par[1])
    t2 <- round(par[2])

    if (t1 >= t2) {
        return(9999)
    }

    T_mean <- apply(data$Ti, 2, function(Ti) {
        mean(Ti[t1:t2])
    })
    mod_lm <- lm(data$transition_dates ~ T_mean)
    doy <- fitted(mod_lm)

    doy <- StandardizeOutput(doy)

    return(doy)
}

#' Thermal time or spring warming model
#' @param par Ordered parameter vector including `t0`, `T_base`, `F_crit`. 
#' `t0` is the first day to count; `T_base` is the base temperature threshold; 
#' `F_crit` is the critical threshold reaching which a spring phenophase 
#' transition would occur.
#' @param data Data containing spring phenophase transition dates and the 
#' corresponding temperature records
#' @param diagnose (logical) `TRUE` will return diagnostic information. 
#' In this case, calculated forcing degree days will be returned.
#' @return Fitted phenophase transition dates, in day of year (DOY) format.
#' @export
ThermalTimeModel <- function(par, data, diagnose = FALSE) {
    # extract the parameter values from the
    # par argument for readability
    t0 <- round(par[1])
    T_base <- par[2]
    F_crit <- par[3]

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


#' Thermal Time model with fixed t0 and Tb
#' 
#' @param par Ordered parameter vector but only `F_crit` should be included.
#' @param data Data containing spring phenophase transition dates and the
#'   corresponding temperature records
#' @return itted phenophase transition dates, in day of year (DOY) format.
 #' @export
TTfix <- function(par, data, t0, T_base) {
    if (length(par) != 1) {
        stop("model parameter(s) out of range (too many, too few)")
    }

    F_crit <- par[1]

    # simple degree day sum setup
    Rf <- data$Ti - T_base
    Rf[Rf < 0] <- 0
    Rf[1:t0, ] <- 0

    # DOY of budburst criterium
    doy <- apply(Rf, 2, function(xt) {
        data$doy[which(cumsum(xt) >= F_crit)[1]]
    })

    doy <- StandardizeOutput(doy)
    
    return(doy)
}


#' Parallel model, which assumes forcing and chilling happen concurrently. 
#' @param par Ordered parameter vector including `t0`, `t0_chill`, `T_base`, 
#' `T_opt`, `T_min`, `T_max`, `C_min`, `F_crit`, `C_req`
#' @param data Data containing spring phenophase transition dates and the 
#' corresponding temperature records
#' @param diagnose (logical) `TRUE` will return diagnostic information. In this 
#' case, calculated forcing degree days, chilling units and so on will be 
#' returned.
#' @return Fitted phenophase transition dates, in day of year (DOY) format.
#' @export
ParallelModel <- function(par, data, diagnose = FALSE) {
    # exit the routine as some parameters are missing
    if (length(par) != 9) {
        stop("model parameter(s) out of range (too many, too few)")
    }

    # get the parameter values
    t0 <- round(par[1]) # first day to count forcing
    t0_chill <- round(par[2]) # first day to count chilling
    T_base <- par[3] # base temperature
    T_opt <- par[4] # optimal temperature
    T_min <- par[5] # daily minimum temperature
    T_max <- par[6] # daily maximum temperature
    C_min <- par[7] # minimum growth competence that unchilled buds can have
    F_crit <- par[8] # forcing criteria
    C_req <- par[9] # chilling requirement

    if (t0_chill > t0) {
        return(rep(9999, ncol(data$Ti)))
    }

    # chilling
    Rc <- triangular_temperature_response(data$Ti,
        T_opt = T_opt,
        T_min = T_min,
        T_max = T_max
    )
    Rc[1:t0_chill, ] <- 0
    Sc <- apply(Rc, 2, cumsum)

    # growth competence function
    k <- C_min + Sc * (1 - C_min) / C_req
    k[Sc >= C_req] <- 1

    # forcing
    Rf <- data$Ti - T_base
    Rf[Rf < 0] <- 0
    Rf <- Rf * k
    Rf[1:t0, ] <- 0

    # DOY of budburst criterium
    doy <- apply(Rf, 2, function(xt) {
        data$doy[which(cumsum(xt) >= F_crit)[1]]
    })

    doy <- StandardizeOutput(doy)

    if (diagnose == TRUE) {
        return(list(doy = doy, k = k, Rf = Rf, Sc = Sc))
    } else {
       return(doy)
    }
}

ParallelModelFixNov1 <- function(par, data, diagnose = FALSE) {
    # exit the routine as some parameters are missing
    if (length(par) != 8) {
        stop("model parameter(s) out of range (too many, too few)")
    }

    # get the parameter values
    # t0 <- round(par[1]) # first day to count forcing
    t0 <- 305 - 110 # Fix to Nov 1st.
    t0_chill <- round(par[1]) # first day to count chilling
    T_base <- par[2] # base temperature
    T_opt <- par[3] # optimal temperature
    T_min <- par[4] # daily minimum temperature
    T_max <- par[5] # daily maximum temperature
    C_min <- par[6] # minimum growth competence that unchilled buds can have
    F_crit <- par[7] # forcing criteria
    C_req <- par[8] # chilling requirement

    if (t0_chill > t0) {
        return(rep(9999, ncol(data$Ti)))
    }

    # chilling
    Rc <- triangular_temperature_response(data$Ti,
        T_opt = T_opt,
        T_min = T_min,
        T_max = T_max
    )
    Rc[1:t0_chill, ] <- 0
    Sc <- apply(Rc, 2, cumsum)

    # growth competence function
    k <- C_min + Sc * (1 - C_min) / C_req
    k[Sc >= C_req] <- 1

    # forcing
    Rf <- data$Ti - T_base
    Rf[Rf < 0] <- 0
    Rf <- Rf * k
    Rf[1:t0, ] <- 0

    # DOY of budburst criterium
    doy <- apply(Rf, 2, function(xt) {
        data$doy[which(cumsum(xt) >= F_crit)[1]]
    })

    doy <- StandardizeOutput(doy)

    if (diagnose == TRUE) {
        return(list(doy = doy, k = k, Rf = Rf, Sc = Sc))
    } else {
        return(doy)
    }
}


#' Sequential model from  `phenor`
#' @param par Ordered parameter vector including: `t0`, `t0_chill`, `T_base`, 
#' `T_opt`, `T_min`, `T_max`, `F_crit`, `C_req`.
#' @param data Data containing spring phenophase transition dates and the 
#' corresponding temperature records.
#' @param diagnose (logical) `TRUE` will return diagnostic information.
#' @return Fitted phenophase transition dates, in day of year (DOY) format.
#' @export
SequentialModel <- function(par, data, diagnose = FALSE) {
    # exit the routine as some parameters are missing
    if (length(par) != 8) {
        stop("model parameter(s) out of range (too many, too few)")
    }

    # extract the parameter values from the par argument for readability
    t0 <- round(par[1])
    t0_chill <- round(par[2])
    T_base <- par[3]
    T_opt <- par[4]
    T_min <- par[5]
    T_max <- par[6]
    F_crit <- par[7]
    C_req <- par[8]

    # sanity check t0 always comes after t0_chill
    if (t0 <= t0_chill) {
        return(rep(NA, ncol(data$Ti)))
    }

    # chilling
    Rc <- triangular_temperature_response(
        data$Ti,
        T_opt = T_opt,
        T_min = T_min,
        T_max = T_max
    )

    Rc[1:t0_chill, ] <- 0
    Rc[t0:nrow(Rc), ] <- 0
    Sc <- apply(Rc, 2, cumsum)

    # chilling requirement has to be met before
    # accumulation starts (binary choice)
    k <- as.numeric(Sc >= C_req)

    # forcing
    Rf <- data$Ti - T_base
    Rf[Rf < 0] = 0
    Rf <- Rf * k
    Rf[1:t0, ] <- 0 # CHECK THIS IN LITERATURE

    # DOY of budburst criterium
    doy <- apply(Rf, 2, function(xt) {
        data$doy[which(cumsum(xt) >= F_crit)[1]]
    })

    doy <- StandardizeOutput(doy)

    if (diagnose == TRUE) {
        return(list(doy = doy, Sc = Sc, Rf = Rf, k = k))
    }
    
    return(doy)
}

SequentialModelFixNov1 <- function(par, data, diagnose = FALSE) {
    # exit the routine as some parameters are missing
    if (length(par) != 7) {
        stop("model parameter(s) out of range (too many, too few)")
    }

    # extract the parameter values from the par argument for readability
    t0 <- 305 - 110
    t0_chill <- round(par[1])
    T_base <- par[2]
    T_opt <- par[3]
    T_min <- par[4]
    T_max <- par[5]
    F_crit <- par[6]
    C_req <- par[7]

    # sanity check t0 always comes after t0_chill
    if (t0 <= t0_chill) {
        return(rep(NA, ncol(data$Ti)))
    }

    # chilling
    Rc <- triangular_temperature_response(
        data$Ti,
        T_opt = T_opt,
        T_min = T_min,
        T_max = T_max
    )

    Rc[1:t0_chill, ] <- 0
    Rc[t0:nrow(Rc), ] <- 0
    Sc <- apply(Rc, 2, cumsum)

    # chilling requirement has to be met before
    # accumulation starts (binary choice)
    k <- as.numeric(Sc >= C_req)

    # forcing
    Rf <- data$Ti - T_base
    Rf[Rf < 0] = 0
    Rf <- Rf * k
    Rf[1:t0, ] <- 0 # CHECK THIS IN LITERATURE

    # DOY of budburst criterium
    doy <- apply(Rf, 2, function(xt) {
        data$doy[which(cumsum(xt) >= F_crit)[1]]
    })

    doy <- StandardizeOutput(doy)

    if (diagnose == TRUE) {
        return(list(doy = doy, Sc = Sc, Rf = Rf, k = k))
    }
    
    return(doy)
}

#' Sequential model defined in Hanninen et al 1990.
#' @param par Ordered parameter vector including: `t0`, `t0_chill`, `T_base`, 
#' `T_opt`, `T_min`, `T_max`, `F_crit`, `C_req`.
#' @param data Data containing spring phenophase transition dates and the 
#' corresponding temperature records.
#' @param diagnose (logical) `TRUE` will return diagnostic information.
#' @return Fitted phenophase transition dates, in day of year (DOY) format.
#' @export
SequentialModel_J <- function(par, data, diagnose = FALSE) {
    # exit the routine as some parameters are missing
    if (length(par) != 7) {
        stop("model parameter(s) out of range (too many, too few)")
    }

    # extract the parameter values from the par argument for readability
    t0 <- round(par[1])
    C_req <- par[2]
    T_base <- par[3]
    T_opt <- par[4]
    T_min <- par[5]
    T_max <- par[6]
    F_crit <- par[7]

    # chilling
    Rc <- triangular_temperature_response(
        data$Ti,
        T_opt = T_opt,
        T_min = T_min,
        T_max = T_max
    )

    Rc[1:t0, ] <- 0
    Sc <- apply(Rc, 2, cumsum)
    t0_frc <- apply(Sc, 2, function(Sci) {
        # When chilling is met, start to accumulate forcing
        which(Sci >= C_req)[1] 
    })
    
    # forcing
    Rf <- data$Ti - T_base
    k <- matrix(NA, nrow = 365, ncol = ncol(Sc))
    for (i in 1:length(t0_frc)) {
        if (is.na(t0_frc[i])) { # chilling hasn't met
            Rf[, i] <- 0
        } else { # chilling has met
            # chilling requirement has to be met before
            # accumulation starts (binary choice)
            k[, i] <- ifelse(Sc[, i] >= C_req, 1, 0)

            Rf[Rf < 0] <- 0
            Rf[1:t0_frc[i], i] <- 0 # CHECK THIS IN LITERATURE
            Rf[, i] <- Rf[, i] * k[, i]
        }
    }
    
    # DOY of budburst criterium
    doy <- apply(Rf, 2, function(xt) {
        data$doy[which(cumsum(xt) >= F_crit)[1]]
    })

    doy <- StandardizeOutput(doy)

    if (diagnose == TRUE) {
        return(list(doy = doy, Sc = Sc, Rf = Rf, k = k, t0_frc = t0_frc))
    }
    
    return(doy)
}


#' Alternating model.
#' 
#' @param par Ordered parameter vector including: `t0`, `T_base`, `a`, `b`, `c`.
#' @param data Data containing spring phenophase transition dates and the
#' corresponding temperature records.
#' @param diagnose (logical) `TRUE` will return diagnostic information.
#' @return Fitted phenophase transition dates, in day of year (DOY) format.
#' @export
AlternatingModel <- function(par, data, diagnose = FALSE) {
    # exit the routine as some parameters are missing
    if (length(par) != 5) {
        stop("model parameter(s) out of range (too many, too few)")
    }

    # extract the parameter values
    t0 <- round(par[1])
    T_base <- par[2]
    a <- par[3]
    b <- par[4]
    c <- par[5]

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


#' The Unified model by chuine 2000.
#' 
#' @param par Ordered parameter vector including: 
#' @param data Data containing spring phenophase transition dates and the 
#' corresponding temperature records.
#' @param diagnose (logical) `TRUE` will return diagnostic information.
#' @return Fitted phenophase transition dates, in day of year (DOY) format.
#' @export
UnifiedModel <- function(par, data, diagnose = FALSE) {
    # This is an effort to reproduce the Unified Model
    # of Chuine 2000 in full form (not simplified)

    # exit the routine if parameters are missing
    if (length(par) != 9) {
        stop("model parameter(s) out of range (too many, too few)")
    }

    CF <- function(x, a_c, b_c, c_c) {
        1 / (1 + exp(a_c * (x - c_c)^2 + b_c * (x - c_c)))
    }

    # extract the parameter values from the
    # par argument in a more human readable form
    tc <- round(par[1]) # doy until when to accumulate chilling days
    a_c <- par[2] # sigmoid function chilling parameter a
    b_c <- par[3] # sigmoid function chilling parameter b
    c_c <- par[4] # sigmoid function chilling parameter c
    b_f <- par[5] # sigmoid function forcing parameter b
    c_f <- par[6] # sigmoid function forcing parameter c
    w <- par[7] # F* parameter w
    k <- par[8] # F* parameter k
    C_req <- par[9] # Chilling degree threshold requirement
    # i.e. C* whatever it is called

    # chilling accumulation using the
    # triangular shaped temperature response
    # basically convert normal temperatures to
    # what is called "chilling units" in the paper
    Rc <- CF(x = data$Ti, a_c, b_c, c_c)

    # Set values  < 0 to 0 (shouldn't count)
    # set NA values to 0 (when the output of the
    # triangular response is "empty" i.e. NA)
    Rc[is.na(Rc)] <- 0
    Rc[Rc < 0] <- 0

    # accumulate the chilling units, this is a
    # cummulative sum so you have all values along
    # the time axis (this saves time / iterations)
    Sc <- apply(Rc, 2, cumsum)

    # chilling requirement has to be met before
    # accumulation starts (binary choice) basically
    # binary mask to be applied to the Forcing temperature
    # data (sets anything before C_req to 0)
    m <- apply(Sc >= C_req, 2, as.numeric)

    # calculates when (row number) C_req is met
    row_loc <- apply(m, 2, function(x) which(x == 1)[1])

    # if any row_loc is NA (C_req not met) or all of m == 0 
    # (same deal / fallback to be sure) skip the rest as you won't be able to 
    # set C_tot which by default should be > C_req.
    if (any(is.na(row_loc)) | all(m == 0)) {
        return(StandardizeOutput(doy = rep(9999, ncol(Sc))))
    }

    # if all columns have valid values and the associate time 
    # location (row value), check if the tc value which determines the total 
    # chilling degree day accumulation exceeds the maximum value, if not 
    # C_tot < C_req which is not allowed the total is always equal to or greater
    # than C_req. Skip if the condition is not met
    if (tc < max(row_loc)) {
        return(StandardizeOutput(doy = rep(9999, ncol(Sc))))
    }

    # if above conditions are met
    # select the row which defines C_tot
    # this is non-dynamic across all sites / years
    # (as far as I can deduce from the Chuine paper)
    C_tot <- Sc[tc, ]

    # apply the unified CF function
    # with a parameter set to 0
    Rf <- CF(x = data$Ti, 0, b_f, c_f)

    # Apply the chilling mask to forcing
    # temperature values
    Rf <- Rf * m

    # cummulate the forcing values
    Sfc <- apply(Rf, 2, cumsum)

    # calculate the Forcing requirement
    # based upon the C_tot value
    F_req <- w * exp(k * C_tot)

    # Trap invalid F* values, no need to waste additional cycles
    if (any(is.na(F_req)) | any(is.infinite(F_req))) {
        return(StandardizeOutput(doy = rep(9999, ncol(Sc))))
    }

    # take the difference between the
    # forcing matrix and one filled with
    # the required F* values, where it
    # exceeds 0 first is the day of
    # leaf development
    Sfc <- sweep(Sfc, 2, F_req, FUN = "-")
    doy <- apply(Sfc, 2, function(x) {
        data$doy[which(x > 0)[1]]
    })

    # set export format, either a rasterLayer or a vector
    return(StandardizeOutput(doy))
}

UnifiedModel2 <- function(par, data, diagnose = FALSE) {
    # This is an effort to reproduce the Unified Model
    # of Chuine 2000 in full form (not simplified)

    # exit the routine if parameters are missing
    if (length(par) != 8) {
        stop("model parameter(s) out of range (too many, too few)")
    }

    CF <- function(x, a_c, b_c, c_c) {
        1 / (1 + exp(a_c * (x - c_c)^2 + b_c * (x - c_c)))
    }

    # extract the parameter values from the
    # par argument in a more human readable form
    tc <- round(par[1]) # doy until when to accumulate chilling days
    a_c <- par[2] # sigmoid function chilling parameter a
    b_c <- 0 # sigmoid function chilling parameter b
    c_c <- par[3] # sigmoid function chilling parameter c
    b_f <- par[4] # sigmoid function forcing parameter b
    c_f <- par[5] # sigmoid function forcing parameter c
    w <- par[6] # F* parameter w
    k <- par[7] # F* parameter k
    C_req <- par[8] # Chilling degree threshold requirement
    # i.e. C* whatever it is called

    # chilling accumulation using the
    # triangular shaped temperature response
    # basically convert normal temperatures to
    # what is called "chilling units" in the paper
    Rc <- CF(x = data$Ti, a_c, b_c, c_c) * 2

    # Set values  < 0 to 0 (shouldn't count)
    # set NA values to 0 (when the output of the
    # triangular response is "empty" i.e. NA)
    Rc[is.na(Rc)] <- 0
    Rc[Rc < 0] <- 0

    # accumulate the chilling units, this is a
    # cummulative sum so you have all values along
    # the time axis (this saves time / iterations)
    Sc <- apply(Rc, 2, cumsum)

    # chilling requirement has to be met before
    # accumulation starts (binary choice) basically
    # binary mask to be applied to the Forcing temperature
    # data (sets anything before C_req to 0)
    m <- apply(Sc >= C_req, 2, as.numeric)

    # calculates when (row number) C_req is met
    row_loc <- apply(m, 2, function(x) which(x == 1)[1])

    # if any row_loc is NA (C_req not met) or all of m == 0 
    # (same deal / fallback to be sure) skip the rest as you won't be able to 
    # set C_tot which by default should be > C_req.
    if (any(is.na(row_loc)) | all(m == 0)) {
        return(StandardizeOutput(doy = rep(9999, ncol(Sc))))
    }

    # if all columns have valid values and the associate time 
    # location (row value), check if the tc value which determines the total 
    # chilling degree day accumulation exceeds the maximum value, if not 
    # C_tot < C_req which is not allowed the total is always equal to or greater
    # than C_req. Skip if the condition is not met
    if (tc < max(row_loc)) {
        return(StandardizeOutput(doy = rep(9999, ncol(Sc))))
    }

    # if above conditions are met
    # select the row which defines C_tot
    # this is non-dynamic across all sites / years
    # (as far as I can deduce from the Chuine paper)
    C_tot <- Sc[tc, ]

    # apply the unified CF function
    # with a parameter set to 0
    Rf <- CF(x = data$Ti, 0, b_f, c_f)

    # Apply the chilling mask to forcing
    # temperature values
    Rf <- Rf * m

    # cummulate the forcing values
    Sfc <- apply(Rf, 2, cumsum)

    # calculate the Forcing requirement
    # based upon the C_tot value
    F_req <- w * exp(k * C_tot)

    # Trap invalid F* values, no need to waste additional cycles
    if (any(is.na(F_req)) | any(is.infinite(F_req))) {
        return(StandardizeOutput(doy = rep(9999, ncol(Sc))))
    }

    # take the difference between the
    # forcing matrix and one filled with
    # the required F* values, where it
    # exceeds 0 first is the day of
    # leaf development
    Sfc <- sweep(Sfc, 2, F_req, FUN = "-")
    doy <- apply(Sfc, 2, function(x) {
        data$doy[which(x > 0)[1]]
    })

    # set export format, either a rasterLayer or a vector
    return(StandardizeOutput(doy))
}


#' The modified unified model as Zhang et al 2022.
#' The only difference of this version is that it models the start date of 
#' chilling cumulation, while the original version models total chilling 
#' accumulation days (see `tc` vs `t0`).
#' 
#' @param par Ordered parameter vector including:
#' @param data Data containing spring phenophase transition dates and the
#' corresponding temperature records.
#' @return Fitted phenophase transition dates, in day of year (DOY) format.
#' @export
UnifiedModel2 <- function(par, data) {
    # exit the routine as some parameters are missing
    if (length(par) != 9) {
        stop("model parameter(s) out of range (too many, too few)")
    }

    # get the parameter values
    t0 <- round(par[1]) # chilling start cumulating date
    a_c <- par[2] # sigmoid function chilling parameter a
    b_c <- par[3] # sigmoid function chilling parameter b
    c_c <- par[4] # sigmoid function chilling parameter c
    b_f <- par[5] # sigmoid function forcing parameter b
    c_f <- par[6] # sigmoid function forcing parameter c
    w <- par[7] # F* parameter w
    k <- par[8] # F* parameter k
    C_req <- par[9] # Chilling degree threshold requirement, i.e. C*

    CF <- function(x, a, b, c) {
        1 / (1 + exp(a * (x - c)^2 + b * (x - c)))
    }

    # chilling accumulation using the triangular shaped temperature response.
    # basically convert normal temperatures to what is called "chilling units"
    # in the paper.
    Rc <- CF(x = data$Ti, a_c, b_c, c_c)

    # set values < 0 to 0 (shouldn't count)
    # set NA values to 0 (when the output of the triangular response is "empty"
    # i.e. NA)
    Rc[is.na(Rc)] <- 0
    Rc[Rc < 0] <- 0
    Rc[1:t0] <- 0

    # accumulate the chilling units
    Sc <- apply(Rc, 2, cumsum)

    # chilling requirement has to be met before accumulation starts 
    # (binary choice). basically binary mask to be applied to the Forcing 
    # temperature data (sets anything before C_req to 0)
    m <- apply(Sc >= C_req, 2, as.numeric)

    # calculates when (row number) C_req is met
    row_loc <- apply(m, 2, function(x) which(x == 1)[1])

    # if any row_loc is NA (C_req not met) or all of m == 0 
    # (same deal / fallback to be sure) skip the rest as you won't be able to 
    # set C_tot which by default should be > C_req.
    if (any(is.na(row_loc)) | all(m == 0)) {
        return(StandardizeOutput(doy = rep(9999, ncol(Sc))))
    }

    # if all columns have valid values and the associate time location 
    # (row value), check if the tc value which determines the total chilling 
    # degree day accumulation exceeds the maximum value, if not C_tot < C_req 
    # which is not allowed the total is always equal to or greater than C_req. 
    # Skip if the condition is not met
    # if (tc < max(row_loc)) {
    #     return(StandardizeOutput(doy = rep(9999, ncol(Sc))))
    # }

    # if above conditions are met
    # select the row which defines C_tot
    # this is non-dynamic across all sites / years
    # (as far as I can deduce from the Chuine paper)
    C_tot <- Sc

    # Apply the chilling mask to forcing
    # temperature values, set anything < 0 to 0 (see )
    Rf <- data$Ti * m

    # apply the unified CF function with a parameter set to 0
    Rf <- CF(x = Rf, 0, b_f, c_f)

    # cummulate the forcing values
    Sfc <- apply(Rf, 2, cumsum)

    # calculate the Forcing requirement
    # based upon the C_tot value
    F_req <- w * exp(k * C_tot)

    # Trap invalid F* values, no need to waste additional cycles
    if (any(is.na(F_req)) | any(is.infinite(F_req))) {
        return(StandardizeOutput(doy = rep(9999, ncol(Sc))))
    }

    # take the difference between the forcing matrix and one filled with the
    # required F* values, where it >= 0 first is the day of leaf development
    Sfc <- Sfc - F_req
    doy <- apply(Sfc, 2, function(x) {
        data$doy[which(x > 0)[1]]
    })

    # set export format, either a rasterLayer or a vector
    return(StandardizeOutput(doy))
}




# ~ Bayesian models ####
# ~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The model fitting procedure for Bayesian models is different from the normal
# models, so they are put in a separated section.


#' Minkyu Moon's model on Moon et al 2021, he drived the model from
#' Clark et al 2014. The funtion only returns the model object string.
#' 
#' @export
CDSOMModel <- function() {
    # The model
    cdsom_model <- "model {
        # observation level
        for(i in 1:n){
            Y[i] ~ dbern(p[i])                      # binary outcome
            yp[i] ~ dbern(p[i])                     # predictions to validate
            logit(p[i]) <- kappa + lambda * h[i]    # logit link
        }
        # development state for head nodes (i.e. start of the driving forces)
        for (i in head_nodes){
            h[i] <- 0
        }
        # development state for main nodes (i.e. all nodes except for the head nodes)
        for (i in main_nodes){
            h[i] <- h[i - 1] + max(0, X[i - 1,] %*% beta) * (1 - h[i - 1] / hmax)
        }
        kappa ~ dnorm(0, 0.001)T(, 0)

        # priors for tb and beta's
        for(i in 1:np){
            beta[i] ~ dnorm(0, 0.001) # vague prior for beta
        }
    }"

    return(cdsom_model)
}


#' The survival model devised in Elmendorf et al 2019.
#' Currently, it conly contians a logistic model without random effects.
#' 
#' @return The model object string
#' @export
ElSAModel <- function() {
    elsa_model <- "model{
        # observation level
        for(i in 1:n){
            Y[i] ~ dbern(p[i])                   # binary outcome
            yp[i] ~ dbern(p[i])                  # predictions to validate
            logit(p[i]) <- X[i, ] %*% beta
        }
        # priors
        for(i in 1:np){
            beta[i] ~ dnorm(0, 0.0001) # vague prior for beta
        }
    }"

    return(elsa_model)
}
