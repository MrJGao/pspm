#'******************************************************************************
#' Description: Chilling calculating methods
#'******************************************************************************



#' Chill days. Number of days temperature is below a threshold.
#' 
#' @param temp_vec Temperature vector.
#' @param Tb Base temperature threshold.
#' @param above0 Logical, whether only consider temperatures above 0. 
#' @return Numeric vector. Chill days
ChillDays <- function(temp_vec, Tb, above0 = FALSE) {
    if (above0 == FALSE) {
        return(as.numeric(temp_vec < Tb))
    } else {
       return(as.numeric(temp_vec < Tb & temp_vec > 0))
    }
}


#' Chilling units, the units below a temperature threshold.
#' 
#' @param temp_vec Temperature vector.
#' @param Tb Base temperature.
#' @return Numeric vector. Chilling units.
ChillUnit <- function(temp_vec, Tb) {
    cu <- temp_vec
    cu <- ifelse(cu < Tb & cu > 0, Tb - cu, 0)
    
    return(cu)
}


#' Triangular chilling response.
#' 
#' @param temp_vec Temperature vector.
#' @param Topt,Tmin,Tmax Temperature thresholds in the triangular function.
#' @return Numeric vector. The chilling response unit.
TriChill <- function(temp_vec, Topt, Tmin, Tmax) {
    # sanity checks
    if (Topt >= Tmax || Topt <= Tmin || Tmax <= Tmin) {
    return(NA)
    }

    # find locations of rising and falling
    # part of the triangular function
    loc_rising <- which(temp_vec < Topt & temp_vec >= Tmin)
    loc_falling <- which(temp_vec <= Tmax & temp_vec >= Topt)

    # fill this vector according to a triangular
    # ruleset function

    # set out of range values to 0
    temp_vec[temp_vec < Tmin | temp_vec > Tmax] <- 0

    # convert temperatures
    temp_vec[loc_rising] <- (temp_vec[loc_rising] - Tmin) / (Topt - Tmin)
    temp_vec[loc_falling] <- 1 - (temp_vec[loc_falling] - Topt) / 
        (Tmax - Topt)

    # return the coverted temperature vector
    return(temp_vec)
}


#' Chilling accumulation based on fixed temperature thresholds.
#' 
#' @param temp_vec Temperature vector.
#' @param method The chosen method to use. Possible values are "a", "b", "c", ...
#' @return Numeric vector. The chilling response unit.
ChillFix <- function(temp_vec, method) {
    if (method == "a") {
        return(as.numeric(temp_vec < 0))
    } else if (method == "b") {
        return(as.numeric(temp_vec < 5))
    } else if (method == "c") {
        return(as.numeric(temp_vec > -10 & temp_vec < 5))
    } else if (method == "d") {
        return(as.numeric(temp_vec > 0 & temp_vec < 5))
    } else if (method == "e") {
        return(as.numeric(temp_vec < 7))
    } else if (method == "f") {
        return(as.numeric(temp_vec > -10 & temp_vec < 7))
    } else if (method == "g") {
        return(as.numeric(temp_vec > 0 & temp_vec < 7))
    } else if (method == "h") { # The Utah model
        chill_units <- sapply(temp_vec, function(i) {
            if (i <= 1.4) {
                return(0)
            }
            if (i > 1.4 & i <= 2.4) {
                return(0.5)
            }
            if (i > 2.4 & i <= 9.1) {
                return(1)
            }
            if (i > 9.1 & i <= 12.4) {
                return(0.5)
            }
            if (i > 12.4 & i <= 15.9) {
                return(0)
            }
            if (i > 15.9 & i <= 18) {
                return(-0.5)
            }
            if (i > 18) {
                return(0)
            }
        })
        return(chill_units)
    }
}



