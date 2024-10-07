library(R6)
library(data.table)
library(lubridate)


PhenoData <- R6Class("PhenoData", 

    public = list(
        
        # Data description
        description = "",

        initialize = function(dat = NULL, pheno_dt = NULL, temp_dt = NULL) {
            # Use data.table or phenor list to construct the object
            if (!is.null(dat)) {
                if (is.element("data.table", class(dat))) {
                    valid_cols <- c("siteID", "Date", 
                        "PhenoYear", "Tmean", "SOS"
                    )
                    find_match <- valid_cols %in% names(dat)
                    if (all(find_match == TRUE)) {
                        private$m_pheno_dt <- dat
                    } else {
                        stop("Data does not contain all needed columns!")
                    }
                } else if (is.element("list", class(dat))) {
                    tryCatch(
                        {
                            private$m_pheno_dt <- self$Conv2Dt(dat)
                        },
                        error = function(x) {
                            stop("Error initialize")
                        }
                    )
                }
            } else if (!is.null(pheno_dt) & !is.null(temp_dt)) {
                stopifnot(all(
                    c("siteID", "PhenoYear", "SOS") %in% names(pheno_dt)
                ))
                stopifnot(all(
                    c("siteID", "Date", "Tmean", "Tmax", "Tmin") %in% 
                    names(temp_dt)
                ))
                
                private$m_pheno_dt <- self$BuildTempPheno(pheno_dt, temp_dt)
            }
        },

        #' Export phenology data table to a csv file
        #' @param outfname Output filename.
        exportDT = function(outfname) {
            fwrite(private$m_pheno_dt, outfname) 
        },

        # Convert the formatted list back to data.table
        Conv2Dt = function(data_li) {
            data_dt <- data.table()
            for (i in 1:length(data_li$site)) {
                site <- data_li$site[i]
                year <- data_li$year[i]
                sos <- data_li$transition_dates[i]
                Date <- as_date(paste0(year, "-01-01")) + data_li$doy
                Tmean <- data_li$Ti[, i]
                Tmin <- data_li$T_min[, i]
                Tmax <- data_li$T_max[, i]
                data_dt <- rbind(data_dt, data.table(
                    siteID = rep(site, 365),
                    Date = Date,
                    Tmax = Tmax,
                    Tmin = Tmin,
                    Tmean = Tmean,
                    PhenoYear = year,
                    SOS = sos
                ))
            }

            return(data_dt)
        },

        #' Build phenology and temperature data for modeling
        #' 
        #' @param pheno_dt Phenology dataset.
        #' @param temp_dt Temperature dataset.
        #'
        #' @return Combined data.table
        BuildTempPheno = function(pheno_dt, temp_dt) {
            temp_dt[, Date := as_date(Date)]
            # Format TopoWx data with BLSP data
            output_dt <- lapply(1:nrow(pheno_dt), function(i) {
                cur_sy <- pheno_dt[i, ]
                # R uses a 0 based index for dates
                end_date <- as_date(
                    255 - 1,
                    origin = paste0(cur_sy[["PhenoYear"]], "-01-01")
                )
                start_date <- end_date - 364
                # Check data length
                data_len <- length(seq(start_date, end_date, by = "day"))
                stopifnot("temperature data is not 365 rows!" = data_len == 365)

                # Temperature
                cur_tp <- temp_dt[siteID == cur_sy$siteID & 
                    between(Date, start_date, end_date),
                ]

                # Phenology
                cur_tp$SOS <- pheno_dt[
                    siteID == cur_sy$siteID & 
                    PhenoYear == cur_sy$PhenoYear,
                    SOS
                ][1]

                # Pheno year
                cur_tp$PhenoYear <- cur_sy$PhenoYear

                return(cur_tp)
            }) %>%
                do.call(rbind, .)

            return(output_dt)
        },

        #' Format data table to a list object that is more convenient for
        #' fitting sping phenology models 
        #' @param dt A data table containing`siteID`, `PhenoYear`, `SOS`, 
        #' `Tmean`, `Tmax`, and `Tmin` 
        #' @return A list object.
        Conv2PhenorLi = function(dt, doy = c((-110):(-1), 1:255)) {
            # get pheno obs
            pheno_dt <- unique(dt[, .(siteID, PhenoYear, SOS)])

            # get daily mean temperature matrix
            temp_mean_mat <- apply(pheno_dt, 1, function(sy) {
                temp <- dt[
                    siteID == sy[["siteID"]] & 
                    PhenoYear == sy[["PhenoYear"]], 
                ]
                if (nrow(temp) != 365) {
                    return(rep(NA, 365))
                }
                return(temp$Tmean)
            })

            # get daily min temperature matrix
            temp_min_mat <- apply(pheno_dt, 1, function(sy) {
                temp <- dt[
                    siteID == sy[["siteID"]] & 
                    PhenoYear == sy[["PhenoYear"]], 
                ]
                if (nrow(temp) != 365) {
                    return(rep(NA, 365))
                }
                return(temp$Tmin)
            })

            # get daily max temperature matrix
            temp_max_mat <- apply(pheno_dt, 1, function(sy) {
                temp <- dt[
                    siteID == sy[["siteID"]] & 
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
        }
        
    ),

    private = list(
        
        # Phenology data in data.table format
        m_pheno_dt = NULL,
        # Phenology data in `phenor` format
        m_pheno_li = NULL

    ),

    active = list(

        get_pheno_dt = function() {
            return(private$m_pheno_dt)
        },
        
        get_pheno_li = function() {
            if (is.null(private$m_pheno_li)) {
                private$m_pheno_li <- self$Conv2PhenorLi(private$m_pheno_dt)
            }
            return(private$m_pheno_li)
        },
        
        set_pheno_li = function(value) {
            stopifnot(!missing(value))
            private$m_pheno_li <- value
            private$m_pheno_dt <- self$Conv2Dt(value)
        }

    )
)


