#'******************************************************************************
#' Description: Extract temperature seasonality metrics from temperature time
#' series.
#'******************************************************************************
require(data.table)
require(lubridate)




#' For the point time sereis, compute average temperature time series across years,
#' then smooth it and extract the corresponding seasonal metrics.
#' 
#' @param dt Temperature time series data table.
#' @param tc_th,t0_th Thresholds ranging from [0,1] to extract critical dates.
TempSeason <- function(dt, tc_th, t0_th) {
    # Calculate DOY
    doy <- ifelse(year(dt$Date) < dt$PhenoYear,
        dt$Date - as.Date(paste0(dt$PhenoYear, "-01-01")),
        dt$Date - as.Date(paste0(dt$PhenoYear, "-01-01")) + 1
    )
    dt$doy <- doy

    # Avg temperature time series
    avg_temp_dt <- dt[, mean(Tmean, na.rm = TRUE), by = doy]
    # Smooth temperature
    avg_temp_dt$smooth <- filter(avg_temp_dt$V1, 
        filter = rep(1 / 30, 30), 
        circular = TRUE
    )

    # Minimum
    min_idx <- which.min(avg_temp_dt$smooth)
    # Autumn temperature
    autumn_temp <- avg_temp_dt$smooth[1:min_idx]
    # Autumn max and min
    autumn_max <- max(autumn_temp, na.rm = TRUE)
    autumn_min <- min(autumn_temp, na.rm = TRUE)
    # Autumn amplitude
    autumn_amp <- autumn_max - autumn_min

    # Chilling start
    t0chill <- avg_temp_dt[which.min(
        abs(autumn_min + autumn_amp * tc_th - autumn_temp)
    ), doy]


    # Spring temperature
    spring_temp <- avg_temp_dt$smooth[min_idx:nrow(avg_temp_dt)]
    # Spring max, min, and amplitude
    spring_max <- max(spring_temp, na.rm = TRUE)
    spring_min <- min(spring_temp, na.rm = TRUE)
    # Spring amplitude
    spring_amp <- spring_max - spring_min

    # Forcing start
    t0forc <- avg_temp_dt[min_idx + which.min(
        abs(spring_min + spring_amp * t0_th - spring_temp)
    ), doy]


    res_row <- data.table(
        siteID = dt$siteID[1],
        t0chill = t0chill,
        autumn_max = autumn_max,
        autumn_min = autumn_min,
        autumn_amp = autumn_amp,
        t0forc = t0forc,
        spring_max = spring_max,
        spring_min = spring_min,
        spring_amp = spring_amp
    )

    return(res_row)
}



#' Plot illustratiof temperature seasonality
#' 
#' @param dt Temperature time series data table.
#' @param tc_th,t0_th Thresholds ranging from [0,1] to extract critical dates.
PlotTempSeason <- function(dt, tc_th, t0_th) {
    temp_sea <- TempSeason(dt, tc_th, t0_th)

    # Plot unsmoothed and smoothed temperatures
    # Calculate DOY
    doy <- ifelse(year(dt$Date) < dt$PhenoYear,
        dt$Date - as.Date(paste0(dt$PhenoYear, "-01-01")),
        dt$Date - as.Date(paste0(dt$PhenoYear, "-01-01")) + 1
    )
    dt$doy <- doy

    # Unsmoothed
    rg <- range(dt$Tmean)
    plot(NA,
        xlim = c(-110, 255), ylim = rg,
        bty = "L", xlab = "DOY", 
        ylab = expression("Temperature" ~ (degree ~ C)),
        mgp = c(1.5, 0.5, 0)
    )
    null <- by(dt, dt$PhenoYear, function(yr_dt) {
        lines(yr_dt[, .(doy, Tmean)], type = "l", col = "grey")
    })

    # Smoothed
    # Avg temperature time series
    avg_temp_dt <- dt[, mean(Tmean, na.rm = TRUE), by = doy]
    avg_temp_dt$smooth <- filter(avg_temp_dt$V1,
        filter = rep(1 / 30, 30),
        circular = TRUE
    )
    lines(avg_temp_dt[, .(doy, smooth)], lwd = 2)

    # t0chill
    segments(temp_sea$t0chill, grconvertY(0, "npc", "user"),
        temp_sea$t0chill, avg_temp_dt[doy == temp_sea$t0chill, smooth],
        lty = 2
    )
    text(temp_sea$t0chill, grconvertY(0, "npc", "user"),
        adj = c(1.1, -1),
        labels = expression(t[c])
    )

    # t0
    segments(temp_sea$t0forc, grconvertY(0, "npc", "user"),
        temp_sea$t0forc, avg_temp_dt[doy == temp_sea$t0forc, smooth],
        lty = 2
    )
    text(temp_sea$t0forc, grconvertY(0, "npc", "user"),
        adj = c(1.1, -1),
        labels = expression(t[0])
    )

    # T_base
    T_base <- max(0, as.vector(avg_temp_dt[doy == temp_sea$t0forc, smooth]))
    segments(temp_sea$t0forc, T_base,
        grconvertX(0.9, "npc", "user"), T_base,
        lty = 2
    )
    text(temp_sea$t0forc, T_base,
        adj = c(-5, -1),
        labels = expression(T["base"])
    )
}
