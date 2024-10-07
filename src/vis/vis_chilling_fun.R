#'******************************************************************************
#' Description: Visualize chilling functions.
#'******************************************************************************
source("src/hlp/hlp_chilling_calc.R")

# Visualize Utal model
plot(NA, 
    xlim = c(0, 20), xlab = expression(Temperature~(degree~C)),
    xaxt = "n",
    ylim = c(-1, 1), ylab = "Chilling Units",
    bty = "l"
)
axis(side = 1, at = c(1.4, 2.4, 9.1, 12.4, 15.9, 18))
title("Utah chilling model")

abline(h = 0, lty = 2)
segments(0, 0, 1.4, 0, lwd = 2)
segments(1.4, 0.5, 2.4, 0.5, lwd = 2)
segments(2.4, 1, 9.1, 1, lwd = 2)
segments(9.1, 0.5, 12.4, 0.5, lwd = 2)
segments(12.4, 0, 15.9, 0, lwd = 2)
segments(15.9, -0.5, 18, -0.5, lwd = 2)
segments(18, -1, 25, -1, lwd = 2)





