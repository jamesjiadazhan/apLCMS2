#' Plot the data in the m/z and retention time plane.
#'
#' This is a diagnostic function. It takes the original CDF file, as well as the detected feature table, and plots the data in the m/z - retention time plane, using a user-defined range. The entire data is too big to plot, thus the main purpose is to focus on small subregions of the data and check the peak detection results.
#'
#' @param rawname The CDF file name.
#' @param f The output object of prof.to.feature().
#' @param mzlim The m/z range to plot.
#' @param timelim The retention time range to plot.
#' @param lwd Line width parameter, to be passed on to the function line().
#'
#' @return There is no return value.
#'
#' @references Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
plot.cdf.2d <- function(rawname, f, mzlim, timelim, lwd = 1) {
    # 1. Load raw LC/MS data from CDF file
    this.raw <- load.lcms(rawname)
    # 2. Extract m/z values
    masses <- this.raw$masses
    # 3. Extract intensity values
    intensi <- this.raw$intensi
    # 4. Extract retention time labels
    labels <- this.raw$labels
    # 5. Extract actual times
    times <- this.raw$times
    # 6. Release raw object to save memory
    rm(this.raw)

    # 7. Sort time values for plotting
    times <- times[order(times)]
    # 8. Generate base curve of unique times
    base.curve <- unique(times)
    # 9. Ensure base curve is ordered
    base.curve <- base.curve[order(base.curve)]
    # 10. Create placeholder second column for plotting baseline
    base.curve <- cbind(base.curve, base.curve * 0)

    # 11. Order all measurements by m/z for consistency
    curr.order <- order(masses)
    intensi <- intensi[curr.order]
    labels <- labels[curr.order]
    masses <- masses[curr.order]

    # 12. Select points within requested time and m/z ranges
    sel <- which(labels >= timelim[1] & labels <= timelim[2] & masses >= mzlim[1] & masses <= mzlim[2])
    # 13. Build matrix of selected time and m/z values
    ex.nc <- cbind(labels[sel], masses[sel])
    # 14. Scatter plot of raw points in selected region
    plot(ex.nc[, 1:2], cex = .1)

    # 15. Keep only features within m/z window
    f <- f[f[, 1] >= mzlim[1] & f[, 1] <= mzlim[2], ]
    # 16. Draw horizontal lines for each feature's retention window
    for (i in 1:nrow(f)) {
        # 17. m/z value for current feature
        this.mz <- f[i, 1]
        # 18. Compute start and end times using center and spread
        this.time <- c(f[i, 2] - f[i, 3], f[i, 2] + f[i, 3])
        # 19. Add line representing detected feature
        lines(this.time, rep(this.mz, 2), col = "red", lwd = lwd)
    }
}
