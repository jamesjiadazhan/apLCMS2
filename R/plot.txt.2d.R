#' Plot the data in the m/z and retention time plane.
#'
#' This is a diagnostic function. It takes the original text file, as well as the detected feature table, and plots the data in the m/z - retention time plane, using a user-defined range. The entire data is too big to plot, thus the main purpose is to focus on small subregions of the data and check the peak detection results.
#'
#' @param rawname The text file name.
#' @param f The output object of prof.to.feature().
#' @param mzlim The m/z range to plot.
#' @param timelim The retention time range to plot.
#' @param lwd Line width parameter, to be passed on to the function line().
#'
#' @details The columns in the text file need to be separated by tab. The first column needs to be the retention time, the second column the m/z values, and the third column the intensity values. The first row needs to be the column labels, rather than values.
#'
#' @return There is no return value.
#'
#' @references Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
plot.txt.2d <- function(rawname, f, mzlim, timelim, lwd = 1) {
    # 1. Read tab-delimited raw LC/MS data
    ex.nc <- read.table(rawname, header = TRUE, sep = "\t")
    # 2. Convert to matrix for fast indexing
    ex.nc <- as.matrix(ex.nc)
    # 3. Remove zero-intensity points
    ex.nc <- ex.nc[ex.nc[, 3] != 0, ]

    # 4. Filter entries to requested time and m/z ranges
    ex.nc <- ex.nc[ex.nc[, 1] >= timelim[1] & ex.nc[, 1] <= timelim[2] & ex.nc[, 2] >= mzlim[1] & ex.nc[, 2] <= mzlim[2], ]
    # 5. Plot remaining raw data points
    plot(ex.nc[, 1:2], cex = .1)

    # 6. Restrict feature list to m/z window
    f <- f[f[, 1] >= mzlim[1] & f[, 1] <= mzlim[2], ]
    # 7. Add retention-time spans for each detected feature
    for (i in 1:nrow(f)) {
        # 8. Current feature's m/z
        this.mz <- f[i, 1]
        # 9. Compute start and end times using center and spread
        this.time <- c(f[i, 2] - f[i, 3], f[i, 2] + f[i, 3])
        # 10. Draw horizontal line representing feature
        lines(this.time, rep(this.mz, 2), col = "red", lwd = lwd)
    }
}
