#' An internal function that is not supposed to be directly accessed by the user. Find m/z tolerance level.
#'
#' The function finds the tolerance level in m/z from a given vector of observed m/z values.
#'
#' @param a The vector of observed m/z values.
#' @param uppermost Consider only m/z diffs smaller than this value.
#' @param aver.bin.size The average bin size to determine the number of equally spaced points in the kernel density estimation.
#' @param min.bins The minimum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too few observations are present.
#' @param max.bins The maximum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too many observations are present.
#'
#' @details The method assumes a mixture model: an unknown distribution of m/z variations in the same peak, and an exponential distribution of between-peak diffs. The parameter of the exponential distribution is estimated by the upper 75% of the sorted data, and the cutoff is selected where the density of the empirical distribution is >1.5 times the density of the exponential distribution.
#'
#' @return The tolerance level is returned.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @examples
#' data(prof)
#' find.tol(prof[[1]][,1])
#'
#' @keywords models
find.tol <- function(a, uppermost = 1e-4, aver.bin.size = 4000, min.bins = 50, max.bins = 200) {
    # 1. Sort m/z values to compute adjacent differences
    a <- a[order(a)]
    # 2. Determine number of observed m/z values
    l <- length(a)
    # 3. Compute relative differences between neighbors
    da <- (a[2:l] - a[1:(l - 1)]) / ((a[2:l] + a[1:(l - 1)]) / 2)
    # 4. Keep only small differences below threshold
    da <- da[da < uppermost]
    # 5. Choose kernel density resolution
    n <- min(max.bins, max(round(length(da) / aver.bin.size), min.bins))
    # 6. Estimate density of differences
    des <- density(da, kernel = "gaussian", n = n, bw = uppermost / n * 2, from = 0)
    # 7. Ignore zero to avoid singularities
    y <- des$y[des$x > 0]
    # 8. Corresponding x-values for positive densities
    x <- des$x[des$x > 0]

    # 9. Use tail of distribution for exponential fit
    to.use <- da[da > max(da) / 4] - max(da) / 4
    # 10. Load MASS for maximum likelihood fitting
    library(MASS)
    # 11. Estimate rate of exponential component
    this.rate <- fitdistr(to.use, "exponential")$estimate
    # 12. Compute exponential density
    exp.y <- dexp(x, rate = this.rate)
    # 13. Scale exponential to match tail area
    exp.y <- exp.y * sum(y[x > max(da) / 4]) / sum(exp.y[x > max(da) / 4])

    # 14. Plot empirical density
    plot(x, y, xlab = "Delta", ylab = "Density", main = "find m/z tolerance", cex = .25)
    # 15. Overlay fitted exponential for comparison
    lines(x, exp.y, col = "red")

    # 16. Cumulative density of empirical distribution
    yy <- cumsum(y)
    # 17. Cumulative density of exponential distribution (unused but kept for clarity)
    y2 <- cumsum(exp.y)
    # 18. Identify where empirical density exceeds 1.5x exponential
    yy <- cumsum(y > 1.5 * exp.y)
    # 19. Sequence index for threshold detection
    yi <- 1:length(yy)
    # 20. Determine last index meeting criterion
    sel <- min(which(yy < yi)) - 1

    # 21. Mark tolerance cutoff on plot
    abline(v = x[sel], col = "blue")
    # 22. Return estimated m/z tolerance level
    return(x[sel])
}
