#' Interpolate missing intensities and calculate the area for a single EIC.
#'
#' This is an internal function that's not supposed to be called directly by the user.
#'
#' @param x The positions of x(retention time) where non-NA y is observed.
#' @param y The observed intensities.
#' @param all.x All possible x(retention time) in the LCMS profile.
#' @param all.w The "footprint" of each measured retention time, used as weight for the corresponding y.
#'
#' @details This is an internal function. It interpolates missing y using linear interpolation, and then calculates the area under the curve.
#'
#' @return The area is returned.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
# x is retention time, all.w is weight of each all.x, y is observed intensities
interpol.area <- function(x, y, all.x, all.w) {
    r <- range(x)
    x.sel <- which(all.x >= r[1] & all.x <= r[2])

    x2 <- all.x[x.sel]
    w <- all.w[x.sel]

    y2 <- approx(x, y, xout = x2, method = "linear")$y
    return(sum(w * y2))
}
