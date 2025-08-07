#' Find peaks and valleys of a curve.
#'
#' This is an internal function which is not supposed to be directly accessed by the user. Finds the peaks and valleys of a smooth curve.
#'
#' @param y The y values of a curve in x-y plane.
#'
#' @details
#'
#' @return A list object:
#'   \item{pks}{The peak positions.}
#'   \item{vlys}{The valley positions}
#'
#' @references Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
find.turn.point <- function(y) {
    # 1. Detect local maxima and minima positions
    peaks2 <- function(x, ties.method) {
        # 2. Embed sequence to evaluate neighbors
        z <- embed(rev(as.vector(c(-Inf, x, -Inf))), dim = 3)
        # 3. Restore original order after embedding
        z <- z[rev(seq(nrow(z))), ]
        # 4. Identify central points greater than neighbors
        v <- max.col(z, ties.method = ties.method) == 2
        # 5. Return logical vector of extrema
        v
    }
    # 6. Wrapper to return maxima and minima indices
    msExtrema <- function(x) {
        # 7. Length of input sequence
        l <- length(x)
        # 8. Detect peaks using first-tie method
        index1 <- peaks2(x, ties.method = "first")
        # 9. Detect valleys using last-tie method on negated data
        index2 <- peaks2(-x, ties.method = "last")
        # 10. Locations that are maxima but not minima
        index.max <- index1 & !index2
        # 11. Locations that are minima but not maxima
        index.min <- index2 & !index1
        # 12. Return list of maxima and minima indices
        list(index.max = index.max, index.min = index.min)
    }

    # 13. Remove missing values from input series
    y <- y[!is.na(y)]
    # 14. If the series is constant, define middle as peak and ends as valleys
    if (length(unique(y)) == 1) {
        # 15. Peak at center of sequence
        pks <- round(length(y) / 2)
        # 16. Valleys at boundaries
        vlys <- c(1, length(y))
        # 17. Create list structure for return
        x <- new("list")
        # 18. Store peak positions
        x$pks <- pks
        # 19. Store valley positions
        x$vlys <- vlys
        # 20. Return result
        return(x)
    }

    # 21. Find extrema for non-constant series
    b <- msExtrema(y)
    # 22. Extract peak indices
    pks <- which(b$index.max)
    # 23. Extract valley indices
    vlys <- which(b$index.min)
    # 24. Include start as valley if first point is not a peak
    if (pks[1] != 1) vlys <- c(1, vlys)
    # 25. Include end as valley if last point is not a peak
    if (pks[length(pks)] != length(y)) vlys <- c(vlys, length(y))

    # 26. Ensure at least two valley points when only one peak exists
    if (length(pks) == 1) vlys <- c(1, length(y))
    # 27. Prepare return list
    x <- new("list")
    # 28. Store peaks
    x$pks <- pks
    # 29. Store valleys
    x$vlys <- vlys
    # 30. Return peak and valley information
    return(x)
}
