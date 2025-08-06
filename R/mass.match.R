#' Internal function: Match m/z values to known metabolites
#'
#' This function matches m/z values to known metabolites within a specified tolerance.
#'
#' @param x Vector of m/z values to be matched.
#' @param known.mz Vector of known m/z values.
#' @param match.tol.ppm The tolerance level in ppm for matching.
#'
#' @return A vector of indicators (0 or 1) indicating whether each m/z value matches to known metabolites.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
mass.match <- function(x, known.mz, match.tol.ppm = 5) {
    mass.matched.pos <- rep(0, length(x))
    for (i in 1:length(x)) {
        this.d <- abs((x[i] - known.mz) / x[i])
        if (min(this.d) < match.tol.ppm / 1e6) mass.matched.pos[i] <- 1
    }
    return(mass.matched.pos)
}
