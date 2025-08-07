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
    # 1. Prepare logical vector recording if each m/z has a match
    mass.matched.pos <- rep(0, length(x))
    # 2. Loop over every input mass
    for (i in 1:length(x)) {
        # 3. Compute fractional distance to each known m/z value
        this.d <- abs((x[i] - known.mz) / x[i])
        # 4. Flag as matched when smallest distance is within ppm tolerance
        if (min(this.d) < match.tol.ppm / 1e6) mass.matched.pos[i] <- 1
    }
    # 5. Return vector indicating which masses have matches
    return(mass.matched.pos)
}
