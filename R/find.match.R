#' Internal function: finding the best match between a set of detected features and a set of known features.
#'
#' Given a small matrix of distances, find the best column-row pairing that minimize the sum of distances of the matched pairs.
#'
#' @param a A matrix of distances.
#' @param unacceptable A distance larger than which cannot be accepted as pairs.
#'
#' @details
#'
#' @return A matrix the same dimension as the input matrix, with matched position taking value 1, and all other positions taking value 0.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
find.match <- function(a, unacceptable = 4) {
    # 1. Helper to locate smallest value in matrix
    find.min.pos <- function(d) {
        # 2. Linear index of minimum distance
        pos <- which(d == min(d))[1]
        # 3. Convert to row index
        pos.x <- pos %% nrow(d)
        # 4. Correct when modulus is zero
        if (pos.x == 0) pos.x <- nrow(d)
        # 5. Compute column index
        pos.y <- floor((pos - 1) / nrow(d) + 1)
        # 6. Return indices as vector
        pos <- c(pos.x, pos.y)
        # 7. Provide position of minimum distance
        return(pos)
    }

    # 8. Initialize assignment matrix with zeros
    b <- a * 0
    # 9. Handle case with single known feature
    if (ncol(a) == 1) {
        # 10. Locate smallest distance
        sel <- which(a[, 1] == min(a[, 1]))[1]
        # 11. Mark match if within tolerance
        if (a[sel, 1] <= unacceptable) b[sel, 1] <- 1
    # 12. Handle case with single detected feature
    } else if (nrow(a) == 1) {
        # 13. Locate smallest distance
        sel <- which(a[1, ] == min(a[1, ]))[1]
        # 14. Mark match if within tolerance
        if (a[1, sel] <= unacceptable) b[1, sel] <- 1
    } else {
        # 15. Find closest pair among all combinations
        p <- find.min.pos(a)
        # 16. Repeat until all remaining distances exceed tolerance
        while (a[p[1], p[2]] <= unacceptable) {
            # 17. Record matched pair
            b[p[1], p[2]] <- 1
            # 18. Exclude matched row from further consideration
            a[p[1], ] <- 1e10
            # 19. Exclude matched column from further consideration
            a[, p[2]] <- 1e10
            # 20. Find next closest pair
            p <- find.min.pos(a)
        }
    }
    # 21. Return matrix indicating optimal matches
    return(b)
}
