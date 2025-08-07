#' Removing long ridges at the same m/z.
#'
#' This is an internal function. It substracts a background estimated through kernel smoothing when an EIC continuously span more than half the retention time range.
#'
#' @param x Retetion time vector.
#' @param y2 Intensity vector.
#' @param bw Bandwidth for the kernel smoother. A very wide one is used here.
#'
#' @return A vector of intensity value is returned.
#'
#' @references Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
rm.ridge <- function(x, y2, bw) {
    # 1. Identify low-intensity points as baseline candidates
    sel <- which(y2 < quantile(y2, 0.75))
    # 2. Record extreme retention times within baseline region
    max.x.sel <- max(x[sel])
    min.x.sel <- min(x[sel])

    # 3. Partition retention times relative to baseline range
    in.sel <- which(x >= min.x.sel & x <= max.x.sel)
    over.sel <- which(x > max.x.sel)
    under.sel <- which(x < min.x.sel)

    # 4. Smooth baseline across central region with wide kernel
    this.s <- ksmooth(x[sel], y2[sel], x.points = x[in.sel], kernel = "normal", bandwidth = bw)
    # 5. Abort if smoothing fails to produce values
    if (sum(is.na(this.s$y)) > 0) return(y2)

    # 6. Subtract baseline from central region
    y2[in.sel] <- y2[in.sel] - this.s$y
    # 7. Adjust over-range points using end of baseline
    y2[over.sel] <- y2[over.sel] - this.s$y[which(this.s$x == max(this.s$x))[1]]
    # 8. Adjust under-range points using start of baseline
    y2[under.sel] <- y2[under.sel] - this.s$y[which(this.s$x == min(this.s$x))[1]]

    # 9. Force negative intensities to zero
    y2[y2 < 0] <- 0
    # 10. Return ridge-removed intensities
    return(y2)
}
