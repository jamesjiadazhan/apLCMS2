#' An internal function.
#'
#' This is a internal function. It shouldn't be called by the end user.
#'
#' @param a vector of retention time.
#' @param mz vector of m/z ratio.
#' @param inte vector of signal strength.
#'
#' @references Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
merge.seq.3 <- function(a, mz, inte) {
    # 1. Determine total number of time points
    l <- length(a)
    # 2. Find boundaries where retention time changes
    breaks <- c(0, which(a[1:(l - 1)] != a[2:l]), l)
    # 3. Prepare containers for merged intensities and m/z
    new.int <- new.mz <- rep(0, length(breaks) - 1)

    # 4. Iterate over contiguous time segments
    for (i in 1:(length(breaks) - 1)) {
        # 5. Slice intensity values for this segment
        this.int <- inte[(breaks[i] + 1):breaks[i + 1]]
        # 6. Slice m/z values for this segment
        this.mz <- mz[(breaks[i] + 1):breaks[i + 1]]
        # 7. Sum intensities within the segment
        new.int[i] <- sum(this.int)
        # 8. Record median m/z of peak with maximum intensity
        new.mz[i] <- median(this.mz[which(this.int == max(this.int))])
    }
    # 9. Unique retention times for output
    new.a <- unique(a)
    # 10. Return merged m/z, time, and intensity
    return(cbind(new.mz, new.a, new.int))
}
