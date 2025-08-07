#' An internal function that is not supposed to be directly accessed by the user. Find elution time tolerance level.
#'
#' This function finds the time tolerance level. Also, it returns the grouping information given the time tolerance.
#'
#' @param mz mz value of all peaks in all profiles in the study.
#' @param chr retention time of all peaks in all profiles in the study.
#' @param lab label of all peaks in all profiles in the study.
#' @param num.exp The number of spectra in this analysis.
#' @param mz.tol m/z tolerance level for the grouping of signals into peaks. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level.
#' @param chr.tol the elution time tolerance. If NA, the function finds the tolerance level first. If a numerical value is given, the function directly goes to the second step - grouping peaks based on the tolerance.
#' @param aver.bin.size The average bin size to determine the number of equally spaced points in the kernel density estimation.
#' @param min.bins the minimum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too few observations are present.
#' @param max.bins the maximum number of bins to use in the kernel density estimation. It overrides aver.bin.size when too many observations are present.
#' @param max.mz.diff As the m/z tolerance in alignment is expressed in relative terms (ppm), it may not be suitable when the m/z range is wide. This parameter limits the tolerance in absolute terms. It mostly influences feature matching in higher m/z range.
#' @param max.num.segments the maximum number of segments.
#'
#' @details The peaks are first ordered by m/z, and split into groups by the m/z tolerance. Then within every peak group, the pairwise elution time difference is calculated. All the pairwise elution time differences within groups are merged into a single vector. A mixture model (unknown distribution for distance between peaks from the same feature, and a triangle-shaped distribution for distance between peaks from different features) is fit to find the elution time tolerance level. The elution times within each peak group are then ordered. If a gap between consecutive retention times is larger than the elution time tolerance level, the group is further split at the gap. Grouping information is returned, as well as the elution time tolerance level.
#'
#' @return A list object is returned:
#'   \item{chr.tol}{ The elution time tolerance level.}
#'   \item{comp2 }{A matrix with six columns. Every row corrsponds to a peak in one of the spectrum. The columns are: m/z, elution time, spread, signal strength, spectrum label, and peak group label. The rows are ordered by the median m/z of each peak group, and with each peak group the rows are ordered by the elution time.}
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
find.tol.time <- function(mz, chr, lab, num.exp, mz.tol = 2e-5, chr.tol = NA, aver.bin.size = 200, min.bins = 50, max.bins = 100, max.mz.diff = 0.01, max.num.segments = 10000) {
    # 1. Determine order of m/z values
    o <- order(mz)
    # 2. Reorder m/z based on sorted indices
    mz <- mz[o]
    # 3. Reorder retention times accordingly
    chr <- chr[o]
    # 4. Reorder labels accordingly
    lab <- lab[o]
    # 5. Remove temporary variable
    rm(o)

    # 6. Record total number of peaks
    l <- length(mz)

    # 7. Identify breaks between m/z segments using tolerance limits
    breaks <- c(0, which((mz[2:l] - mz[1:(l - 1)]) > min(max.mz.diff, mz.tol * ((mz[2:l] + mz[1:(l - 1)]) / 2))), l)

    # 8. Within each m/z segment, order peaks by retention time
    for (i in 2:length(breaks)) {
        # 9. Order retention times within current segment
        this.o <- order(chr[(breaks[i - 1] + 1):breaks[i]])
        # 10. Adjust indices to global positions
        this.o <- this.o + breaks[i - 1]
        # 11. Reorder m/z values within segment
        mz[(breaks[i - 1] + 1):breaks[i]] <- mz[this.o]
        # 12. Reorder retention times within segment
        chr[(breaks[i - 1] + 1):breaks[i]] <- chr[this.o]
        # 13. Reorder labels within segment
        lab[(breaks[i - 1] + 1):breaks[i]] <- lab[this.o]
    }

    # 14. Remove sentinel values at ends of break vector
    breaks <- breaks[c(-1, -length(breaks))]
    # 15. If retention tolerance not provided, estimate it
    if (is.na(chr.tol)) {
        # 16. Initialize vector of time differences
        da <- 0
        # 17. Limit number of segments processed
        if (length(breaks) > max.num.segments) {
            # 18. Select evenly spaced subset of segments
            s <- floor(seq(2, length(breaks), length.out = max.num.segments))
        } else {
            # 19. Otherwise use all segments
            s <- 2:length(breaks)
        }

        # 20. Iterate through selected segments
        for (i in s) {
            # 21. Indices for peaks in current segment
            this.sel <- (breaks[i - 1] + 1):breaks[i]

            # 22. Consider segments with limited peaks for pairwise differences
            if (length(this.sel) <= 3 * num.exp) {
                # 23. Compute all pairwise retention time differences
                this.d <- as.vector(dist(chr[this.sel]))
                # 24. Sample at most 100 differences for efficiency
                if (length(this.d) > 100) this.d <- sample(this.d, 100)
                # 25. Accumulate differences across segments
                da <- c(da, this.d)
            }
        }

        # 26. Remove missing values from differences
        da <- da[!is.na(da)]
        # 27. Maximum observed difference defines search range
        uppermost <- max(da)
        # 28. Choose density estimation resolution
        n = min(max.bins, max(min.bins, round(length(da) / aver.bin.size)))
        # 29. Estimate density of time differences
        des <- density(da, kernel = "gaussian", n = n, bw = uppermost / n * 2, from = 0)
        # 30. Use positive densities only
        y <- des$y[des$x > 0]
        # 31. Corresponding x-values
        x <- des$x[des$x > 0]

        # 32. Fit linear model to upper tail for baseline estimation
        this.l <- lm(y[x > uppermost / 4] ~ x[x > uppermost / 4])

        # 33. Predict baseline densities across x
        exp.y <- this.l$coef[1] + this.l$coef[2] * x

        # 34. Plot empirical density and baseline
        plot(x, y, main = "find retention time tolerance", xlab = "Delta", ylab = "Density", cex = 0.25)
        # 35. Add baseline line to plot
        lines(x, exp.y, col = "red")
        # 36. Prepare shifted densities for smoothing
        y2 <- y[1:(length(y) - 1)]
        y3 <- y[2:(length(y))]
        # 37. Replace dips with following values to enforce monotonicity
        y2[which(y2 < y3)] <- y3[which(y2 < y3)]
        # 38. Update original density with adjusted values
        y[1:(length(y) - 1)] <- y2

        # 39. Identify where empirical density exceeds baseline
        yy <- cumsum(y > 1.5 * exp.y)
        # 40. Sequence indices for comparison
        yi <- 1:length(yy)
        # 41. Determine cutoff index for tolerance
        sel <- min(which(yy < yi)) - 1

        # 42. Mark tolerance on plot
        abline(v = x[sel], col = "blue")
        # 43. Save estimated retention time tolerance
        chr.tol <- x[sel]
    }

    # 44. Compute consecutive retention time differences
    da <- chr[2:l] - chr[1:(l - 1)]
    # 45. Identify breaks exceeding tolerance
    breaks.2 <- which(da > chr.tol)
    # 46. Combine original breaks with new ones and endpoints
    all.breaks <- c(0, unique(c(breaks, breaks.2)), l)
    # 47. Sort combined breakpoints
    all.breaks <- all.breaks[order(all.breaks)]

    # 48. Initialize group labels for each peak
    grps <- 1:length(mz)
    # 49. Assign group numbers based on breakpoints
    for (i in 2:length(all.breaks)) {
        grps[(all.breaks[i - 1] + 1):all.breaks[i]] <- i
    }

    # 50. Prepare result with reordered data and tolerance
    to.return <- list(mz = mz, chr = chr, lab = lab, grps = grps, chr.tol = chr.tol, mz.tol = mz.tol)
    # 51. Return calculated values
    return(to.return)
}
