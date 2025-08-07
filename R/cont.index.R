#' Continuity index
#'
#' This is an internal function. It uses continuity index (or "run filter") to select putative peaks from EIC.
#'
#' @param newprof The matrix containing m/z, retention time, intensity, and EIC label as columns.
#' @param min.pres Run filter parameter. The minimum proportion of presence in the time period for a series of signals grouped by m/z to be considered a peak.
#' @param min.run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.
#'
#' @details This is the run filter described in Yu et al Bioinformatics 2009.
#'
#' @return A list is returned.
#'   \item{new.rec}{ The matrix containing m/z, retention time, intensity, and EIC label as columns after applying the run filter. }
#'   \item{height.rec}{ The vector of peak heights. }
#'   \item{time.range.rec}{ The vector of peak retention time span. }
#'   \item{mz/pres.rec}{ The vector of proportion of non-missing m/z. }
#'
#' @references Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
cont.index <- function(newprof, min.pres = 0.6, min.run = 5) {
    # 1. Internal helper to collapse rows sharing identical times
    collapse <- function(a) {
        # 2. Sort rows by retention time
        a <- a[order(a[, 2]), ]
        # 3. Total number of rows
        l <- nrow(a)
        # 4. Locate boundaries where time changes
        a.breaks <- which(a[2:l, 2] != a[1:(l - 1), 2])
        # 5. Include first and last positions
        a.breaks <- c(0, a.breaks, l)

        # 6. Placeholder for collapsed results
        newa <- c(0, 0, 0, 0)
        # 7. Loop through each time block
        for (i in 2:length(a.breaks)) {
            # 8. Indices for current block
            sel <- (a.breaks[i - 1] + 1):a.breaks[i]
            # 9. Aggregate by median m/z, representative time, summed intensity and group
            newa <- rbind(newa, c(median(a[a.breaks[i], 1]), a[a.breaks[i], 2], sum(a[sel, 3]), a[a.breaks[i], 4]))
        }

        # 10. Drop placeholder and return collapsed matrix
        return(newa[-1, ])
    }

    # 11. Extract retention-time labels
    labels <- newprof[, 2]
    # 12. Unique time values
    times <- unique(labels)
    # 13. Order time points
    times <- times[order(times)]
    # 14. Count total unique times
    time.points <- length(times)

    # 15. Replace times with index numbers
    for (i in 1:length(times)) labels[which(newprof[, 2] == times[i])] <- i
    # 16. Update profile with indexed times
    newprof[, 2] <- labels

    # 17. Keep ordered time vector
    times <- times[order(times)]
    # 18. Number of observations
    l <- nrow(newprof)
    # 19. Initialize zero timeline for smoothing
    timeline <- rep(0, time.points)
    # 20. Placeholder indices
    i = 1
    num.stack <- 1
    # 21. Convert min.run to count scale
    min.count.run <- min.run * time.points / (max(times) - min(times))
    # 22. Average time distance between points
    aver.time.range <- (max(times) - min(times)) / time.points

    # 23. Extract EIC group labels
    grps <- newprof[, 4]
    # 24. Unique groups
    uniq.grp <- unique(grps)
    # 25. Counter for output records
    curr.label <- 1

    # 26. Count occurrences per group
    ttt <- table(grps)
    # 27. Retain groups meeting minimal presence requirement
    ttt <- ttt[ttt >= max(min.count.run * min.pres, 2)]
    # 28. Convert names back to numeric
    uniq.grp <- as.numeric(names(ttt))

    # 29. Keep only qualifying groups
    newprof <- newprof[newprof[, 4] %in% uniq.grp, ]
    # 30. Order by group then m/z
    newprof <- newprof[order(newprof[, 4], newprof[, 1]), ]
    # 31. Record number of rows after filtering
    r.newprof <- nrow(newprof)
    # 32. Mark boundaries between groups
    breaks <- c(0, which(newprof[1:(r.newprof - 1), 4] != newprof[2:r.newprof, 4]), r.newprof)

    # 33. Placeholder for filtered records
    new.rec <- newprof * 0
    # 34. Pointer to next write position
    rec.pointer <- 1

    # 35. Initialize vectors for peak attributes
    height.rec <- mz.pres.rec <- time.range.rec <- rep(0, length(breaks))
    # 36. Pointer for m/z presence records
    mz.pres.ptr <- 1

    # 37. Round minimum run length to integer
    min.run <- round(min.count.run)

    # 38. Iterate through each group
    for (m in 2:length(breaks)) {
        # 39. Extract current profile subset
        this.prof <- newprof[(breaks[m - 1] + 1):breaks[m], ]

        # 40. Order signals by time
        this.prof <- this.prof[order(this.prof[, 2]), ]
        # 41. Collect time indices
        this.times <- this.prof[, 2]
        # 42. Corresponding intensities
        this.intensi <- this.prof[, 3]
        # 43. Associated m/z values
        this.mass <- this.prof[, 1]
        # 44. Group label
        this.grp <- this.prof[1, 4]

        # 45. Start timeline copy for smoothing
        this.timeline <- timeline
        # 46. Mark present time points
        this.timeline[this.times] <- 1

        # 47. Initialize keep vector
        to.keep <- this.times * 0

        # 48. Smooth presence pattern using box kernel
        dens <- ksmooth(seq(-min.run + 1, length(this.timeline) + min.run), c(rep(0, min.run), this.timeline, rep(0, min.run)), kernel = "box", bandwidth = min.run, x.points = 1:length(this.timeline))
        # 49. Extract smoothed densities
        dens <- dens$y

        # 50. Proceed if any region exceeds presence threshold
        if (max(dens) >= min.pres) {
            # 51. Initialize measured and good-point indicators
            measured.points <- good.points <- timeline
            # 52. Mark measured time points
            measured.points[this.times] <- 1

            # 53. Candidate regions above threshold
            good.sel <- which(dens >= min.pres)
            # 54. Mark good regions
            good.points[good.sel] <- 1
            # 55. Extend good region around each candidate
            for (j in (-min.run):min.run) {
                curr.sel <- good.sel + j
                curr.sel <- curr.sel[curr.sel > 0 & curr.sel <= length(times)]
                good.points[curr.sel] <- 1
            }

            # 56. Retain points measured within good regions
            measured.points <- measured.points * good.points
            # 57. Flag time indices to keep
            to.keep[which(this.times %in% which(measured.points == 1))] <- 1
        }

        # 58. If any points survive filtering
        if (sum(to.keep) > 0) {
            # 59. Select indices marked to keep
            this.sel <- which(to.keep == 1)
            # 60. Assemble filtered profile rows
            this.new <- cbind(this.mass[this.sel], this.times[this.sel], this.intensi[this.sel], rep(this.grp, length(this.sel)))
            # 61. Number of retained rows
            r.new <- nrow(this.new)
            # 62. Record peak height
            height.rec[curr.label] <- max(this.intensi[this.sel])
            # 63. Record time span of peak
            time.range.rec[curr.label] <- times[max(this.times[this.sel])] - times[min(this.times[this.sel])]
            # 64. Record proportion of observed points
            mz.pres.rec[curr.label] <- length(this.sel) / (max(this.times[this.sel]) - min(this.times[this.sel]) + 1)
            # 65. Advance record counter
            curr.label <- curr.label + 1

            # 66. Store filtered rows into output matrix
            new.rec[rec.pointer:(rec.pointer + r.new - 1), ] <- this.new
            # 67. Move pointer for next insertion
            rec.pointer <- rec.pointer + r.new
        }
    }
    # 68. Trim excess preallocated rows
    new.rec <- new.rec[1:(rec.pointer - 1), ]
    # 69. Convert time indices back to actual times
    new.rec[, 2] <- times[new.rec[, 2]]
    # 70. Collect results in list
    results <- new('list')
    # 71. Store filtered records
    results$new.rec <- new.rec
    # 72. Store peak heights
    results$height.rec <- height.rec[1:(curr.label - 1)]
    # 73. Store time ranges
    results$time.range.rec <- time.range.rec[1:(curr.label - 1)]
    # 74. Store m/z presence proportions
    results$mz.pres.rec <- mz.pres.rec[1:(curr.label - 1)]
    # 75. Return run-filter results
    return(results)
}
