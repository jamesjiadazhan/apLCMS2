#' Adaptive binning specifically for the machine learning approach.
#'
#' This is an internal function. It creates EICs using adaptive binning procedure
#'
#' @param x A matrix with columns of m/z, retention time, intensity.
#' @param tol m/z tolerance level for the grouping of data points. This value is expressed as the fraction of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level. The recommended value is the machine's nominal accuracy level. Divide the ppm value by 1e6. For FTMS, 1e-5 is recommended.
#' @param ridge.smoother.window The size of the smoother window used by the kernel smoother to remove long ridge noise from the EIC.
#' @param baseline.correct After grouping the observations, the highest intensity in each group is found. If the highest is lower than this value, the entire group will be deleted. The default value is NA, in which case the program uses the 75th percentile of the height of the noise groups.
#'
#' @details It uses repeated smoothing and splitting to separate EICs. The details are described in the reference and flowchart.
#'
#' @return A list is returned.
#'   \item{height.rec}{The records of the height of each EIC.}
#'   \item{masses}{ The vector of m/z values after binning. }
#'   \item{labels}{ The vector of retention time after binning. }
#'   \item{intensi}{ The vector of intensity values after binning. }
#'   \item{grps}{ The EIC labels, i.e. which EIC each observed data point belongs to. }
#'   \item{times}{ All the unique retention time values, ordered. }
#'   \item{tol}{ The m/z tolerance level. }
#'
#' @references Bioinformatics. 30(20): 2941-2948. Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
adaptive.bin.2 <- function(x, tol, ridge.smoother.window = 50, baseline.correct) {
    # 1. Pull out mass, time, and intensity vectors
    masses <- x$masses
    labels <- x$labels
    intensi <- x$intensi
    times <- x$times

    # 2. Remove original object to free memory
    rm(x)

    # 3. Prepare base curve of unique times
    base.curve <- unique(times)
    base.curve <- base.curve[order(base.curve)]
    base.curve <- cbind(base.curve, base.curve * 0)

    # 4. Sort observations by mass
    curr.order <- order(masses)
    intensi <- intensi[curr.order]
    labels <- labels[curr.order]
    masses <- masses[curr.order]

    rm(curr.order)

    # 5. Estimate global mass density and locate valleys
    cat(c("m/z tolerance is: ", tol, "\n"))
    l <- length(masses)
    curr.bw <- 0.5 * tol * max(masses)
    all.mass.den <- density(masses, weights = intensi / sum(intensi), bw = curr.bw, n = 2^min(15, floor(log2(l)) - 2))
    all.mass.turns <- find.turn.point(all.mass.den$y)
    all.mass.vlys <- all.mass.den$x[all.mass.turns$vlys]
    breaks <- c(0, unique(round(approx(masses, 1:l, xout = all.mass.vlys, rule = 2, ties = 'ordered')$y))[-1])

    # 6. Initialize bookkeeping
    grps <- masses * 0
    times <- unique(labels)
    times <- times[order(times)]
    curr.label <- 1
    time.range <- diff(range(times))
    aver.time.range <- (max(labels) - min(labels)) / length(times)
    if (ridge.smoother.window < time.range / 10) ridge.smoother.window <- time.range / 10

    newprof <- matrix(0, nrow = length(masses), ncol = 4)
    height.rec <- matrix(0, nrow = length(masses), ncol = 8)
    colnames(height.rec) <- c("group.label", "time.points", "maximum.height", "normalized.mz.range", "normalized.mz.sd.weighted", "normalized.mz.range.before.merge", "normalized.mz.sd.weighted.before.merge", "RT.peak.location")
    prof.pointer <- 1
    height.pointer <- 1

    # 7. Iterate through mass-density valleys
    for (i in 1:(length(breaks) - 1)) {
        this.labels <- labels[(breaks[i] + 1):breaks[i + 1]]
        this.masses <- masses[(breaks[i] + 1):breaks[i + 1]]
        this.intensi <- intensi[(breaks[i] + 1):breaks[i + 1]]

        curr.order <- order(this.labels)
        this.masses <- this.masses[curr.order]
        this.intensi <- this.intensi[curr.order]
        this.labels <- this.labels[curr.order]

        this.bw = 0.5 * tol * median(this.masses)
        mass.den <- density(this.masses, weights = this.intensi / sum(this.intensi), bw = this.bw)
        mass.den$y[mass.den$y < min(this.intensi) / 10] <- 0
        mass.turns <- find.turn.point(mass.den$y)
        mass.pks <- mass.den$x[mass.turns$pks]
        mass.vlys <- c(-Inf, mass.den$x[mass.turns$vlys], Inf)

        # 8. Process each peak within the slice
        for (j in 1:length(mass.pks)) {
            mass.lower <- max(mass.vlys[mass.vlys < mass.pks[j]])
            mass.upper <- min(mass.vlys[mass.vlys > mass.pks[j]])

            if (length(mass.pks) == 1) mass.lower <- mass.lower - 1
            mass.sel <- which(this.masses > mass.lower & this.masses <= mass.upper)

            if (length(mass.sel) > 0) {
                that.labels <- this.labels[mass.sel]
                that.masses <- this.masses[mass.sel]
                that.intensi <- this.intensi[mass.sel]

                # 9. Record pre-merge mass spread statistics
                that.masses.mean <- weighted.mean(that.masses, that.intensi)
                mz.range.before.merge <- abs(diff(range(that.masses))) / that.masses.mean
                mz.sd.before.merge <- sqrt(weighted.mean((that.masses - that.masses.mean)^2, that.intensi)) / that.masses.mean

                # 10. Merge contiguous time segments
                that.merged <- merge.seq.3(that.labels, that.masses, that.intensi)
                if (nrow(that.merged) == 1) {
                    new.merged <- that.merged
                } else {
                    new.merged <- that.merged[order(that.merged[, 1]), ]
                }

                that.labels <- new.merged[, 2]
                that.masses <- new.merged[, 1]
                that.intensi <- new.merged[, 3]
                that.range <- diff(range(that.labels))
                that.RT.peak.loc <- that.labels[which(that.intensi == max(that.intensi))[1]]

                # 11. Remove ridge background if peak spans long time
                if (that.range > 0.5 * time.range & length(that.labels) > that.range / aver.time.range * 0.6) {
                    that.intensi <- rm.ridge(that.labels, that.intensi, bw = max(ridge.smoother.window, that.range / 2))

                    that.sel <- which(that.intensi != 0)
                    that.labels <- that.labels[that.sel]
                    that.masses <- that.masses[that.sel]
                    that.intensi <- that.intensi[that.sel]
                }

                # 12. Store cleaned segment and statistics
                that.n <- length(that.masses)
                newprof[prof.pointer:(prof.pointer + that.n - 1), ] <- cbind(that.masses, that.labels, that.intensi, rep(curr.label, that.n))
                prof.pointer <- prof.pointer + that.n
                that.masses.mean <- weighted.mean(that.masses, that.intensi)
                height.rec[height.pointer, ] <- c(curr.label, that.n, max(that.intensi), abs(diff(range(that.masses))) / that.masses.mean, sqrt(weighted.mean((that.masses - that.masses.mean)^2, that.intensi)) / that.masses.mean, mz.range.before.merge, mz.sd.before.merge, that.RT.peak.loc)
                height.pointer <- height.pointer + 1
                curr.label <- curr.label + 1
            }
        }
    }

    # 13. Trim matrices and order results
    newprof <- newprof[1:(prof.pointer - 1), ]
    height.rec <- height.rec[1:(height.pointer - 1), ]
    newprof <- newprof[order(newprof[, 1], newprof[, 2]), ]

    # 14. Assemble and return profile list
    raw.prof <- new("list")
    raw.prof$height.rec <- height.rec
    raw.prof$masses <- newprof[, 1]
    raw.prof$labels <- newprof[, 2]
    raw.prof$intensi <- newprof[, 3]
    raw.prof$grps <- newprof[, 4]
    raw.prof$times <- times
    raw.prof$tol <- tol

    return(raw.prof)
}
