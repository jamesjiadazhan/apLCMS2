#' Adaptive binning
#'
#' This is an internal function. It creates EICs using adaptive binning procedure
#'
#' @param x A matrix with columns of m/z, retention time, intensity.
#' @param min.run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.
#' @param min.pres Run filter parameter. The minimum proportion of presence in the time period for a series of signals grouped by m/z to be considered a peak.
#' @param tol m/z tolerance level for the grouping of data points. This value is expressed as the fraction of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level. The recommended value is the machine's nominal accuracy level. Divide the ppm value by 1e6. For FTMS, 1e-5 is recommended.
#' @param baseline.correct After grouping the observations, the highest intensity in each group is found. If the highest is lower than this value, the entire group will be deleted. The default value is NA, in which case the program uses the 75th percentile of the height of the noise groups.
#' @param weighted Whether to weight the local density by signal intensities.
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
#'   \item{min.count.run}{ The minimum number of elution time points for a series of signals grouped by m/z to be considered a peak.}
#'   \item{weighted}{Whether to weight the local density by signal intensities.}
#'
#' @references Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
adaptive.bin <- function(x, min.run, min.pres, tol, baseline.correct, weighted = FALSE) {
    # 1. Extract input vectors for mass, time, and intensity
    masses <- x$masses
    labels <- x$labels
    intensi <- x$intensi
    times <- x$times

    # 2. Drop original list to free memory
    rm(x)

    # 3. Prepare base curve of unique times for later reference
    base.curve <- unique(times)
    base.curve <- base.curve[order(base.curve)]
    base.curve <- cbind(base.curve, base.curve * 0)

    # 4. Order observations by mass
    curr.order <- order(masses)
    intensi <- intensi[curr.order]
    labels <- labels[curr.order]
    masses <- masses[curr.order]

    rm(curr.order)

    # 5. Report tolerance and set up mass density estimate
    cat(c("m/z tolerance is: ", tol, "\n"))
    l <- length(masses)
    curr.bw <- 0.5 * tol * max(masses)
    if (weighted) {
        all.mass.den <- density(masses, weights = intensi / sum(intensi), bw = curr.bw, n = 2^min(15, floor(log2(l)) - 2))
    } else {
        all.mass.den <- density(masses, bw = curr.bw, n = 2^min(15, floor(log2(l)) - 2))
    }
    all.mass.turns <- find.turn.point(all.mass.den$y)
    all.mass.vlys <- all.mass.den$x[all.mass.turns$vlys]
    breaks <- c(0, unique(round(approx(masses, 1:l, xout = all.mass.vlys, rule = 2, ties = 'ordered')$y))[-1])

    # 6. Initialize bookkeeping structures
    grps <- masses * 0 # this is which peak group the signal belongs to, not time group
    min.lab <- min(labels)
    max.lab <- max(labels)
    times <- unique(labels)
    times <- times[order(times)]
    curr.label <- 1
    min.count.run <- min.run * length(times) / (max(times) - min(times))
    time.range <- diff(range(times))
    aver.time.range <- (max(labels) - min(labels)) / length(times)

    newprof <- matrix(0, nrow = length(masses), ncol = 4)
    height.rec <- matrix(0, nrow = length(masses), ncol = 3)
    prof.pointer <- 1
    height.pointer <- 1

    # 7. Iterate over mass-density valleys to isolate EIC candidates
    for (i in 1:(length(breaks) - 1)) {
        this.labels <- labels[(breaks[i] + 1):breaks[i + 1]]
        if (length(unique(this.labels)) >= min.count.run * min.pres) {
            # 8. Extract masses and intensities for current slice
            this.masses <- masses[(breaks[i] + 1):breaks[i + 1]]
            this.intensi <- intensi[(breaks[i] + 1):breaks[i + 1]]

            curr.order <- order(this.labels)
            this.masses <- this.masses[curr.order]
            this.intensi <- this.intensi[curr.order]
            this.labels <- this.labels[curr.order]

            # 9. Find peaks within this mass slice
            this.bw = 0.5 * tol * median(this.masses)
            if (weighted) {
                mass.den <- density(this.masses, weights = this.intensi / sum(this.intensi), bw = this.bw)
            } else {
                mass.den <- density(this.masses, bw = this.bw)
            }
            mass.den$y[mass.den$y < min(this.intensi) / 10] <- 0
            mass.turns <- find.turn.point(mass.den$y)
            mass.pks <- mass.den$x[mass.turns$pks]
            mass.vlys <- c(-Inf, mass.den$x[mass.turns$vlys], Inf)

            # 10. Process each detected mass peak
            for (j in 1:length(mass.pks)) {
                mass.lower <- max(mass.vlys[mass.vlys < mass.pks[j]])
                mass.upper <- min(mass.vlys[mass.vlys > mass.pks[j]])

                if (length(mass.pks) == 1) mass.lower <- mass.lower - 1
                mass.sel <- which(this.masses > mass.lower & this.masses <= mass.upper)

                if (length(mass.sel) > 0) {
                    that.labels <- this.labels[mass.sel]
                    that.masses <- this.masses[mass.sel]
                    that.intensi <- this.intensi[mass.sel]

                    # 11. Merge contiguous time sequences
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

                    # 12. Remove long baseline ridges when necessary
                    if (that.range > 0.5 * time.range & length(that.labels) > that.range * min.pres & length(that.labels) / (diff(range(that.labels)) / aver.time.range) > min.pres) {
                        that.intensi <- rm.ridge(that.labels, that.intensi, bw = max(10 * min.run, that.range / 2))

                        that.sel <- which(that.intensi != 0)
                        that.labels <- that.labels[that.sel]
                        that.masses <- that.masses[that.sel]
                        that.intensi <- that.intensi[that.sel]
                    }
                    # 13. Record cleaned profile segment
                    that.n <- length(that.masses)
                    newprof[prof.pointer:(prof.pointer + that.n - 1), ] <- cbind(that.masses, that.labels, that.intensi, rep(curr.label, that.n))
                    prof.pointer <- prof.pointer + that.n
                    height.rec[height.pointer, ] <- c(curr.label, that.n, max(that.intensi))
                    height.pointer <- height.pointer + 1
                    curr.label <- curr.label + 1
                }
            }
        } else {
            # 14. Randomly keep small segments to maintain background estimate
            if (runif(1) < 0.05) {
                this.masses <- masses[(breaks[i] + 1):breaks[i + 1]]
                this.intensi <- intensi[(breaks[i] + 1):breaks[i + 1]]
                curr.order <- order(this.labels)
                this.masses <- this.masses[curr.order]
                this.intensi <- this.intensi[curr.order]
                this.labels <- this.labels[curr.order]

                that.merged <- merge.seq.3(this.labels, this.masses, this.intensi)
                that.n <- nrow(that.merged)
                newprof[prof.pointer:(prof.pointer + that.n - 1), ] <- cbind(that.merged, rep(curr.label, that.n))
                prof.pointer <- prof.pointer + that.n

                height.rec[height.pointer, ] <- c(curr.label, that.n, max(that.merged[, 3]))
                height.pointer <- height.pointer + 1
                curr.label <- curr.label + 1
            }
        }
    }

    # 15. Trim unused rows and order by m/z and time
    newprof <- newprof[1:(prof.pointer - 1), ]
    height.rec <- height.rec[1:(height.pointer - 1), ]
    newprof <- newprof[order(newprof[, 1], newprof[, 2]), ]

    # 16. Assemble output list
    raw.prof <- new("list")
    raw.prof$height.rec <- height.rec
    raw.prof$masses <- newprof[, 1]
    raw.prof$labels <- newprof[, 2]
    raw.prof$intensi <- newprof[, 3]
    raw.prof$grps <- newprof[, 4]
    raw.prof$times <- times
    raw.prof$tol <- tol
    raw.prof$min.count.run <- min.count.run

    # 17. Return profile summary
    return(raw.prof)
}
