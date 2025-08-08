#' Align peaks from spectra into a feature table.
#'
#' Identifies which of the peaks from the profiles correspond to the same feature.
#'
#' @param features A list object. Each component is a matrix which is the output from proc.to.feature().
#' @param min.exp A feature has to show up in at least this number of profiles to be included in the final result.
#' @param mz.tol The m/z tolerance level for peak alignment. The default is NA, which allows the program to search for the tolerance level based on the data. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level.
#' @param chr.tol The retention time tolerance level for peak alignment. The default is NA, which allows the program to search for the tolerance level based on the data.
#' @param find.tol.max.d Argument passed to find.tol(). Consider only m/z diffs smaller than this value.This is only used when the mz.tol is NA.
#' @param max.align.mz.diff As the m/z tolerance is expressed in relative terms (ppm), it may not be suitable when the m/z range is wide. This parameter limits the tolerance in absolute terms. It mostly influences feature matching in higher m/z range.
#' @param nodes The number of CPU cores to be used through parallel processing.
#'
#' @details The function first searches for the m/z tolerance level using a mixture model. After the mz.tol is obtained, the peaks are grouped based on it. Consecutive peaks with m/z value difference smaller than the tolerance level are considered to belong to the same peak group. Non-parametric density estimation within each peak group is used to further split peak groups.
#' The function then searches for the retention time tolerance level. Because the peaks are grouped using m/z, only metabolites that share m/z require this parameter. A rather lenient retention time tolerance level is found using a mixture model. After splitting the peak groups by this value, non-parametric density estimation is used to further split peak groups. Peaks belonging to one group are considered to correspond to the same feature.
#'
#' @return Returns a list object with the following objects in it:
#'   \item{aligned.ftrs}{A matrix, with columns of m/z values, elution times, signal strengths in each spectrum.}
#'   \item{pk.times}{A matrix, with columns of m/z, median elution time, and elution times in each spectrum.}
#'   \item{mz.tol}{The m/z tolerance used in the alignment.}
#'   \item{chr.tol}{The elution time tolerance in the alignment.}
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @seealso proc.to.feature
#'
#' @examples
#' data(features)
#' features.2<-adjust.time(features)
#' this.aligned<-feature.align(features,min.exp=2)
#' summary(this.aligned)
#' this.aligned$aligned.ftrs[1:5,]
#' this.aligned$pk.times[1:5,]
#'
#' @keywords models
# 1. High-level alignment routine: find tolerances, group by m/z and RT, then construct per-feature intensity/time matrices
feature.align <- function(features, min.exp = 2, mz.tol = NA, chr.tol = NA, find.tol.max.d = 1e-4, max.align.mz.diff = 0.01, nodes = 1) {
    # 2. Arrange plotting panels to display progress messages and QC histograms
    par(mfrow = c(3, 2))
    # 3. Blank canvas for on-screen status
    plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
    # 4. Title text indicating current stage
    text(x = 0, y = 0, "Feature alignment", cex = 2)
    # 5. Prepare a second blank plot for upcoming status messages
    plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
    # 6. Helper to fold multiple peaks belonging to the same feature into a single row using either sum or median per experiment
    to.attach <- function(this.pick, num.exp, use = "sum") {

        # 7. Initialize output strengths across experiments
        this.strengths <- rep(0, num.exp)
        if (is.null(nrow(this.pick))) {

            # 8. Single observation: place intensity in its experiment slot; copy mz and chr bounds
            this.strengths[this.pick[6]] <- this.pick[5]
            return(c(this.pick[1], this.pick[2], this.pick[1], this.pick[1], this.strengths))
        } else {
            for (m in 1:length(this.strengths)) {
                cat(m, use)
                # 9. Aggregate intensities within the same experiment either by sum or median
                if (use == "sum") this.strengths[m] <- sum(this.pick[this.pick[, 6] == m, 5])
                if (use == "median") this.strengths[m] <- median(this.pick[this.pick[, 6] == m, 5])
            }
            # 10. Return representative mz/chr and per-experiment intensities
            return(c(mean(this.pick[, 1]), mean(this.pick[, 2]), min(this.pick[, 1]), max(this.pick[, 1]), this.strengths))
        }
    }
    # 11. Determine number of experiments; alignment requires at least 2
    num.exp <- nrow(summary(features))
    if (num.exp > 1) {
        # 12. Flatten the list of feature matrices into concatenated mz/chr vectors with lab indices
        a <- summary(features)

        sizes <- as.numeric(a[, 1]) / ncol(features[[1]])
        sizes <- cumsum(sizes)
        sel <- length(sizes)

        masses <- chr <- lab <- rep(0, sizes[sel])
        sizes <- c(0, sizes)

        for (i in 1:sel) {
            masses[(sizes[i] + 1):sizes[i + 1]] <- features[[i]][, 1]
            chr[(sizes[i] + 1):sizes[i + 1]] <- features[[i]][, 2]
            lab[(sizes[i] + 1):sizes[i + 1]] <- i
        }
        # 13. Sort by m/z (then RT) to enable tolerance-based segmentation
        o <- order(masses, chr)
        masses <- masses[o]
        chr <- chr[o]
        lab <- lab[o]
        l <- length(masses)
        # 14. Find/confirm m/z tolerance: either estimate from data or display provided value
        if (is.na(mz.tol)) {
            mz.tol <- find.tol(masses, uppermost = find.tol.max.d)
            if (length(mz.tol) == 0) {
                mz.tol <- 1e-5
                warning("Automatic tolerance finding failed, 10 ppm was assigned. May need to manually assign alignment mz tolerance level.")
            }
        } else {
            plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", main = "alignment m/z tolerance level given", axes = FALSE)
            text(x = 0, y = 0, mz.tol, cex = 1.2)
        }

        # 15. Estimate/confirm RT tolerance using find.tol.time(); display provided values if given
        if (!is.na(chr.tol)) {
            plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", main = "retention time \n tolerance level given", axes = FALSE)
            text(x = 0, y = 0, chr.tol, cex = 1.2)
        }

        # 16. Call helper to segment peaks by m/z and RT; returns group labels and final chr.tol
        all.ft <- find.tol.time(masses, chr, lab, num.exp = num.exp, mz.tol = mz.tol, chr.tol = chr.tol, max.mz.diff = max.align.mz.diff)
        chr.tol <- all.ft$chr.tol

        message("**** performing feature alignment ****")
        message(paste("m/z tolerance level: ", mz.tol))
        message(paste("time tolerance level:", chr.tol))
        aligned.ftrs <- pk.times <- rep(0, 4 + num.exp)
        mz.sd.rec <- 0

        # 17. Group labels to iterate over candidate features
        labels <- unique(all.ft$grps)

        area <- grps <- masses


        # 18. For each experiment, attach lab/group information; also keep area
        for (i in 1:num.exp) {        
            this <- features[[i]]
            sel <- which(all.ft$lab == i)
            that <- cbind(all.ft$mz[sel], all.ft$chr[sel], all.ft$grps[sel])
            this <- this[order(this[, 1], this[, 2]), ]
            that <- that[order(that[, 1], that[, 2]), ]

            masses[(sizes[i] + 1):sizes[i + 1]] <- this[, 1]
            chr[(sizes[i] + 1):sizes[i + 1]] <- this[, 2]
            area[(sizes[i] + 1):sizes[i + 1]] <- this[, 5]
            grps[(sizes[i] + 1):sizes[i + 1]] <- that[, 3]
            lab[(sizes[i] + 1):sizes[i + 1]] <- i
        }
        # 19. Filter to groups present in >= min.exp experiments; pre-size output
        ttt <- table(all.ft$grps)
        curr.row <- sum(ttt >= min.exp) * 3
        mz.sd.rec <- rep(0, curr.row)
        curr.row <- 1

        sel.labels <- as.numeric(names(ttt)[ttt >= min.exp])

        # 20. Parallelize loop over selected groups (FORK type on Unix); compute per-feature summaries
        cl <- parallel::makeCluster(8, type = "FORK")

        aligned.ftrs <- parallel::parLapply(cl, 1:length(sel.labels), function(i) {

            if (i %% 100 == 0) gc()
            this.return <- NULL
            sel <- which(grps == sel.labels[i])

            if (length(sel) > 1) {
                this <- cbind(masses[sel], chr[sel], chr[sel], chr[sel], area[sel], lab[sel])
                if (length(unique(this[, 6])) >= min.exp) {
                    # 21. Within a group, split by m/z density peaks to avoid merging nearby isotopes/adducts
                    this.den <- density(this[, 1], bw = mz.tol * median(this[, 1]))
                    turns <- find.turn.point(this.den$y)
                    pks <- this.den$x[turns$pks]
                    vlys <- this.den$x[turns$vlys]
                    for (j in 1:length(pks)) {

                        this.lower <- max(vlys[vlys < pks[j]])
                        this.upper <- min(vlys[vlys > pks[j]])
                        this.sel <- which(this[, 1] > this.lower & this[, 1] <= this.upper)

                        that <- this[this.sel, ]

                        if (!is.null(nrow(that))) {
                            if (length(unique(that[, 6])) >= min.exp) {
                                # 22. Further split by RT density peaks; for each sub-cluster attach intensities and times
                                that.den <- density(na.omit(that[, 2]), bw = chr.tol / 1.414)
                                that.turns <- find.turn.point(that.den$y)
                                that.pks <- that.den$x[that.turns$pks]
                                that.vlys <- that.den$x[that.turns$vlys]
                                for (k in 1:length(that.pks)) {
                                    that.lower <- max(that.vlys[that.vlys < that.pks[k]])
                                    that.upper <- min(that.vlys[that.vlys > that.pks[k]])

                                    thee <- that[that[, 2] > that.lower & that[, 2] <= that.upper, ]
                                    if (!is.null(nrow(thee))) {
                                        if (length(unique(thee[, 6])) >= min.exp) {
                                            # 23. Build one row per feature using sum of intensities; and a corresponding row using median RT per experiment
                                            this.return <- c(to.attach(thee, num.exp, use = "sum"), to.attach(thee[, c(1, 2, 3, 4, 2, 6)], num.exp, use = "median"), sd(thee[, 1], na.rm = TRUE))
                                            cat(i, j, k, "surprised if we're here\n") 
                                        }
                                    }
                                }
                            } else {}
                        }
                    }
                } else {}
            } else {
                # 24. Degenerate case: group contains a single observation; only valid when min.exp == 1
                if (min.exp == 1) {
                    thee <- c(masses[sel], chr[sel], chr[sel], chr[sel], area[sel], lab[sel])
                    this.return <- c(to.attach(thee, num.exp, use = "sum"), to.attach(thee[c(1, 2, 3, 4, 2, 6)], num.exp, use = "median"), NA)
                }
            }

            sink()
            rm(list = ls()[ls() != 'this.return'])
            return(this.return)
        })
        parallel::stopCluster(cl)
        
        # 25. Combine parallel results; split into aligned.ftrs (intensities) and pk.times (times)
        aligned.ftrs <- do.call(rbind, aligned.ftrs)
        pk.times <- aligned.ftrs[, (5 + num.exp):(2 * (4 + num.exp))]
        mz.sd.rec <- aligned.ftrs[, ncol(aligned.ftrs)]
        aligned.ftrs <- aligned.ftrs[, 1:(4 + num.exp)]

        colnames(aligned.ftrs) <- c("mz", "chr", "min.mz", "max.mz", paste("exp", 1:num.exp))
        colnames(pk.times) <- c("mz", "chr", "min.mz", "max.mz", paste("exp", 1:num.exp))

        rec <- new("list")
        rec$aligned.ftrs <- aligned.ftrs
        rec$pk.times <- pk.times
        rec$mz.tol <- mz.tol
        rec$chr.tol <- chr.tol

        # 26. QC histograms: m/z dispersion within features; RT dispersion across experiments
        hist(mz.sd.rec, xlab = "m/z SD", ylab = "Frequency", main = "m/z SD distribution")
        hist(apply(pk.times[, -1:-4], 1, sd, na.rm = TRUE), xlab = "Retention time SD", ylab = "Frequency", main = "Retention time SD distribution")
        return(rec)
    } else {
        # 27. Nothing to align when only one experiment is present
        message("There is but one experiment.  What are you trying to align?")
        return(0)
    }
}
