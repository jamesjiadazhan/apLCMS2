#' Internal function: Extract data feature from EIC.
#'
#' The function extracts data features after applying different smoother settings.
#'
#' @param raw.prof The data after adaptive binning, i.e. the output from adaptive.bin.2().
#' @param smoother.window The smoother window sizes to use for data feature extraction.
#'
#' @details We take a number of data characteristic measurements from each EIC, including m/z span, m/z standard deviation, retention time (RT) span, RT peak location, and summary statistics on the raw intensity values of the EIC. We also centroid the data in each EIC such that it becomes two-dimensional data (intensity v.s. RT). We then apply different smoothers (shape/window size) in combination of different weighting schemes (unweighted, weighted with intensity, weighted with log intensity) to each EIC. At each smoothing setting, we record summary statistics of smoothed data.
#'
#' @return A matrix. Every row corresponds to an EIC. Every column corresponds to a data feature.
#'
#' @references Bioinformatics. 30(20): 2941-2948.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
eic.disect <- function(raw.prof, smoother.window = c(1, 5, 10)) {
    # 1. Combine mass, time, intensity, and group into matrix
    newprof <- cbind(raw.prof$masses, raw.prof$labels, raw.prof$intensi, raw.prof$grps)
    # 2. Preserve precomputed height features for later merge
    height.rec <- raw.prof$height.rec

    # 3. Extract time labels for relabeling
    labels <- newprof[, 2]
    # 4. Identify unique times present in profile
    times <- unique(labels)
    # 5. Sort unique times
    times <- times[order(times)]
    # 6. Count number of time points
    time.points <- length(times)

    # 7. Replace original labels with sequential indices
    for (i in 1:length(times)) labels[which(newprof[, 2] == times[i])] <- i # now labels is the index of time points
    # 8. Update profile with indexed time labels
    newprof[, 2] <- labels

    # 9. Re-sort times (ensures order after indexing)
    times <- times[order(times)]
    # 10. Determine number of observations
    l <- nrow(newprof)
    # 11. Initialize timeline placeholder
    timeline <- rep(0, time.points)
    # 12. Average time range per point for reference
    aver.time.range <- (max(times) - min(times)) / time.points

    # 13. Retrieve group identifiers
    grps <- newprof[, 4]
    # 14. Current label tracker (unused but kept)
    curr.label <- 1

    # 15. Order profile by group then m/z
    newprof <- newprof[order(newprof[, 4], newprof[, 1]), ]
    # 16. Number of ordered rows
    r.newprof <- nrow(newprof)
    # 17. Determine slice boundaries for each group
    breaks <- c(0, which(newprof[1:(r.newprof - 1), 4] != newprof[2:r.newprof, 4]), r.newprof)

    # 18. Helper to estimate density at specified times
    my.dens <- function(x, times, w = NULL, bw, kernel = "gaussian") {
        # 19. Default to equal weights when none provided
        if (is.null(w[1])) w <- rep(1, length(x))
        # 20. Compute weighted kernel density estimate
        d <- density(x, weights = w / sum(w), bw = bw, kernel = kernel)
        # 21. Rescale densities by total weight
        d$y <- d$y * sum(w)
        # 22. Interpolate densities onto target time grid
        d2 <- approx(d$x, d$y, xout = times, rule = 2)
        # 23. Return interpolated density
        d2
    }

    # 24. Create index vector for time grid
    new.times <- 1:length(times)
    # 25. Save original smoother settings to loop over
    all.smoother.window <- smoother.window

    # 26. Iterate through requested smoothing windows
    for (mm in 1:length(all.smoother.window)) {
        # 27. Select current window size
        smoother.window <- all.smoother.window[mm]
        # 28. Initialize matrix to record EIC features
        eic.rec <- matrix(0, nrow = length(breaks) - 1, ncol = 44)
        # 29. Define column names for recorded features
        this.names <- c("grp.label", "mz", "run.length", "raw.density.0", "raw.density.25", "raw.density.50", "raw.density.75", "raw.density.100", "raw.density.100.overall", "b.raw.density.0", "b.raw.density.25", "b.raw.density.50", "b.raw.density.75", "b.raw.density.100", "b.raw.density.100.overall", "signal.height.0", "signal.height.25", "signal.height.50", "signal.height.75", "signal.height.100", "wtd.density.0", "wtd.density.25", "wtd.density.50", "wtd.density.75", "wtd.density.100", "wtd.density.100.overall", "b.wtd.density.0", "b.wtd.density.25", "b.wtd.density.50", "b.wtd.density.75", "b.wtd.density.100", "b.wtd.density.100.overall", "log.w.density.0", "log.w.density.25", "log.w.density.50", "log.w.density.75", "log.w.density.100", "log.w.density.100.overall", "b.log.w.density.0", "b.log.w.density.25", "b.log.w.density.50", "b.log.w.density.75", "b.log.w.density.100", "b.log.w.density.100.overall")
        # 30. Annotate names with bandwidth where applicable
        this.names[-c(1, 2, 3, 16:20)] <- paste("bw", smoother.window, this.names[-c(1, 2, 3, 16:20)], sep = ".")
        # 31. Assign column names to feature matrix
        colnames(eic.rec) <- this.names

        # 32. Choose central time point for single-point density reference
        this.time <- round(median(new.times))
        # 33. Compute density at that point using Gaussian kernel
        single.point.dens <- my.dens(this.time, times = new.times, bw = smoother.window, kernel = "gaussian")
        # 34. Summarize density distribution statistics
        single.point.dens <- c(quantile(single.point.dens$y[this.time], c(0, 0.25, 0.5, 0.75, 1)), max(single.point.dens$y))

        # 35. Repeat computation using rectangular kernel
        single.point.dens.2 <- my.dens(this.time, times = new.times, bw = smoother.window, kernel = "rectangular")
        # 36. Summarize rectangular density statistics
        single.point.dens.2 <- c(quantile(single.point.dens.2$y[this.time], c(0, 0.25, 0.5, 0.75, 1)), max(single.point.dens.2$y))

        # 37. Iterate through each EIC slice
        for (m in 2:length(breaks)) {
            # 38. Extract one group's profile slice
            this.prof <- newprof[(breaks[m - 1] + 1):breaks[m], ]
            # 39. Handle case with single observation
            if (is.null(nrow(this.prof))) {
                # 40. Record features directly for single point
                eic.rec[m - 1, ] <- c(this.prof[4], this.prof[1], 1, single.point.dens, single.point.dens.2, rep(this.prof[3], 5), single.point.dens * this.prof[3], single.point.dens.2 * this.prof[3], single.point.dens * log10(1 + this.prof[3]), single.point.dens.2 * log10(1 + this.prof[3]))
            } else {
                # 41. Order multiple points by time
                this.prof <- this.prof[order(this.prof[, 2]), ]
                # 42. Extract times from profile slice
                this.times <- this.prof[, 2]
                # 43. Extract intensities
                this.intensi <- this.prof[, 3]
                # 44. Extract m/z values
                this.mass <- this.prof[, 1]
                # 45. Placeholder label variable
                this.label <- 1

                # 46. Density of time points without weighting
                dens <- my.dens(x = this.times, times = new.times, bw = smoother.window)
                # 47. Summaries of unweighted density
                dens <- c(quantile(dens$y[this.times], c(0, 0.25, 0.5, 0.75, 1)), max(dens$y))
                # 48. Density using box kernel
                dens.box <- my.dens(x = this.times, times = new.times, bw = smoother.window, kernel = "rectangular")
                # 49. Summaries of box-kernel density
                dens.box <- c(quantile(dens.box$y[this.times], c(0, 0.25, 0.5, 0.75, 1)), max(dens.box$y))

                # 50. Use raw intensity as weights
                y <- this.intensi
                # 51. Density weighted by intensity
                dens2 <- my.dens(x = this.times, times = new.times, bw = smoother.window, w = y)
                # 52. Summaries of intensity-weighted density
                dens2 <- c(quantile(dens2$y[this.times], c(0, 0.25, 0.5, 0.75, 1)), max(dens2$y))
                # 53. Box-kernel version of intensity weighting
                dens2.box <- my.dens(x = this.times, times = new.times, bw = smoother.window, w = y, kernel = "rectangular")
                # 54. Summaries of box-weighted density
                dens2.box <- c(quantile(dens2.box$y[this.times], c(0, 0.25, 0.5, 0.75, 1)), max(dens2.box$y))

                # 55. Apply log transform to intensity
                y <- log10(1 + this.intensi)
                # 56. Density weighted by log intensity
                dens3 <- my.dens(x = this.times, times = new.times, bw = smoother.window, w = y)
                # 57. Summaries of log-weighted density
                dens3 <- c(quantile(dens3$y[this.times], c(0, 0.25, 0.5, 0.75, 1)), max(dens3$y))
                # 58. Box-kernel density with log weights
                dens3.box <- my.dens(x = this.times, times = new.times, bw = smoother.window, w = y, kernel = "rectangular")
                # 59. Summaries of log box-kernel density
                dens3.box <- c(quantile(dens3.box$y[this.times], c(0, 0.25, 0.5, 0.75, 1)), max(dens3.box$y))

                # 60. Record all features for this slice
                eic.rec[m - 1, ] <- c(this.prof[1, 4], median(this.prof[, 1]), diff(range(this.prof[, 2])) + 1, dens, dens.box, quantile(this.prof[, 3], c(0, 0.25, 0.5, 0.75, 1)), dens2, dens2.box, dens3, dens3.box)
            }
        }
        # 61. Initialize or append to full record
        if (mm == 1) {
            # 62. First smoothing window simply stores record
            all.eic.rec <- eic.rec
        } else {
            # 63. Subsequent windows append columns excluding duplicates
            all.eic.rec <- cbind(all.eic.rec, eic.rec[, -c(1, 2, 3, 16:20)])
        }
    }

    # 64. Convert indexed times back to actual retention times
    newprof[, 2] <- times[newprof[, 2]]
    # 65. Append height features to EIC record
    all.eic.rec <- cbind(all.eic.rec, height.rec[, c(-1, -3)])
    # 66. Return processed profile and extracted features
    return(list(prof = newprof, eic.ftrs = all.eic.rec))
}
