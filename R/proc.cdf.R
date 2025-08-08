#' Filter noise and detect peaks from LC/MS data in CDF format
#'
#' This function applies the run filter to remove noise. Data points are grouped into EICs in this step.
#'
#' @param filename The CDF file name. If the file is not in the working directory, the path needs to be given.
#' @param min.pres Run filter parameter. The minimum proportion of presence in the time period for a series of signals grouped by m/z to be considered a peak.
#' @param min.run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.
#' @param tol m/z tolerance level for the grouping of data points. This value is expressed as the fraction of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level. The recommended value is the machine's nominal accuracy level. Divide the ppm value by 1e6. For FTMS, 1e-5 is recommended.
#' @param baseline.correct After grouping the observations, the highest intensity in each group is found. If the highest is lower than this value, the entire group will be deleted. The default value is NA, in which case the program uses the 75th percentile of the height of the noise groups.
#' @param baseline.correct.noise.percentile The perenctile of signal strength of those EIC that don't pass the run filter, to be used as the baseline threshold of signal strength.
#' @param do.plot Whether to produce diagnostic plots.
#' @param intensity.weighted Whether to weight the local density by signal intensities.
#'
#' @details The m/z are grouped based on the tolerance level using multi-stage smoothing and peak finding. Non-parametric density estimation is used in both m/z dimension and elution time dimension to fine-tune the signal grouping. A run filter is applied, which requires a "true peak" to have a minimum length in the retention time dimension (parameter: min.run), as well as being detected at or higher than a proportion of the time points within the time period (parameter: min.pres).
#'
#' @return A matrix with four columns: m/z value, retention time, intensity, and group number.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
# 1. Main function for pre-processing a single CDF profile: binning into EICs, run filter, baseline thresholding, and QC plots
proc.cdf <- function(filename, min.pres = 0.5, min.run = 12, tol = 1e-5, baseline.correct = 0, baseline.correct.noise.percentile = 0, do.plot = TRUE, intensity.weighted = FALSE) {
    # 2. Set up plotting window if requested to show progression and QC histograms
    if (do.plot) par(mfrow = c(2, 2))
    # 3. Construct cache file name based on parameters for reproducible reruns
    this.name <- paste(strsplit(tolower(filename), "\\.")[[1]][1], "_", min.run, "_", min.pres, "_", tol, ".rawprof", sep = "")
    all.files <- dir()
    is.there <- all.files[which(all.files == this.name)]

    # 4. Load cached profile if available to skip recomputation (also displays tol)
    if (length(is.there) > 0) {
        load(this.name)
        if (do.plot) plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", main = "tolerance level loaded", axes = FALSE)
        if (do.plot) text(x = 0, y = 0, tol, cex = 1.2)
    } else {
        # 5. Read spectrum, drop NA rows, warn, then create EICs via adaptive.bin()
        this <- load.lcms(filename)
        na.sel <- c(which(is.na(this$masses)), which(is.na(this$labels)), which(is.na(this$intensi)))
        if (length(na.sel) > 0) {
            na.sel <- unique(na.sel)
            this$masses <- this$masses[-na.sel]
            this$labels <- this$labels[-na.sel]
            this$intensi <- this$intensi[-na.sel]
        }
        message("WARNING: there are NA values in the m/z or intensity. Check the file:")
        message(filename)
        # 6. Build grouped raw profile; optionally weight by intensity during mass density estimation
        raw.prof <- adaptive.bin(this, min.run = min.run, min.pres = min.pres, tol = tol, baseline.correct = baseline.correct, weighted = intensity.weighted)
        save(raw.prof, file = this.name)
    }

    # 7. Combine EIC profile data into matrix form (mz, time, intensity, group)
    newprof <- cbind(raw.prof$masses, raw.prof$labels, raw.prof$intensi, raw.prof$grps)
    # 8. Summarize group maximum heights for noise vs. signal diagnostics and baseline selection
    h.1 <- log10(raw.prof$height.rec[raw.prof$height.rec[, 2] <= max(2, raw.prof$min.count.run * min.pres / 2), 3])
    h.2 <- log10(raw.prof$height.rec[raw.prof$height.rec[, 2] >= raw.prof$min.count.run * min.pres, 3])
    if (do.plot) {
        # 9. Display distribution of (log) group heights, or fall back to rug-like lines when too few noise groups
        if (length(h.1) > 50) {
            plot(density(h.1), xlab = "maximum height of group (log scale)", xlim = range(c(h.1, h.2)), main = "Black - noise groups \n Blue - selected groups")
        } else {
            plot(NA, NA, xlab = "maximum height of group (log scale)", xlim = range(c(h.1, h.2)), ylim = c(0, 1), main = "Black - noise groups \n Blue - selected groups")
            if (length(h.1) > 0) abline(v = h.1)
        }
    }

    # 10. Determine baseline cutoff: either auto by noise percentile or user-provided; draw cutoff
    if (is.na(baseline.correct)) {
        baseline.correct <- 10^quantile(h.1, baseline.correct.noise.percentile)
        message(c("maximal height cut is automatically set at the", baseline.correct.noise.percentile, "percentile of noise group heights: ", baseline.correct))
    } else {
        message(c("maximal height cut is provided by user: ", baseline.correct))
    }
    if (do.plot) abline(v = log10(baseline.correct), col = "red")

    # 11. Keep groups that meet minimal presence and exceed baseline; then apply continuity (run) filter for final EICs
    if (is.na(baseline.correct)) baseline.correct <- 0
    run.sel <- raw.prof$height.rec[which(raw.prof$height.rec[, 2] >= raw.prof$min.count.run * min.pres & raw.prof$height.rec[, 3] > baseline.correct), 1]
    newprof <- newprof[newprof[, 4] %in% run.sel, ]
    # 12. Continuity index refines EICs; store QC lines/hists when plotting is on
    new.prof <- cont.index(newprof, min.pres = min.pres, min.run = min.run)
    if (do.plot) lines(density(log10(new.prof$height.rec)), col = "blue")
    time.range.rec <- new.prof$time.range.rec
    mz.pres.rec <- new.prof$mz.pres.rec
    if (do.plot) hist(time.range.rec, xlab = "Range of retention time in the same group", ylab = "Density", freq = FALSE, nclass = 100, main = "Group retention time range distribution")
    if (do.plot) hist(mz.pres.rec, xlab = "% signal present in the same group", ylab = "Density", freq = FALSE, nclass = 20, main = "Group % present signal distribution")

    # 13. Return the final matrix of filtered EICs (mz, time, intensity, group)
    return(new.prof$new.rec)
}
