#' Plot extracted ion chromatograms
#'
#' Given an output object from the function cdf.to.ftr(), this function plots the EICs selected by the user.
#'
#' @param aligned An output object from cdf.to.ftr().
#' @param rows A numeric vector selecting which rows of the aligned feature table to be plotted.
#' @param colors The colors (one per profile) the user wishes to use for the plots. The default is NA, in which case a default color set is used.
#' @param transform There are four possible values. "none": the original intensity data is plotted; "log": the intensity data is transformed by log(x+1); "sqrt": the intensity data is square root transformed; "cuberoot": the intensity data is cube root transformed.
#' @param subset The user can choose a subset of the profiles for which the EICs are plotted. It is given as a vector of profile indecies. The default is NA, in which case the EICs from all the profiles are plotted.
#' @param min.run The min.run parameter used in the proc.cdf() step.
#' @param min.pres The min.pres parameter used in the proc.cdf() step.
#' @param max.spline.time.points The maximum time points to use in spline fit.
#'
#' @details The EICs are plotted as overlaid line plots. The graphic device is divided into four parts, each of which is used to plot one EIC. When all four parts are occupied, the function calls x11() to open another graphic device. The colors used (one per profile) is printed in the command window.
#'
#' @return There is no return value.
#'
#' @references Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
EIC.plot <- function(aligned, rows = NA, colors = NA, transform = "none", subset = NA, min.run, min.pres, max.spline.time.points = 1000) {
    # 1. Validate transformation option
    if (!(transform %in% c("none", "log", "sqrt", "cuberoot"))) {
        # 2. Warn user for invalid choice
        message("Invalid transformation. It has to be from: none, log, sqrt, and cuberoot")
        break;
    }
    # 3. Proceed only if rows specified
    if (!is.na(rows[1])) {
        # 4. Use spline functions for time adjustment
        library(splines)
        # 5. Number of experiments/profiles
        num.exp <- nrow(summary(aligned$features))
        # 6. Default to all profiles if subset not given
        if (is.na(subset[1])) subset <- 1:num.exp
        # 7. Assign default colors if none provided
        if (is.na(colors[1])) colors <- c("red", "blue", "dark blue", "orange", "green", "yellow", "cyan", "pink", "violet", "bisque", "azure", "brown", "chocolate", rep("grey", length(subset)))
        # 8. Extract file names for profiles
        files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
        # 9. Construct raw profile file names based on parameters
        rawprof.names <- paste(unlist(strsplit(files, "\\."))[seq(1, 2 * num.exp, by = 2)], "_", min.run, "_", min.pres, "_", aligned$mz.tol, ".rawprof", sep = "")

        # 10. Store adjusted time vectors for each profile
        adj.times <- new("list")
        for (i in subset) {
            # 11. Load raw profile for current experiment
            load(rawprof.names[i])
            # 12. Retrieve and sort unique time points
            times <- unique(raw.prof$labels)
            times <- times[order(times)]

            # 13. Extract original and adjusted times where intensities exist
            orig.time <- aligned$features[[i]][!is.na(aligned$features[[i]][, 3]), 2]
            adjusted.time <- aligned$features2[[i]][!is.na(aligned$features2[[i]][, 3]), 2]
            # 14. Round times to avoid floating discrepancies
            orig.time <- round(orig.time, 8)
            adjusted.time <- round(adjusted.time, 8)
            # 15. Frequency of adjusted and original times
            ttt <- table(adjusted.time)
            ttt2 <- table(orig.time)
            # 16. Keep time points unique in both series
            to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt == 1]) & orig.time %in% as.numeric(names(ttt2)[ttt2 == 1]))
            # 17. Limit number of points for spline fitting
            if (length(to.use) > max.spline.time.points) to.use <- sample(to.use, max.spline.time.points, replace = FALSE)
            orig.time <- orig.time[to.use]
            adjusted.time <- adjusted.time[to.use]

            # 18. Fit spline between original and adjusted times
            sp <- interpSpline(adjusted.time ~ orig.time)
            # 19. Predict adjusted times for all reference times
            adj.times[[i]] <- predict(sp, times)$y
        }
        # 20. Iterate over requested rows to plot
        for (n in 1:length(rows)) {
            # 21. Reset plotting panel every four plots
            if (n %% 4 == 1) {
                par(mfrow = c(2, 2))
            }
            # 22. Row index currently processed
            this.row <- rows[n]
            # 23. Intensities across profiles for this feature
            this.intensi <- aligned$final.ftrs[this.row, 5:(num.exp + 4)]
            # 24. Corresponding times per profile
            this.times <- aligned$final.times[this.row, 5:(num.exp + 4)]
            # 25. Feature m/z value
            this.mz <- aligned$final.ftrs[this.row, 1]
            # 26. Container for plotted data
            to.plot <- new("list")
            # 27. Track max intensity and time range for scaling
            y.max <- 0
            x.max <- 0
            x.min <- Inf

            # 28. Loop through each selected profile
            for (iii in 1:length(subset)) {
                i <- subset[iii]
                # 29. Only proceed if feature detected in profile
                if (this.intensi[i] != 0) {
                    # 30. Load raw profile
                    load(rawprof.names[i])
                    times <- unique(raw.prof$labels)
                    times <- times[order(times)]

                    # 31. Find matching feature within aligned data
                    mz.diff <- abs(aligned$features2[[i]][, 1] - this.mz)
                    time.diff <- abs(aligned$features2[[i]][, 2] - this.times[i])
                    sel <- which(time.diff < 1e-17)
                    if (length(sel) > 1) sel <- sel[which(mz.diff[sel] == min(mz.diff[sel]))[1]]

                    # 32. Handle case of merged peaks lacking exact match
                    if (length(sel) == 0) {
                        mz.lim <- aligned$final.ftrs[this.row, c(3, 4)]
                        sel <- which(aligned$features2[[i]][, 1] >= mz.lim[1] & aligned$features2[[i]][, 1] <= mz.lim[2] & !is.na(aligned$features2[[i]][, 6]))
                        sub.features <- aligned$features2[[i]][sel, ]
                        sub.time.diff <- abs(sub.features[, 2] - this.times[i])
                        sel <- sel[which(sub.time.diff == min(sub.time.diff))][1]
                    }

                    # 33. Retrieve target m/z and original time
                    target.mz <- aligned$features2[[i]][sel, 1]
                    target.time <- aligned$features[[i]][sel, 2]
                    time.adjust <- aligned$features2[[i]][sel, 2] - aligned$features[[i]][sel, 2]

                    # 34. Find slice in raw profile closest to target m/z
                    mz.diff <- abs(raw.prof$masses - target.mz)
                    sel.slice <- raw.prof$grps[which(mz.diff == min(mz.diff))[1]]
                    sel.time.range <- range(raw.prof$labels[raw.prof$grps == sel.slice])

                    # 35. Ensure target time falls within slice; otherwise search neighboring slices
                    while (target.time < sel.time.range[1] | target.time > sel.time.range[2]) {
                        mz.diff[raw.prof$grps == sel.slice] <- 100
                        sel.slice <- raw.prof$grps[which(mz.diff == min(mz.diff))[1]]
                        sel.time.range <- range(raw.prof$labels[raw.prof$grps == sel.slice])
                    }

                    # 36. Extract time and intensity series for selected slice
                    sel.time <- raw.prof$labels[raw.prof$grps == sel.slice]
                    sel.intensi <- raw.prof$intensi[raw.prof$grps == sel.slice]
                    sel.intensi <- sel.intensi[order(sel.time)]
                    sel.time <- sel.time[order(sel.time)]

                    # 37. Align intensities to full time grid
                    all.intensi <- times * 0
                    all.intensi[times %in% sel.time] <- sel.intensi
                    all.times <- adj.times[[i]]

                    # 38. Apply requested intensity transformation
                    if (transform == "log") all.intensi <- log10(all.intensi + 1)
                    if (transform == "sqrt") all.intensi <- sqrt(all.intensi)
                    if (transform == "cuberoot") all.intensi <- all.intensi^(1 / 3)

                    # 39. Update plot scaling limits
                    if (max(all.intensi) > y.max) y.max <- max(all.intensi)
                    if (max(all.times[all.intensi > 0]) > x.max) x.max <- max(all.times[all.intensi > 0])
                    if (min(all.times[all.intensi > 0]) < x.min) x.min <- min(all.times[all.intensi > 0])
                    # 40. Store series for later plotting
                    to.plot[[iii]] <- cbind(all.times, all.intensi)
                }
            }

            # 41. Plot compiled EICs if any were found
            if (!is.null(nrow(summary(to.plot)))) {
                for (iii in 1:min(length(subset), nrow(summary(to.plot)))) {
                    i <- subset[iii]
                    if (iii == 1) {
                        # 42. Draw first profile establishing axes
                        plot(to.plot[[iii]], xlim = c(x.min, x.max), ylim = c(0, y.max), type = "l", col = colors[iii], xlab = "retention time", ylab = "intensity", main = paste("EIC for row", rows[n], ",m/z", round(this.mz, 4)))
                    } else {
                        # 43. Add subsequent profiles
                        lines(to.plot[[iii]], col = colors[iii])
                    }
                }
            }
        }
        # 44. Print color assignments to console
        message("the colors used are:")
        print(paste(subset, ": ", colors[1:length(subset)], sep = ""))
    }
}
