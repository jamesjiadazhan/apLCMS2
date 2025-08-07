#' Plot extracted ion chromatograms based on the machine learning method output
#'
#' Given an output object from the function semi.sup.learn(), this function plots the EICs selected by the user.
#'
#' @param aligned An output object from cdf.to.ftr().
#' @param rows A numeric vector selecting which rows of the aligned feature table to be plotted.
#' @param colors The colors (one per profile) the user wishes to use for the plots. The default is NA, in which case a default color set is used.
#' @param transform There are four possible values. "none": the original intensity data is plotted; "log": the intensity data is transformed by log(x+1); "sqrt": the intensity data is square root transformed; "cuberoot": the intensity data is cube root transformed.
#' @param subset The user can choose a subset of the profiles for which the EICs are plotted. It is given as a vector of profile indecies. The default is NA, in which case the EICs from all the profiles are plotted.
#' @param tol The mz tolerance level used in learn.cdf().
#' @param ridge.smoother.window The ridge.smoother.window parameter value used in learn.cdf().
#' @param baseline.correct The baseline.correct parameter value used in learn.cdf().
#' @param max.spline.time.points The maximum number of points to use in the spline fit along the retention time axis.
#'
#' @details The function plots a single EIC. It plots intensity against retention time. It uses different color for different profiles.
#'
#' @return There is no return value.
#'
#' @references Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
EIC.plot.learn <- function(aligned, rows = NA, colors = NA, transform = "none", subset = NA, tol = 2.5e-5, ridge.smoother.window = 50, baseline.correct = 0, max.spline.time.points = 1000) {
        # 1. Validate that requested transformation is allowed
        if (!(transform %in% c("none", "log", "sqrt", "cuberoot"))) {
            # 2. Inform user when transformation is invalid
            message("Invalid transformation. It has to be from: none, log, sqrt, and cuberoot")
            # 3. Exit the function early on invalid input
            break;
        }
        # 4. Proceed only when rows are specified
        if (!is.na(rows[1])) {
            # 5. Load spline functions for time adjustment
            library(splines)
            # 6. Count number of experiments in the aligned object
            num.exp <- nrow(summary(aligned$features))
            # 7. Default to using all profiles if subset not provided
            if (is.na(subset[1])) subset <- 1:num.exp
            # 8. Provide default color palette when colors not specified
            if (is.na(colors[1])) colors <- c("red", "blue", "dark blue", "orange", "green", "yellow", "cyan", "pink", "violet", "bisque", "azure", "brown", "chocolate", rep("grey", length(subset)))
            # 9. Retrieve raw data filenames from aligned features
            files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
            # 10. Compose expected raw profile filenames
            rawprof.names <- paste(unlist(strsplit(tolower(files), "\\."))[seq(1, 2 * num.exp, by = 2)], "_", tol, "_", ridge.smoother.window, "_", baseline.correct, ".rawlearn", sep = "")

            # 11. Prepare container for adjusted time axes
            adj.times <- new("list")
            # 12. Loop through selected profiles
            for (i in subset) {
                # 13. Load raw profile for this experiment
                load(rawprof.names[i])
                # 14. Extract unique retention times and sort them
                times <- unique(raw.prof$labels)
                times <- times[order(times)]

                # 15. Obtain original and adjusted times from features
                orig.time <- aligned$features[[i]][!is.na(aligned$features[[i]][, 3]), 2]
                adjusted.time <- aligned$features2[[i]][!is.na(aligned$features2[[i]][, 3]), 2]
                # 16. Round times for spline calculation stability
                orig.time <- round(orig.time, 8)
                adjusted.time <- round(adjusted.time, 8)
                # 17. Tabulate occurrences of adjusted and original times
                ttt <- table(adjusted.time)
                ttt2 <- table(orig.time)
                # 18. Keep unique time pairs shared by both tables
                to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt == 1]) & orig.time %in% as.numeric(names(ttt2)[ttt2 == 1]))
                # 19. Sample time points when too many are available
                if (length(to.use) > max.spline.time.points) to.use <- sample(to.use, max.spline.time.points, replace = FALSE)
                # 20. Restrict times to selected subset for spline fitting
                orig.time <- orig.time[to.use]
                adjusted.time <- adjusted.time[to.use]

                # 21. Fit spline mapping from original to adjusted time
                sp <- interpSpline(adjusted.time ~ orig.time)
                # 22. Predict adjusted time axis for all raw times
                adj.times[[i]] <- predict(sp, times)$y
            }
            # 23. Iterate through requested feature rows
            for (n in 1:length(rows)) {
                # 24. Determine row index being plotted
                this.row <- rows[n]
                # 25. Extract intensities across profiles for this feature
                this.intensi <- aligned$final.ftrs[this.row, 5:(num.exp + 4)]
                # 26. Extract retention times for same feature
                this.times <- aligned$final.times[this.row, 5:(num.exp + 4)]
                # 27. Record m/z value for labeling
                this.mz <- aligned$final.ftrs[this.row, 1]
                # 28. Prepare list to store processed EIC traces
                to.plot <- new("list")
                # 29. Track maximal intensity for plot scaling
                y.max <- 0
                # 30. Track maximum time for x-axis
                x.max <- 0
                # 31. Track minimum time for x-axis
                x.min <- Inf

                # 32. Loop through subset of profiles for plotting
                for (iii in 1:length(subset)) {
                    # 33. Map subset index to actual profile index
                    i <- subset[iii]
                    # 34. Skip profiles with zero intensity
                    if (this.intensi[i] != 0) {
                        # 35. Load raw profile data for this experiment
                        load(rawprof.names[i])
                        # 36. Obtain and sort unique time points
                        times <- unique(raw.prof$labels)
                        times <- times[order(times)]

                        # 37. Compute m/z and time differences to target feature
                        mz.diff <- abs(aligned$features2[[i]][, 1] - this.mz)
                        time.diff <- abs(aligned$features2[[i]][, 2] - this.times[i])
                        # 38. Locate matching peak by retention time
                        sel <- which(time.diff < 1e-17)
                        # 39. When multiple matches, choose closest m/z
                        if (length(sel) > 1) sel <- sel[which(mz.diff[sel] == min(mz.diff[sel]))[1]]
                        # 40. Note alignment merges via shared labels
                        # this means two peaks are merged at the alignment step
                        # they share labels in the aligned$features2 object
                        # 41. Handle missing match by searching within m/z limits
                        if (length(sel) == 0) {
                                # 42. Retrieve m/z range for current feature
                                mz.lim <- aligned$final.ftrs[this.row, c(3, 4)]
                                # 43. Find candidate peaks within m/z window and valid labels
                                sel <- which(aligned$features2[[i]][, 1] >= mz.lim[1] & aligned$features2[[i]][, 1] <= mz.lim[2] & !is.na(aligned$features2[[i]][, 6]))
                                # 44. Subset to candidate feature rows
                                sub.features <- aligned$features2[[i]][sel, ]
                                # 45. Choose candidate closest in time
                                sub.time.diff <- abs(sub.features[, 2] - this.times[i])
                                sel <- sel[which(sub.time.diff == min(sub.time.diff))][1]
                            }

                        # 46. Record final m/z and time from selected peak
                        target.mz <- aligned$features2[[i]][sel, 1]
                        target.time <- aligned$features[[i]][sel, 2]
                        # 47. Compute time adjustment between aligned and raw
                        time.adjust <- aligned$features2[[i]][sel, 2] - aligned$features[[i]][sel, 2]

                        # 48. Locate slice in raw profile matching target m/z
                        mz.diff <- abs(raw.prof$masses - target.mz)
                        sel.slice <- raw.prof$grps[which(mz.diff == min(mz.diff))[1]]
                        sel.time.range <- range(raw.prof$labels[raw.prof$grps == sel.slice])

                        # 49. Ensure slice covers the target time by adjusting selection
                        while (target.time < sel.time.range[1] | target.time > sel.time.range[2]) {
                            # 50. Mark current slice as invalid and search next closest
                            mz.diff[raw.prof$grps == sel.slice] <- 100
                            sel.slice <- raw.prof$grps[which(mz.diff == min(mz.diff))[1]]
                            sel.time.range <- range(raw.prof$labels[raw.prof$grps == sel.slice])
                        }

                        # 51. Extract times and intensities for selected slice
                        sel.time <- raw.prof$labels[raw.prof$grps == sel.slice]
                        sel.intensi <- raw.prof$intensi[raw.prof$grps == sel.slice]
                        # 52. Order points by time for plotting
                        sel.intensi <- sel.intensi[order(sel.time)]
                        sel.time <- sel.time[order(sel.time)]

                        # 53. Initialize intensity vector over full time axis
                        all.intensi <- times * 0
                        # 54. Insert observed intensities into full grid
                        all.intensi[times %in% sel.time] <- sel.intensi
                        # 55. Retrieve adjusted times for this profile
                        all.times <- adj.times[[i]]

                        # 56. Apply transformation when requested
                        if (transform == "log") all.intensi <- log10(all.intensi + 1)
                        if (transform == "sqrt") all.intensi <- sqrt(all.intensi)
                        if (transform == "cuberoot") all.intensi <- all.intensi^(1 / 3)

                        # 57. Track global plot limits for all profiles
                        if (max(all.intensi) > y.max) y.max <- max(all.intensi)
                        if (max(all.times[all.intensi > 0]) > x.max) x.max <- max(all.times[all.intensi > 0])
                        if (min(all.times[all.intensi > 0]) < x.min) x.min <- min(all.times[all.intensi > 0])
                        # 58. Store time-intensity pairs for this profile
                        to.plot[[i]] <- cbind(all.times, all.intensi)
                    }
                }

                # 59. Plot collected traces in requested colors
                for (iii in 1:min(length(subset), nrow(summary(to.plot)))) {
                    i <- subset[iii]
                    if (iii == 1) {
                        # 60. Create base plot for first profile
                        plot(to.plot[[i]], xlim = c(x.min, x.max), ylim = c(0, y.max), type = "l", col = colors[iii], xlab = "retention time", ylab = "intensity", main = paste("EIC for row", rows[n], ",m/z", round(this.mz, 4)))
                    } else {
                        # 61. Add additional profiles to the plot
                        lines(to.plot[[i]], col = colors[iii])
                    }
                }
            }
            # 62. Display color assignments used in the plot
            message("the colors used are:")
            # 63. Print mapping of profile indices to colors
            print(paste(subset, ": ", colors[1:length(subset)], sep = ""))
        }
    }
