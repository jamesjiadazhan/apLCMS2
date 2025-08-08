#' Peak detection using the machine learning approach.
#'
#' The procedure uses information of known metabolites, and constructs prediction models to differentiate EICs.
#'
#' @param filename The cdf file name. If the file is not in the working directory, the path needs to be given.
#' @param tol m/z tolerance level for the grouping of data points. This value is expressed as the fraction of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level. The recommended value is the machine's nominal accuracy level. Divide the ppm value by 1e6. For FTMS, 1e-5 is recommended.
#' @param min.run Run filter parameter. The minimum length of elution time for a series of signals grouped by m/z to be considered a peak.
#' @param min.pres Run filter parameter. The minimum proportion of presence in the time period for a series of signals grouped by m/z to be considered a peak.
#' @param baseline.correct After grouping the observations, the highest intensity in each group is found. If the highest is lower than this value, the entire group will be deleted. The default value is NA, in which case the program uses the 75th percentile of the height of the noise groups.
#' @param ridge.smoother.window The size of the smoother window used by the kernel smoother to remove long ridge noise from each EIC.
#' @param smoother.window The smoother windows to use in data feature generation.
#' @param known.mz The m/z values of the known metabolites.
#' @param match.tol.ppm The ppm tolerance to match identified features to known metabolites/features.
#' @param do.plot Whether to produce diagnostic plots.
#' @param pos.confidence The confidence level for the features matched to the known feature list.
#' @param neg.confidence The confidence level for the features not matching to the known feature list.
#' @param max.ftrs.to.use The maximum number of data features to use in a predictive model.
#' @param do.grp.reduce Whether to reduce data features that are similar. It is based on data feature predictability.
#' @param remove.bottom.ftrs The number of worst performing data features to remove before model building.
#' @param max.fpr The proportion of unmatched features to be selected in the feature detection step.
#' @param min.tpr The proportion of matched features to be selected in the feature detection step.
#' @param intensity.weighted Whether to weight the local density by signal intensities.
#'
#' @details The subroutine takes CDF, mxml etc LC/MS profile.
#' First the profile is sliced into EICs using adaptive binning. Then data features are extracted from each EIC. The EICs are classified into two groups: those that have m/z values that match to known m/z values, and those that don't. Classification models are built to separate the two classes, and each EIC is given a score by the classification model. Those with better scores are selected to enter the feature quantification step.
#'
#' @return A matrix with four columns: m/z value, retention time, intensity, and group number.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
learn.cdf <- function(filename, tol = 2e-5, min.run = 4, min.pres = 0.3, baseline.correct = 0, ridge.smoother.window = 50, smoother.window = c(1, 5, 10), known.mz, match.tol.ppm = 5, do.plot = FALSE, pos.confidence = 0.99, neg.confidence = 0.99, max.ftrs.to.use = 10, do.grp.reduce = TRUE, remove.bottom.ftrs = 0, max.fpr = seq(0, 0.6, by = 0.1), min.tpr = seq(0.8, 1, by = 0.1), intensity.weighted = FALSE) {
    # 1. Optionally set up plotting layout
    if (do.plot) par(mfrow = c(2, 2))

    # 2. Construct cache file name based on parameters
    this.name <- paste(strsplit(tolower(filename), "\\.")[[1]][1], "_", tol, "_", ridge.smoother.window, "_", baseline.correct, ".rawlearn", sep = "")
    # 3. List files in working directory
    all.files <- dir()
    # 4. Check if cached profile exists
    is.there <- all.files[which(all.files == this.name)]

    # 5. Load cached profile if available
    if (length(is.there) > 0) {
        load(this.name)
        if (do.plot) {
            plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", main = "tolerance level loaded", axes = FALSE)
            text(x = 0, y = 0, tol, cex = 1.2)
        }
    } else {
        # 6. Load raw LC/MS profile
        this <- load.lcms(filename)
        # 7. Group signals into EICs using adaptive binning
        raw.prof <- adaptive.bin.2(this, tol = tol, ridge.smoother.window = ridge.smoother.window, baseline.correct = baseline.correct, weighted = intensity.weighted)
        # 8. Save processed profile for reuse
        save(raw.prof, file = this.name)
    }

    # 9. Ensure baseline threshold is numeric
    if (is.na(baseline.correct)) baseline.correct <- 0
    # 10. Identify groups meeting run filter and baseline criteria
    run.sel <- raw.prof$height.rec[which(raw.prof$height.rec[, 2] >= min.run * min.pres & raw.prof$height.rec[, 3] > baseline.correct), 1]
    # 11. Assemble profile matrix and keep selected groups
    newprof <- cbind(raw.prof$masses, raw.prof$labels, raw.prof$intensi, raw.prof$grps)
    newprof <- newprof[newprof[, 4] %in% run.sel, ]
    raw.prof$height.rec <- raw.prof$height.rec[raw.prof$height.rec[, 1] %in% newprof[, 4], ]

    # 12. Apply continuity index to further refine peaks
    new.prof <- cont.index(newprof, min.pres = min.pres, min.run = min.run)
    raw.prof$height.rec <- raw.prof$height.rec[raw.prof$height.rec[, 1] %in% new.prof$new.rec[, 4], ]
    raw.prof$masses <- new.prof$new.rec[, 1]
    raw.prof$labels <- new.prof$new.rec[, 2]
    raw.prof$intensi <- new.prof$new.rec[, 3]
    raw.prof$grps <- new.prof$new.rec[, 4]
    raw.prof$height.rec[, 3] <- new.prof$height.rec

    # 13. Generate EIC features and profiles
    eic.rec.0 <- eic.disect(raw.prof, smoother.window = smoother.window)
    prof <- eic.rec.0$prof
    eic.rec.0 <- eic.rec.0$eic.ftrs

    # 14. Predict feature quality using machine learning
    a <- eic.pred(eic.rec.0, known.mz, to.use = max.ftrs.to.use, match.tol.ppm = match.tol.ppm, do.grp.reduce = do.grp.reduce, remove.bottom = remove.bottom.ftrs, max.fpr = max.fpr[1], min.tpr = min.tpr[1], do.plot = do.plot)

    # 15. Enumerate combinations of FPR and TPR thresholds
    fpr.tpr.combo <- expand.grid(max.fpr, min.tpr)
    # 16. Prepare list to store profiles for each threshold pair
    prof.sel <- new("list")

    # 17. Evaluate selection for every FPR/TPR combination
    for (i in 1:nrow(fpr.tpr.combo)) {
        # 18. Initialize selection mask
        chosen <- a$matched * 0
        # 19. Retain matched features meeting TPR cutoff
        chosen[a$matched == 1 & a$tpr <= fpr.tpr.combo[i, 2]] <- 1
        # 20. Retain unmatched features meeting FPR cutoff
        chosen[a$matched == 0 & a$fpr <= fpr.tpr.combo[i, 1]] <- 1

        # 21. Locate selected features
        selected <- which(chosen == 1)
        # 22. Store corresponding profile subset
        prof.sel[[i]] <- prof[which(prof[, 4] %in% raw.prof$height.rec[selected, 1]), ]
        # 23. Name list element using thresholds
        names(prof.sel)[[i]] <- paste("fpr", fpr.tpr.combo[i, 1], "tpr", fpr.tpr.combo[i, 2], sep = "_")
    }
    # 24. Save combo grid and model details in result
    prof.sel$fpr.tpr.combo <- fpr.tpr.combo
    prof.sel$model.detail <- a

    # 25. Return list of selected profiles for each threshold setting
    return(prof.sel)
}
