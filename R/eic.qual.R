#' Internal function: Calculate the single predictor quality.
#'
#' For each column of an EIC data feature matrix, find its predictive power on the m/z match to known metabolites.
#'
#' @param eic.rec The EIC data feature matrix. Each row is an EIC. Each column is a data feature.
#' @param known.mz The known m/z values to be matched to.
#' @param mass.matched A vector of indicators of whether the m/z of each EIC is matched to the known m/z values. The default is NA, in which case it is calculated within the function.
#' @param match.tol.ppm The tolerance level of m/z match.
#' @param do.plot Whether to produce plots of the ROCS.
#' @param pos.confidence The confidence level for the features matched to the known feature list.
#' @param neg.confidence The confidence level for the features not matching to the known feature list.
#'
#' @return A matrix of four columns. The first two columns are the VUS and AUC without uncertainty. The next two columns are the VUS and AUC with uncertainty.
#'
#' @references Bioinformatics. 30(20): 2941-2948.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
eic.qual <- function(eic.rec, known.mz, mass.matched = NA, match.tol.ppm = 5, do.plot = FALSE, pos.confidence = 0.99, neg.confidence = 0.99) {
    # 1. Determine if match indicators need to be computed
    if (is.na(mass.matched[1])) {
        # 2. Calculate m/z matches against known list within tolerance
        y <- mass.match(x = eic.rec[, 2], known.mz = known.mz, match.tol.ppm = match.tol.ppm)
    } else {
        # 3. Use provided match indicators
        y <- mass.matched
    }

    # 4. Prepare containers for ROC results with and without uncertainty
    rec.def <- rec.unc <- new("list")
    # 5. Matrix to store volume-under-surface and AUC metrics
    vus.fcauc <- matrix(nrow = ncol(eic.rec) - 2, ncol = 4)
    # 6. Label columns for deterministic and uncertain estimates
    colnames(vus.fcauc) <- c("vus_1", "fcauc_1", "vus_uncertainty", "fcauc_uncertainty")
    # 7. Assign feature names to rows
    rownames(vus.fcauc) <- colnames(eic.rec)[-1:-2]

    # 8. Evaluate each feature column beyond m/z and retention time
    for (i in 3:ncol(eic.rec)) {
        # 9. Extract current feature values
        x <- eic.rec[, i]
        # 10. Compute ROC under perfect confidence
        r <- rocs.x(x[y == 0], x[y == 1], rep(1, sum(y == 0)), rep(1, sum(y == 1)), n.perm = 1, FDR.cut = 1, do.plot = FALSE)
        # 11. Store deterministic ROC results
        rec.def[[i - 2]] <- r
        # 12. Save VUS and AUC without uncertainty
        vus.fcauc[i - 2, 1:2] <- unlist(r[1:2])

        # 13. Recompute ROC with uncertainty weights
        r <- rocs.x(x[y == 0], x[y == 1], rep(0.99, sum(y == 0)), rep(0.99, sum(y == 1)), n.perm = 1, FDR.cut = 1, do.plot = FALSE)
        # 14. Store uncertainty-adjusted ROC
        rec.unc[[i - 2]] <- r
        # 15. Save VUS and AUC with uncertainty
        vus.fcauc[i - 2, 3:4] <- unlist(r[1:2])
    }

    # 16. Optionally visualize ROC curves for up to nine features
    if (do.plot) {
        # 17. Arrange plotting area
        par(mfrow = c(3, 3))
        # 18. Iterate over a subset of results for display
        for (i in 1:min(9, length(rec.unc))) {
            # 19. Plot ROC curve without uncertainty
            plot(rec.def[[i]]$fp, rec.def[[i]]$tp, type = "l", col = "blue", xlab = "FPR", ylab = "TPR", main = rownames(vus.fcauc)[i])
            # 20. Add ROC curve with uncertainty
            lines(rec.unc[[i]]$fp, rec.unc[[i]]$tp, col = "red")
            # 21. Draw diagonal reference line
            abline(0, 1, col = "grey", lty = 2)
        }
    }
    # 22. Return summary matrix of VUS and AUC metrics
    vus.fcauc
}
