#' Internal function: Updates the information of a feature for the known feature table.
#'
#' The function takes the information about the feature in the known feature table (if available), and updates it using the information found in the current dataset.
#'
#' @param existing.row The existing row in the known feature table (detailed in the semi.sup() document).
#' @param ftrs.row The row of the matched feature in the new aligned feature table.
#' @param chr.row The row of the matched feature in the new retention time table of aligned features.
#'
#' @details The function calculates and updates the mean intensity, variation of intensity, mean retention time etc.
#'
#' @return A vector, the updated row for the known feature table.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
peak.characterize <- function(existing.row = NA, ftrs.row, chr.row) {
    # 1. Helper to merge new observations with existing statistics
    merge.new <- function(mean0, sd0, min0, max0, n, x) {
        # 2. Remove missing values from new data
        x <- x[!is.na(x)]
        # 3. Handle case where prior sample size is one or less
        if (n <= 1) {
            # 4. If exactly one prior observation, include its mean
            if (n == 1) x <- c(x, mean0)
            # 5. Compute basic statistics from combined data
            mean1 <- mean(x)
            sd1 <- sd(x)
            min1 <- min(x)
            max1 <- max(x)
        } else {
            # 6. For larger n, update stats using existing means and sds
            m <- length(x)
            mean1 <- sum(mean0 * n, x) / (n + m)
            sd1 <- sqrt((n * (mean0 - mean1)^2 + sum((x - mean1)^2) + (n - 1) * sd0^2) / (m + n - 1))
            min1 <- min(min0, x)
            max1 <- max(max0, x)
        }
        # 7. Return updated mean, sd, min, max
        return(c(mean1, sd1, min1, max1))
    }

    # 8. Log-transform intensities and handle zeros
    ftrs.row[5:length(ftrs.row)] <- log10(ftrs.row[5:length(ftrs.row)] + 1)
    ftrs.row[ftrs.row == 0] <- NA
    # 9. Initialize existing row if absent
    if (length(existing.row) == 1) {
        existing.row <- rep(NA, 18)
        existing.row[6] <- ftrs.row[1]
    }

    # 10. Determine previous experiment counts
    n <- round(as.numeric(existing.row[7]) * as.numeric(existing.row[8])) # times found in previous experiments
    if (is.na(n)) n <- 0
    # 11. Count current experiment detections
    m <- sum(!is.na(chr.row[5:length(chr.row)])) # times found in current experiment

    # 12. Update total number of experiments and detection frequency
    existing.row[7] <- sum(as.numeric(existing.row[7]), length(ftrs.row) - 4, na.rm = TRUE)
    existing.row[8] <- (n + m) / as.numeric(existing.row[7])
    # 13. Update observed retention time range
    existing.row[9] <- min(as.numeric(existing.row[6]), as.numeric(existing.row[9]), ftrs.row[3], na.rm = TRUE)
    existing.row[10] <- max(as.numeric(existing.row[6]), as.numeric(existing.row[10]), ftrs.row[4], na.rm = TRUE)

    # 14. Merge retention time statistics with new observations
    this <- merge.new(as.numeric(existing.row[11]), as.numeric(existing.row[12]), as.numeric(existing.row[13]), as.numeric(existing.row[14]), n, chr.row[5:length(chr.row)])
    existing.row[11:14] <- this

    # 15. Merge intensity statistics with new data
    this <- merge.new(as.numeric(existing.row[15]), as.numeric(existing.row[16]), as.numeric(existing.row[17]), as.numeric(existing.row[18]), n, ftrs.row[5:length(ftrs.row)])
    existing.row[15:18] <- this

    # 16. Return updated feature characterization
    return(existing.row)
}
