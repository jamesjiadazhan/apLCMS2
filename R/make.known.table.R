#'
#' Given a table of known metabolites with original mass and charge information, and a table of allowable adducts, this function outputs a new table of potential features.
#'
#' @param metabolite.table A table of known metabolites. See the description of the object "metabolite.table" for details.
#' @param adduct.table A table of allowable adducts. See the description of the object "adduct.table" for details.
#' @param ion.mode Character. Either "+" or "-".
#'
#' @details For each allowable ion form, the function produces the m/z of every metabolite given to it. The output table follows the format that is required by the function semi.sup(), so that the user can directly use the table for semi supervised feature detection.
#'
#' @return A data frame containing the known metabolite ions. It contains 18 columns: "chemical_formula": the chemical formula if knonw; "HMDB_ID": HMDB ID if known; "KEGG_compound_ID": KEGG compound ID if known; "neutral.mass": the neutral mass if known; "ion.type": the ion form, such as H+, Na+, ..., if known; "m.z": m/z value, either theoretical for known metabolites, or mean observed value for unknown but previously found features; "Number_profiles_processed": the total number of LC/MS profiles that were used to build this database; "Percent_found": in what percentage was this feature found historically amount all data processed in building this database; "mz_min": the minimum m/z value observed for this feature; "mz_max": the maximum m/z value observed for this feature; "RT_mean": the mean retention time observed for this feature; "RT_sd": the standard deviation of retention time observed for this feature; "RT_min": the minimum retention time observed for this feature; "RT_max": the maximum retention time observed for this feature; "int_mean.log.": the mean log intensity observed for this feature; "int_sd.log.": the standard deviation of log intensity observed for this feature; "int_min.log.": the minimum log intensity observed for this feature; "int_max.log.": the maximum log intensity observed for this feature;
#'
#' @references Yu T, Park Y, Li S, Jones DP (2013) Hybrid feature detection and information accumulation using high-resolution LC-MS metabolomics data. J. Proteome Res. 12(3):1419-27.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @seealso metabolite.table, adduct.table, semi.sup
#'
#' @examples
#' data(metabolite.table)
#' data(adduct.table)
#' known.table.example<-make.known.table(metabolite.table[1001:1020,], adduct.table[1:4,])
#'
#' @keywords models
make.known.table <- function(metabolite.table, adduct.table, ion.mode = "+") {
    # 1. Sort metabolites by neutral mass
    metabolite.table <- metabolite.table[order(metabolite.table[, 4]), ]
    # 2. Ensure first three columns are character vectors
    for (i in 1:3) metabolite.table[, i] <- as.vector(metabolite.table[, i])

    # 3. Append empty columns for ion and statistics information
    metabolite.table <- cbind(metabolite.table, matrix(NA, nrow = nrow(metabolite.table), ncol = 14))
    # 4. Name the newly added columns
    colnames(metabolite.table)[5:ncol(metabolite.table)] <- c("ion.type", "m.z", "Number_profiles_processed", "Percent_found", "mz_min", "mz_max", "RT_mean", "RT_sd", "RT_min", "RT_max", "int_mean(log)", "int_sd(log)", "int_min(log)", "int_max(log)")

    # 5. Identify entries with explicit charge signs
    l <- nchar(metabolite.table[, 1])
    last.char <- substr(metabolite.table[, 1], l, l)
    sel <- which(last.char == "+" | last.char == "-")

    # 6. Initialize result table
    rm.first.row <- TRUE
    new.table <- matrix(NA, nrow = 1, ncol = ncol(metabolite.table))
    colnames(new.table) <- colnames(metabolite.table)

    # 7. Handle metabolites already carrying charge
    if (length(sel) > 0) {
        # 8. Add charged metabolites directly to new table
        new.table <- rbind(new.table, metabolite.table[sel, ])

        # 9. Extract numeric charge counts when present
        l <- nchar(new.table[, 1]) - 1
        num.charge <- substr(new.table[, 1], l, l)
        num.charge <- as.numeric(num.charge)

        # 10. Mark ion type as original and compute m/z
        new.table[, 5] <- "orig"
        new.table[, 6] <- new.table[, 4]
        new.table[which(!is.na(num.charge)), 6] <- new.table[which(!is.na(num.charge)), 6] / num.charge[which(!is.na(num.charge))]

        # 11. Remove charged entries from metabolite table
        metabolite.table <- metabolite.table[-sel, ]
    }

    # 12. Merge metabolites sharing identical neutral masses
    n <- 1
    m <- 2
    to.remove <- rep(0, nrow(metabolite.table))
    while (m <= nrow(metabolite.table)) {
        # 13. If masses nearly equal, combine records
        if (abs(metabolite.table[m, 4] - metabolite.table[n, 4]) < 1e-10) {
            if (metabolite.table[n, 1] != metabolite.table[m, 1]) metabolite.table[n, 1] <- paste(metabolite.table[n, 1], metabolite.table[m, 1], sep = "/")
            for (j in 2:3) metabolite.table[n, j] <- paste(metabolite.table[n, j], metabolite.table[m, j], sep = "/")
            to.remove[m] <- 1
            m <- m + 1
        } else {
            # 14. Move to next comparison pair
            n <- m
            m <- m + 1
        }
        cat("(", n, m, ")")
    }
    # 15. Drop merged rows
    if (sum(to.remove) > 0) metabolite.table <- metabolite.table[-which(to.remove == 1), ]

    # 16. Generate ion adducts for neutral metabolites
    neutral.table <- metabolite.table
    for (n in 1:nrow(adduct.table)) {
        # 17. Copy neutral table for each adduct type
        this <- neutral.table
        # 18. Assign ion type from adduct table
        this[, 5] <- adduct.table[n, 1]
        # 19. Compute resulting m/z based on charge and mass adjustment
        this[, 6] <- neutral.table[, 4] / adduct.table[n, 2] + adduct.table[n, 3]
        # 20. Append to result set
        new.table <- rbind(new.table, this)
    }
    # 21. Remove placeholder row if it exists
    if (rm.first.row) new.table <- new.table[-1, ]

    # 22. Filter by ion mode (positive or negative)
    to.remove <- rep(0, nrow(new.table))
    l <- nchar(new.table[, 1])
    last.char <- substr(new.table[, 1], l, l)
    if (ion.mode == "+") to.remove[which(last.char == "-")] <- 1
    if (ion.mode == "-") to.remove[which(last.char == "+")] <- 1
    if (sum(to.remove) > 0) new.table <- new.table[-which(to.remove == 1), ]

    # 23. Check for identical m/z values across ion forms
    new.table <- new.table[order(new.table[, 6]), ]
    n <- 1
    m <- 2
    to.remove <- rep(0, nrow(new.table))
    while (m <= nrow(new.table)) {
        # 24. Merge rows sharing the same m/z
        if (abs(new.table[m, 6] - new.table[n, 6]) < 1e-10) {
            if (new.table[n, 1] != new.table[m, 1]) new.table[n, 1] <- paste(new.table[n, 1], new.table[m, 1], sep = "/")
            for (j in 2:5) new.table[n, j] <- paste(new.table[n, j], new.table[m, j], sep = "/")
            to.remove[m] <- 1
            m <- m + 1
        } else {
            # 25. Advance to next unique mass
            n <- m
            m <- m + 1
        }
        cat("*(", n, m, ")")
    }
    # 26. Remove duplicated m/z entries
    if (sum(to.remove) > 0) new.table <- new.table[-which(to.remove == 1), ]

    # 27. Return completed known feature table
    return(new.table)
}
