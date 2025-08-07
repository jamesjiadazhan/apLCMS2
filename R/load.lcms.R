#' Loading LC/MS data.
#'
#' This is an internal function. It loads LC/MS data into memory.
#'
#' @param filename The CDF file name.
#'
#' @details The function uses functionality provided by the mzR package from Bioconductor.
#'
#' @return A list is returned.
#'   \item{masses}{ The vector of m/z values. }
#'   \item{labels}{ The vector of retention times. }
#'   \item{intensi}{ The vector of intensity values. }
#'   \item{times}{ The vector of unique time points. }
#'
#' @references Bioinformatics. 25(15):1930-36. BMC Bioinformatics. 11:559.
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @keywords models
load.lcms <- function(filename) {
    # 1. Load mzR package to read mass spectrometry files
    library(mzR)
    # 2. Split filename to check file extension
    splitname <- strsplit(filename, "\\.")[[1]]
    # 3. If file is a CDF file
    if (tolower(splitname[length(splitname)]) == "cdf") {
        # 4. Open using netCDF backend for compatibility
        mz.conn <- openMSfile(filename, backend = "netCDF")
    } else {
        # 5. Otherwise open using default backend
        mz.conn <- openMSfile(filename)
    }

    # 6. Initialize container for m/z values
    masses <- NULL
    # 7. Initialize container for intensity values
    intensi <- NULL
    # 8. Initialize container for retention time labels
    labels <- NULL
    # 9. Extract retention time for each scan
    b <- header(mz.conn)$retentionTime

    # 10. Create segments of 200 scans for batch processing
    segs <- seq(0, length(b), by = 200)
    # 11. Add last segment if scans not divisible by 200
    if ((length(b) %% 200) != 0) segs <- c(segs, length(b))

    # 12. Loop over each segment
    for (n in 2:length(segs)) {
        # 13. Retrieve peak matrices for scans in segment
        a <- peaks(mz.conn, scans = (segs[n - 1] + 1):segs[n])

        # 14. Temporary storage for m/z in current segment
        this.masses <- NULL
        # 15. Temporary storage for intensities in current segment
        this.intensi <- NULL
        # 16. Temporary storage for retention labels in current segment
        this.labels <- NULL

        # 17. Iterate through each scan within segment
        for (i in 1:length(a)) {
            # 18. Extract m/z-intensity matrix for current scan
            this.a <- a[[i]]
            # 19. Proceed only if peaks exist in scan
            if (!is.null(nrow(this.a))) {
                # 20. Remove peaks with near-zero intensity
                this.a <- this.a[this.a[, 2] > 1e-10, ]
                # 21. Ensure data is in matrix form
                if (is.null(nrow(this.a))) this.a <- matrix(this.a, nrow = 1)

                # 22. Append m/z values
                this.masses <- c(this.masses, this.a[, 1])
                # 23. Append intensity values
                this.intensi <- c(this.intensi, this.a[, 2])
                # 24. Record corresponding retention times
                this.labels <- c(this.labels, rep(b[segs[n - 1] + i], nrow(this.a)))
            } else {
                # 25. Mark scans with no data as NA
                b[segs[n - 1] + i] <- NA
            }
        }

        # 26. Accumulate m/z values across segments
        masses <- c(masses, this.masses)
        # 27. Accumulate intensities across segments
        intensi <- c(intensi, this.intensi)
        # 28. Accumulate retention time labels across segments
        labels <- c(labels, this.labels)
    }
    # 29. Remove scans flagged as NA to get valid times
    times <- b[!is.na(b)]
    # 30. Close the mass spectrometry file connection
    close(mz.conn)

    # 31. Return assembled data
    return(list(masses = masses, labels = labels, intensi = intensi, times = times))
}
