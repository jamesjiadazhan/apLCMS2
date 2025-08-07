#'
#' This function takes the matrix output from proc.cdf() and generates a 3D plot of the data. It relies on the rgl library.
#'
#' @param prof The matrix output from the proc.cdf() function.
#' @param fill.holes A lot of peaks have missing values at some time points. If fill.holes is TRUE, the function will fill in the missing values by interpolation.
#' @param transform If the value is "sqrt", the values are transformed by taking square root. If "cuberoot", the values are transformed by taking cubic root.
#' @param time.lim The limit in retention time for the area of spectrum to be plotted. It should be either NA or a vector of two values: the lower limit and the upper limit.
#' @param mz.lim The limit in m/z value for the area of spectrum to be plotted. It should be either NA or a vector of two values: the lower limit and the upper limit.
#' @param box If a box should be drawn.
#' @param axes If the axes should be drawn.
#'
#' @details The function calls the rgl library. Spectrum values within the time.lim and mz.lim range is plotted in 3D.
#'
#' @return There is no return value from this function.
#'
#' @references http://rgl.neoscientists.org/about.shtml
#'
#' @author Tianwei Yu <tyu8@emory.edu>
#'
#' @examples
#' data(prof)
#' present.cdf.3d(prof[[2]],time.lim=c(250,400), mz.lim=c(400,500))
#'
present.cdf.3d <- function(prof, fill.holes = TRUE, transform = "none", time.lim = NA, mz.lim = NA, box = TRUE, axes = TRUE) {
    # 1. Helper to fill missing intensity values by linear interpolation
    fillholes <- function(times, b) {
        # 2. Order observations by time
        b <- b[order(b[, 1]), ]
        # 3. Determine all time points to cover
        newb <- times[times >= min(b[, 1]) & times <= max(b[, 1])]
        # 4. Initialize matrix with NA intensities
        newb <- cbind(newb, rep(NA, length(newb)))
        # 5. Insert observed intensities at known times
        newb[newb[, 1] %in% b[, 1], 2] <- b[, 2]
        # 6. Locate indices with existing values
        exist.b <- which(!is.na(newb[, 2]))
        # 7. Loop across consecutive observed segments
        for (k in 1:(length(exist.b) - 1)) {
            # 8. Check for gaps between observed points
            if (exist.b[k + 1] - exist.b[k] > 1) {
                # 9. Linearly interpolate intensities within gaps
                for (m in (exist.b[k] + 1):(exist.b[k + 1] - 1)) {
                    newb[m, 2] <- (newb[exist.b[k + 1], 2] - newb[exist.b[k], 2]) / (newb[exist.b[k + 1], 1] - newb[exist.b[k], 1]) * (newb[m, 1] - newb[exist.b[k], 1]) + newb[exist.b[k], 2]
                }
            }
        }
        # 10. Return filled matrix
        return(newb)
    }

    # 11. Load rgl library for 3D plotting
    library(rgl)
    # 12. Set initial window size for plot
    r3dDefaults$windowRect <- c(0, 50, 700, 700)
    # 13. Start with profile matrix
    a <- prof

    # 14. Restrict to specified time window if provided
    if (!is.na(time.lim[1])) {
        a <- a[a[, 2] >= time.lim[1] & a[, 2] <= time.lim[2], ]
    }
    # 15. Restrict to specified m/z range if provided
    if (!is.na(mz.lim[1])) {
        a <- a[a[, 1] >= mz.lim[1] & a[, 1] <= mz.lim[2], ]
    }

    # 16. Group by EIC labels and replace m/z with group medians
    a <- a[order(a[, 4]), ]
    n <- nrow(a)
    breaks <- c(0, which(a[1:(n - 1), 4] != a[2:n, 4]), n)
    for (i in 1:(length(breaks) - 1)) {
        a[(breaks[i] + 1):breaks[i + 1], 1] <- median(a[(breaks[i] + 1):breaks[i + 1], 1])
    }

    # 17. Apply optional intensity transformations
    if (transform == "sqrt") a[, 3] <- sqrt(a[, 3])
    if (transform == "cuberoot") a[, 3] <- a[, 3]^(1 / 3)

    # 18. Build sorted vector of time points for grid
    times <- unique(c(time.lim, a[, 2]))
    times <- times[order(times, na.last = FALSE)]
    # 19. Ensure time range covers limits
    if (!is.na(time.lim[1])) {
        if (min(times) > time.lim[1]) times <- c(time.lim[1], min(times) - 0.001, times)
        if (max(times) < time.lim[2]) times <- c(times, max(times) + 0.001, time.lim[2])
    }
    # 20. Drop any NA times
    times <- times[!is.na(times)]
    # 21. Gather unique m/z values for grid
    masses <- unique(a[, 1])
    masses <- masses[order(masses)]

    # 22. Optionally fill missing intensities in each EIC
    if (fill.holes) {
        new.a <- a[1, ]
        for (i in 1:length(masses)) {
            # 23. Extract signals for a single m/z group
            this <- a[a[, 1] == masses[i], ]
            if (!is.null(dim(this))) {
                this.labels <- unique(this[, 4])
                for (j in 1:length(this.labels)) {
                    # 24. Retrieve time-intensity pairs for current label
                    b <- this[this[, 4] == this.labels[j], 2:3]
                    if (!is.null(dim(b))) {
                        # 25. Interpolate missing time points
                        this.new <- fillholes(times, b)
                        # 26. Attach m/z and label back to filled data
                        this.new <- cbind(rep(this[1, 1], nrow(this.new)), this.new, rep(this.labels[j], nrow(this.new)))
                        # 27. Accumulate into result matrix
                        new.a <- rbind(new.a, this.new)
                    }
                }
            }
        }
        # 28. Replace original data with filled version
        a <- new.a[-1, ]
    }

    # 29. Prepare full m/z grid including small offsets for surface drawing
    masses <- unique(c(mz.lim, masses))
    masses <- masses[order(masses)]
    masses <- masses[!is.na(masses)]
    masses <- c(masses, masses - (1e-6), masses + (1e-10), masses + (1e-10 + 1e-6))
    masses <- masses[order(masses)]

    # 30. Initialize intensity matrix for surface
    z <- matrix(0, nrow = length(times), ncol = length(masses))

    # 31. Map m/z and time to matrix indices
    a <- a[order(a[, 1]), ]
    a[, 1] <- as.numeric(as.factor(rank(a[, 1])))
    a <- a[order(a[, 2]), ]
    a[, 2] <- as.numeric(as.factor(rank(a[, 2])))
    a <- a[order(a[, 1]), ]
    a[, 1] <- 4 * (a[, 1] - 1) + 2

    # 32. Fill intensity matrix for plotting
    for (i in 1:nrow(a)) {
        b <- a[i, ]
        z[b[2], b[1]] <- z[b[2], b[1] + 1] <- b[3]
    }

    # 33. Determine color scaling based on intensities
    zlim <- range(z)
    colorlut <- topo.colors(100)
    col <- colorlut[round(sqrt(z) / max(sqrt(z)) * 99) + 1]

    # 34. Render 3D surface of spectrum
    surface3d(times, masses, z, color = col)
    # 35. Define default plotting limits if not supplied
    if (is.na(time.lim[1])) time.lim <- range(times)
    if (is.na(mz.lim[1])) mz.lim <- range(masses)
    # 36. Add axes and box decorations as requested
    decorate3d(aspect = "iso", xlab = "time", ylab = "mz", zlab = "", xlim = time.lim, ylim = mz.lim, box = box, axes = axes)
    # 37. Ensure aspect ratio is equal in all directions
    aspect3d(1, 1, 1)
    # 38. Set initial viewing rotation
    par3d(userMatrix = rotationMatrix(90 * pi / 180, -2, 1, 1.5))
}
