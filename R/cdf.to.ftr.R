#' Convert a number of cdf files in the same directory to a feature table
#'
#' This is a wrapper function, which calls four other functions to convert a number of cdf files to a feature table. All cdf files to be processed must be in a single folder.
#'
#' @param folder The folder where all CDF files to be processed are located. For example "C:/CDF/this_experiment"
#' @param file.pattern The pattern in the names of the files to be processed. The default is ".cdf". Other formats supported by mzR package can also be used, e.g. "mzML" etc.
#' @param n.nodes The number of CPU cores to be used through doSNOW.
#' @param min.exp If a feature is to be included in the final feature table, it must be present in at least this number of spectra.
#' @param min.pres This is a parameter of thr run filter, to be passed to the function proc.cdf(). Please see the help for proc.cdf() for details.
#' @param min.run This is a parameter of thr run filter, to be passed to the function proc.cdf(). Please see the help for proc.cdf() for details.
#' @param mz.tol The user can provide the m/z tolerance level for peak identification. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level. Please see the help for proc.cdf() for details.
#' @param baseline.correct.noise.percentile The perenctile of signal strength of those EIC that don't pass the run filter, to be used as the baseline threshold of signal strength. This parameter is passed to proc.cdf()
#' @param shape.model The mathematical model for the shape of a peak. There are two choices - bi-Gaussian and Gaussian. When the peaks are asymmetric, the bi-Gaussian is better. The default is bi-Gaussian.
#' @param BIC.factor the factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1, models with more peaks are penalized more.
#' @param baseline.correct This is a parameter in peak detection. After grouping the observations, the highest observation in each group is found. If the highest is lower than this value, the entire group will be deleted. The default value is NA, which allows the program to search for the cutoff level. Please see the help for proc.cdf() for details.
#' @param peak.estim.method the bi-Gaussian peak parameter estimation method, to be passed to subroutine prof.to.features. Two possible values: moment and EM.
#' @param min.bw The minimum bandwidth in the smoother in prof.to.features(). Please see the help file for prof.to.features() for details.
#' @param max.bw The maximum bandwidth in the smoother in prof.to.features(). Please see the help file for prof.to.features() for details.
#' @param sd.cut A parameter for the prof.to.features() function. A vector of two. Features with standard deviation outside the range defined by the two numbers are eliminated.
#' @param sigma.ratio.lim A parameter for the prof.to.features() function. A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation and the righ-standard deviation of the bi-Gaussian fuction used to fit the data.
#' @param component.eliminate In fitting mixture of bi-Gaussian (or Gaussian) model of an EIC, when a component accounts for a proportion of intensities less than this value, the component will be ignored.
#' @param moment.power The power parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param subs If not all the CDF files in the folder are to be processed, the user can define a subset using this parameter. For example, subs=15:30, or subs=c(2,4,6,8)
#' @param align.chr.tol The user can provide the elution time tolerance level to override the program's selection. This value is in the same unit as the elution time, normaly seconds. Please see the help for match.time() for details.
#' @param align.mz.tol The user can provide the m/z tolerance level for peak alignment to override the program's selection. This value is expressed as the percentage of the m/z value. This value, multiplied by the m/z value, becomes the cutoff level.Please see the help for feature.align() for details.
#' @param max.align.mz.diff As the m/z tolerance in alignment is expressed in relative terms (ppm), it may not be suitable when the m/z range is wide. This parameter limits the tolerance in absolute terms. It mostly influences feature matching in higher m/z range.
#' @param pre.process Logical. If true, the program will not perform time correction and alignment. It will only generate peak tables for each spectra and save the files. It allows manually dividing the task to multiple machines.
#' @param recover.mz.range A parameter of the recover.weaker() function. The m/z around the feature m/z to search for observations. The default value is NA, in which case 1.5 times the m/z tolerance in the aligned object will be used.
#' @param recover.chr.range A parameter of the recover.weaker() function. The retention time around the feature retention time to search for observations. The default value is NA, in which case 0.5 times the retention time tolerance in the aligned object will be used.
#' @param use.observed.range A parameter of the recover.weaker() function. If the value is TRUE, the actual range of the observed locations of the feature in all the spectra will be used.
#' @param recover.min.count The minimum time point count for a series of point in the EIC for it to be considered a true feature.
#' @param intensity.weighted Whether to weight the local density by signal intensities in the initial peak detection.
#'
#' @details The wrapper function calls five other functions to perform the feature table generation. Every spectrum (cdf file) first goes through proc.cdf() and prof.to.feature() to generate a spectrum-level peak table. The eluction time correction is done by match.time(). Then the peaks are aligned across spectra by feature.align(). For features deteced in a portion of the spectra, weaker signals in other spectra are recovered by recover.weaker().
#' From version 4, the parameter mz.tol can no longer be NA. This is to allow the program better process data other than FTLCMS. It is recommended that the user use the machine's claimed accuracy. For FTMS, 1e-5 is recommended.
#'
#' @return A list is returned.
#'   \item{features}{ A list object, each component of which being the peak table from a single spectrum.}
#'   \item{features2}{A list object, each component of which being the peak table from a single spectrum, after elution time correction.}
#'   \item{aligned.ftrs}{ Feature table BEFORE weak signal recovery.}
#'   \item{final.ftrs}{Feature table after weak signal recovery. This is the end product of the function.}
#'   \item{pk.times}{ Table of feature elution time BEFORE weak signal recovery.}
#'   \item{final.times}{Table of feature elution time after weak signal recovery.}
#'   \item{mz.tol}{The input mz.tol value by the user.}
#'   \item{align.mz.tol}{The m/z tolerance level in the alignment across spectra, either input from the user or automatically selected when the user input is NA.}
#'   \item{align.chr.tol}{The retention time tolerance level in the alignment across spectra, either input from the user or automatically selected when the user input is NA.}
#'
#' @author Tianwei Yu <tyu8@sph.emory.edu>
#'
#' @seealso proc.cdf, prof.to.feature, adjust.time, feature.align, recover.weaker
#'
#' @keywords models
# 1. Top-level pipeline: per-file preprocessing -> feature extraction -> (optional) time correction -> feature alignment -> recovery -> return tables
cdf.to.ftr <- function(folder, file.pattern = ".mzXML", n.nodes = 4, min.exp = 2, min.pres = 0.5, min.run = 12, mz.tol = 1e-5, baseline.correct.noise.percentile = 0.05, shape.model = "bi-Gaussian", BIC.factor = 2, baseline.correct = 0, peak.estim.method = "moment", min.bw = NA, max.bw = NA, sd.cut = c(0.01, 500), sigma.ratio.lim = c(0.01, 100), component.eliminate = 0.01, moment.power = 1, subs = NULL, align.mz.tol = NA, align.chr.tol = NA, max.align.mz.diff = 0.01, pre.process = FALSE, recover.mz.range = NA, recover.chr.range = NA, use.observed.range = TRUE, recover.min.count = 3, intensity.weighted = FALSE) {

    # ---- Parallel plan & progress setup (portable) ----
    if (!requireNamespace("future", quietly = TRUE)) install.packages("future")
    if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
    if (!requireNamespace("progressr", quietly = TRUE)) install.packages("progressr")

    library(future)
    library(future.apply)
    library(progressr)
    # library(doParallel)

    # plan: change to cluster/SLURM via future.batchtools later if you want
    plan(multisession, workers = n.nodes)
    handlers(global = TRUE)                # choose UI: handler_cli(), handler_txtprogressbar(), etc.
    options(progressr.enable = TRUE)

    # 2. Load dependencies; set working folder; enumerate candidate files and (optionally) subset
    library(mzR)

    # set up working directory (where apLCMS reads mzML or mzXML files)
    setwd(folder)

    # save the current function setting
    settings_list <- as.list(environment())
    settings_list <- as.data.frame(settings_list)
    write.table(settings_list, file = file.path(folder, "cdf.to.ftr.settings.txt"), row.names = FALSE)

    # prepare for the log file that may be used for debugging
    logfile <- file.path(folder, paste0("cdf.to.ftr.runlog_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".txt"))
    
    # open ONE connection and use it for both stdout and messages
    log_con <- file(logfile, open = "wt")  # or "at" if you truly want to append across runs
    sink(log_con, split = TRUE)            # stdout
    sink(log_con, type = "message", split = TRUE)        # messages/warnings/errors
    
    # ensure sinks and connection are restored/closed on exit
    on.exit({
      if (sink.number(type = "message") > 0) sink(type = "message")
      while (sink.number() > 0) sink()
      try(close(log_con), silent = TRUE)
    }, add = TRUE)

    
    files <- dir(pattern = file.pattern, ignore.case = TRUE)
    files <- files[order(files)]
    if (!is.null(subs)) {
        if (!is.na(subs[1])) files <- files[subs]
    }

    ###############################################################################################

    # 3. Prepare per-file work identifiers and suffixes for caching; create error bin for problem files
    dir.create("error_files")
    message("***************************** Feature detection and chromatogram building for each sample *****************************")
    suf.prof <- paste(min.pres, min.run, mz.tol, baseline.correct, sep = "_")
    suf <- paste(suf.prof, shape.model, sd.cut[1], sd.cut[2], component.eliminate, moment.power, sep = "_")
    if (shape.model == "bi-Gaussian") suf <- paste(suf, sigma.ratio.lim[1], sigma.ratio.lim[2], sep = "_")

    to.do <- paste(matrix(unlist(strsplit(tolower(files), "\\.")), nrow = 2)[1, ], suf, min.bw, max.bw, ".feature", sep = "_")
    to.do <- which(!(to.do %in% dir()))
    message(c("number of files to process: ", length(to.do)))

    # 4. Parallel block: for each file, produce raw profile via proc.cdf(), then feature table via prof.to.features(); cache outputs
    if (length(to.do) > 0) {

        # helper: do one file index j (keeps body readable and re-usable)
        .process_one <- function(j, files, suf, suf.prof,
                                 min.pres, min.run, mz.tol,
                                 baseline.correct, baseline.correct.noise.percentile,
                                 intensity.weighted, min.bw, max.bw, sd.cut,
                                 shape.model, peak.estim.method, component.eliminate,
                                 moment.power, BIC.factor) {
    
            # ensure required pkg is visible on workers
            if (!"apLCMS" %in% loadedNamespaces()) {
                suppressPackageStartupMessages(require(apLCMS, quietly = TRUE, character.only = TRUE))
            }
    
            feature.name <- paste(strsplit(tolower(files[j]), "\\.")[[1]][1], suf,     min.bw, max.bw, ".feature", sep = "_")
            profile.name <- paste(strsplit(tolower(files[j]), "\\.")[[1]][1], suf.prof,         ".profile",       sep = "_")
    
            # 5. Build raw profile (proc.cdf)
            ok_prof <- TRUE
            res_prof <- try({
                this.prof <- proc.cdf(
                    files[j],
                    min.pres = min.pres, min.run = min.run, tol = mz.tol,
                    baseline.correct = baseline.correct,
                    baseline.correct.noise.percentile = baseline.correct.noise.percentile,
                    do.plot = FALSE, intensity.weighted = intensity.weighted
                )
                # cache the output
                save(this.prof, file = profile.name)
                TRUE
            }, silent = TRUE)
    
            if (inherits(res_prof, "try-error")) {
                ok_prof <- FALSE
                # quarantine problem file
                try(file.copy(from = files[j], to = "error_files"), silent = TRUE)
                try(file.remove(files[j]), silent = TRUE)
            }
    
            # 6. If profile succeeded, run prof.to.features
            if (ok_prof) {
                res_feat <- try({
                    this.feature <- prof.to.features(
                        this.prof, min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut,
                        shape.model = shape.model, estim.method = peak.estim.method,
                        do.plot = FALSE, component.eliminate = component.eliminate,
                        power = moment.power, BIC.factor = BIC.factor
                    )
                    # cache the output
                    save(this.feature, file = feature.name)
                    TRUE
                }, silent = TRUE)
    
                if (inherits(res_feat, "try-error")) {
                    try(file.copy(from = files[j], to = "error_files"), silent = TRUE)
                    try(file.remove(files[j]), silent = TRUE)
                }
            }
    
            # return a tiny status list (keeps memory low)
            list(file = files[j],
                 profile_saved = ok_prof && !inherits(res_prof, "try-error"),
                 feature_saved = ok_prof && !inherits(res_feat, "try-error"))
        }

        # ---- Progress + parallel apply over file indices ----
        with_progress({
            p <- progressor(steps = length(to.do))
    
            invisible(future_lapply(
                to.do,
                function(j) {
                    out <- .process_one(
                        j = j, files = files, suf = suf, suf.prof = suf.prof,
                        min.pres = min.pres, min.run = min.run, mz.tol = mz.tol,
                        baseline.correct = baseline.correct,
                        baseline.correct.noise.percentile = baseline.correct.noise.percentile,
                        intensity.weighted = intensity.weighted,
                        min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut,
                        shape.model = shape.model, peak.estim.method = peak.estim.method,
                        component.eliminate = component.eliminate,
                        moment.power = moment.power, BIC.factor = BIC.factor
                    )
                    p(sprintf("Processed %s", out$file))
                    NULL
                }
            ))
        })
    }
    
    # 8. Restrict to files with available outputs; load all feature tables into a list
    message("Loading the extracted feature tables from each sample into a list")
    all.files <- dir()
    sel <- which(files %in% all.files)
    files <- files[sel]

    features <- new("list")
    for (i in 1:length(files)) {
        feature.name <- paste(strsplit(tolower(files[i]), "\\.")[[1]][1], suf, min.bw, max.bw, ".feature", sep = "_")
        cat(feature.name, "\n")
        load(feature.name)
        features[[i]] <- this.feature
    }

    gc()

    if (!pre.process) {
        ###############################################################################################
        message("****************************** time correction ***************************************")
        # 9. Time correction across profiles using adjust.time(); cache the adjusted feature lists
        suf <- paste(suf, align.mz.tol, align.chr.tol, subs[1], subs[length(subs)], sep = "_")
        time.correction.name <- paste("time_correct_done_", suf, ".bin", sep = "")

        all.files <- dir()
        is.done <- all.files[which(all.files == time.correction.name)]

        if (length(is.done) == 0) {
            if (!requireNamespace("future", quietly = TRUE)) install.packages("future")
            library(future)
        
            res <- value(future({
                # Ensure package is available on the worker
                suppressPackageStartupMessages(require(apLCMS, quietly = TRUE, character.only = TRUE))
                st <- system.time({
                    at <- adjust.time(
                        features,
                        mz.tol = align.mz.tol,
                        chr.tol = align.chr.tol,
                        find.tol.max.d = 10 * mz.tol,
                        max.align.mz.diff = max.align.mz.diff
                    )
                })
                list(f2 = at, elapsed = as.vector(st)[1])
            }))
        
            message("***** correcting time, CPU time (seconds) ", res$elapsed)
            f2 <- res$f2
            save(f2, file = time.correction.name)
        } else {
            load(time.correction.name)
            message("Previously done time corrected feature table is loaded")
        }

        # clear memory
        gc()
        
        ###############################################################################################
        message("****************************  aligning features **************************************")
        # 10. Align features across profiles with feature.align(); cache alignment
        suf <- paste(suf, min.exp, sep = "_")
        feature.alignment.name <- paste("aligned_done_", suf, ".bin", sep = "")
        all.files <- dir()
        is.done <- all.files[which(all.files == feature.alignment.name)]
        
        if (length(is.done) == 0) {
            message(c("***** aligning features, CPU time (seconds): ", as.vector(system.time(aligned <- feature.align(f2, min.exp = min.exp, mz.tol = align.mz.tol, chr.tol = align.chr.tol, find.tol.max.d = 10 * mz.tol, max.align.mz.diff = max.align.mz.diff, n.nodes = n.nodes)))[1]))
            save(aligned, file = feature.alignment.name)

        } else {
            load(feature.alignment.name)
            message("Previously done aligned feature table is loaded")
        }

        # clear memory
        gc()

        ###############################################################################################
        message("**************************** recovering weaker signals *******************************")
        # 11. Weak-signal recovery around aligned features; cache per-file recoveries
        suf <- paste(suf, recover.mz.range, recover.chr.range, use.observed.range, sep = "_")

        worklist <- paste(matrix(unlist(strsplit(tolower(files), "\\.")), nrow = 2)[1, ], suf, ".recover", sep = "_")
        to.do <- which(!(worklist %in% dir()))
        grps <- round(seq(0, length(to.do), length = n.nodes + 1))
        grps <- unique(grps)

        message(c("number of files to process: ", length(to.do)))

        # ---- Recover weaker signals (future.apply + progressr for modern parellelization and debugging) ----
        if (length(to.do) > 0) {

            # To reduce repeated serialization, bind big objects to locals once
            aligned_ftrs  <- aligned$aligned.ftrs
            aligned_times <- aligned$pk.times
            align_mz_tol  <- aligned$mz.tol
            align_chr_tol <- aligned$chr.tol

            # helper: process a single index j
            .recover_one <- function(j, files, suf,
                                     aligned_ftrs, aligned_times, align_mz_tol, align_chr_tol,
                                     features, f2,
                                     recover_mz_range, recover_chr_range, use_observed_range,
                                     mz_tol, min_bw, max_bw, recover_min_count) {
        
                # make sure apLCMS is available on workers
                if (!"apLCMS" %in% loadedNamespaces()) {
                    suppressPackageStartupMessages(require(apLCMS, quietly = TRUE, character.only = TRUE))
                }
        
                feature.recover.name <- paste(strsplit(tolower(files[j]), "\\.")[[1]][1], suf, ".recover", sep = "_")
                cat(feature.recover.name, "\n")
        
                this.recovered <- recover.weaker(
                    filename = files[j], loc = j,
                    aligned.ftrs = aligned_ftrs, pk.times = aligned_times,
                    align.mz.tol = align_mz_tol, align.chr.tol = align_chr_tol,
                    this.f1 = features[[j]], this.f2 = f2[[j]],
                    mz.range = recover_mz_range, chr.range = recover_chr_range,
                    use.observed.range = use_observed_range,
                    orig.tol = mz_tol, min.bw = min_bw, max.bw = max_bw,
                    bandwidth = .5, recover.min.count = recover_min_count
                )
        
                save(this.recovered, file = feature.recover.name)
                feature.recover.name
            }

            # initialize the parellel processing of recover.weaker function
            with_progress({
                p <- progressor(steps = length(to.do))

                invisible(future_lapply(
                    to.do,
                    function(j) {
                        out_name <- .recover_one(
                            j = j, files = files, suf = suf,
                            aligned_ftrs = aligned_ftrs, aligned_times = aligned_times,
                            align_mz_tol = align_mz_tol, align_chr_tol = align_chr_tol,
                            features = features, f2 = f2,
                            recover_mz_range = recover.mz.range, recover_chr_range = recover.chr.range,
                            use_observed_range = use.observed.range,
                            mz_tol = mz.tol, min_bw = min.bw, max_bw = max.bw,
                            recover_min_count = recover.min.count
                        )
                        p(sprintf("Processed %s", out_name))
                        NULL
                    }
                ))
            })

            # clear memory
            gc()
            
        }
                      
        # 12. Build final alignment object by injecting recovered intensities/times; collect outputs for return
        new.aligned <- aligned

        # initial assembly of return list object consolidating aligned and gap filled feature tables

        ## update the feature aligned feature table's column names
        colnames(aligned$aligned.ftrs) <- c("mz", "time", "mz.min", "mz.max", files)
        colnames(aligned$pk.times) <- c("mz", "time", "mz.min", "mz.max", files)

        ## initate a new list
        rec <- new("list")

        ## add data from feature aligned feature table to the list
        rec$aligned.ftrs <- aligned$aligned.ftrs
        rec$pk.times <- aligned$pk.times

        ## remove aligned to save memory
        rm(aligned)
        gc()

        # load all gap filled features
        for (i in 1:length(files)) {
            feature.recover.name <- paste(strsplit(tolower(files[i]), "\\.")[[1]][1], suf, ".recover", sep = "_")
            load(feature.recover.name)
            cat(feature.recover.name, "\n")

            # add the gap filled features to the aligned feature table
            new.aligned$aligned.ftrs[, i + 4] <- this.recovered$this.ftrs
            new.aligned$pk.times[, i + 4] <- this.recovered$this.times
            new.aligned$features[[i]] <- this.recovered$this.f1
            new.aligned$f2[[i]] <- this.recovered$this.f2
        }

        # clear memory
        gc()

        #################################################################################################
        message("Final assembly of return list object consolidating aligned and gap filled feature tables")
        # 13. Final assembly of return list object consolidating aligned and gap filled feature tables

        ## update column names
        colnames(new.aligned$aligned.ftrs) <- c("mz", "time", "mz.min", "mz.max", files)
        colnames(new.aligned$pk.times) <- c("mz", "time", "mz.min", "mz.max", files)

        ## add gap filled data into the list
        rec$features <- new.aligned$features
        rec$features2 <- new.aligned$f2
        rec$final.ftrs <- new.aligned$aligned.ftrs
        rec$final.times <- new.aligned$pk.times
        rec$align.mz.tol <- new.aligned$mz.tol
        rec$align.chr.tol <- new.aligned$chr.tol
        rec$mz.tol <- mz.tol

        ## return the list. final.ftrs is the feature table with both feature alignment and gap filling. aligned.ftrs is the feature table with only feature alignment.
        return(rec)
    }
}
