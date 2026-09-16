# Right heart catheterization data: download and preparation
#
# Prepares the data analysed in section 7.3 of Orihara, Momozaki and Sugasawa
# (2025), "Impact of Right Heart Catheterization on Mortality with Confounder
# Selection". The analysis itself is in scripts/application_rhc.R, which reads
# the file this script writes.
#
# The data accompany Connors et al. (1996) and are distributed by the
# Vanderbilt Department of Biostatistics at https://hbiostat.org/data/ under
# their terms of use, which ask that the source be acknowledged. They are
# downloaded here and never committed: the download and the prepared copy both
# go to a cache directory outside the repository.
#
# Run it with
#   Rscript data-raw/rhc.R
# from the top of the package. Set DRBAYES_DATA_DIR to keep the cache
# somewhere permanent; the default is the session's temporary directory, so a
# fresh session downloads the file again.
#
# References
#   Connors, A. F., et al. (1996). The effectiveness of right heart
#     catheterization in the initial care of critically ill patients. JAMA,
#     276(11), 889-897.
#   Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#     Causal Inference via Posterior Coupling. arXiv:2506.04868.

rhc_url <- "https://hbiostat.org/data/repo/rhc.csv"

#' Directory the download and the prepared copy are cached in
#'
#' The cache lives outside the package so that a downloaded file can never be
#' picked up by R CMD build or committed by accident.
#'
#' @return Path to the cache directory, which is created if it does not exist.
#'   \code{DRBAYES_DATA_DIR} chooses it; without that it is a directory in the
#'   session's temporary space, so a fresh session downloads the file again.
rhc_cache_dir <- function() {
  dir <- Sys.getenv("DRBAYES_DATA_DIR", unset = "")
  if (!nzchar(dir)) {
    dir <- file.path(tempdir(), "drbayes-data")
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}


# ---------------------------------------------------------------------------
# Candidate confounders
# ---------------------------------------------------------------------------
#
# The candidates are the baseline variables of Connors et al. (1996): the
# admission diagnosis, the comorbidity history, the demographics and the
# physiological measurements taken on the first day. Everything recorded after
# admission is left out, because a variable measured after the catheter was
# placed can be affected by it and conditioning on it would block part of the
# effect being estimated. So are the identifiers (ptid), the administrative
# dates (sadmdte, dschdte, dthdte, lstctdte) and the outcomes (death, dth30,
# t3d30, surv2md1 is a baseline prognosis and is kept).
#
# Three variables in the file are dropped for missingness rather than for
# timing: cat2, the secondary disease category, is missing for 79 percent of
# patients, adld3p for 75 percent and urin1 for 53 percent. Section 7.2 selects
# confounders from among the candidates, so a candidate that would force three
# quarters of the sample to be discarded is not a useful one to offer.
#
# Section 7.3 of the paper reports 48 candidates expanding to 58 columns once
# the multi-level factors are dummy coded, and refers to Harada and Taguri
# (2025) for the list rather than restating it. The list below is the standard
# Connors et al. set and gives 50 candidates and 65 columns. Substitute the
# published list here to match the paper's counts exactly.

# Multi-level and two-level factors. Their number of dummy columns is what
# separates the candidate count from the design matrix column count.
rhc_factors <- c(
  "cat1",      # primary disease category, 9 levels
  "ca",        # cancer: no, localised, metastatic
  "sex", "race", "income", "ninsclas",
  "dnr1",      # do-not-resuscitate order on day 1
  # Admission diagnosis categories, one indicator each.
  "resp", "card", "neuro", "gastr", "renal", "meta", "hema", "seps",
  "trauma", "ortho"
)

# Comorbidities recorded as 0/1 in the file, so they need no conversion.
rhc_history <- c(
  "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx", "renalhx",
  "liverhx", "gibledhx", "malighx", "immunhx", "transhx", "amihx"
)

rhc_numeric <- c(
  "age", "edu",
  "surv2md1",   # SUPPORT model estimate of 2-month survival at admission
  "das2d3pc",   # Duke activity status index
  "aps1", "scoma1", "meanbp1", "wblc1", "hrt1", "resp1", "temp1", "pafi1",
  "alb1", "hema1", "bili1", "crea1", "sod1", "pot1", "paco21", "ph1",
  "wtkilo1"
)

rhc_candidates <- c(rhc_factors, rhc_history, rhc_numeric)


#' Download the raw file, or reuse the cached copy
#'
#' @param dir Directory to cache the download in.
#' @param url Address of the CSV file.
#'
#' @return Path to the downloaded file.
rhc_download <- function(dir = rhc_cache_dir(), url = rhc_url) {
  destination <- file.path(dir, "rhc.csv")

  if (file.exists(destination)) {
    message("Using the cached copy at ", destination, ".")
    return(destination)
  }

  message("Downloading ", url, " to ", destination, ".")
  status <- tryCatch(
    utils::download.file(url, destination, mode = "wb", quiet = TRUE),
    error = function(e) e)

  if (inherits(status, "error") || !file.exists(destination)) {
    stop("Could not download ", url, ". Check the network connection, or ",
         "fetch the file by hand and save it as ", destination,
         ", after which this script will use it.", call. = FALSE)
  }

  destination
}


#' Turn the raw file into the analysis data frame
#'
#' @param path Path to the downloaded CSV file.
#'
#' @return A data frame holding the treatment \code{rhc}, the outcome
#'   \code{death30} and the candidate confounders, and nothing else, so that
#'   the candidates can be recovered with
#'   \code{setdiff(names(x), c("rhc", "death30"))}.
rhc_prepare <- function(path) {
  raw <- utils::read.csv(path, stringsAsFactors = FALSE)

  expected <- c("swang1", "dth30", rhc_candidates)
  missing_columns <- setdiff(expected, names(raw))
  if (length(missing_columns) > 0) {
    stop("The downloaded file does not have the columns this script expects: ",
         toString(missing_columns), ". The distributed file has changed; ",
         "compare it against the codebook at https://hbiostat.org/data/ and ",
         "update the candidate lists at the top of this script.",
         call. = FALSE)
  }

  prepared <- raw[, rhc_candidates, drop = FALSE]
  for (variable in rhc_factors) {
    prepared[[variable]] <- factor(prepared[[variable]])
  }

  # The exposure of section 7.1 is receipt of a right heart catheter within
  # the first 24 hours, and the outcome is death within 30 days of admission.
  # Both are stored as text, so they are compared against the level that
  # counts as the event rather than coerced.
  prepared$rhc      <- as.integer(raw$swang1 == "RHC")
  prepared$death30  <- as.integer(raw$dth30 == "Yes")

  prepared <- prepared[, c("rhc", "death30", rhc_candidates), drop = FALSE]

  incomplete <- !stats::complete.cases(prepared)
  if (any(incomplete)) {
    message(sum(incomplete), " of ", nrow(prepared),
            " patients have a missing value among the analysis variables ",
            "and are dropped.")
    prepared <- prepared[!incomplete, , drop = FALSE]
  }

  rownames(prepared) <- NULL
  prepared
}


#' Describe the prepared data the way section 7 does
#'
#' @param data The data frame returned by \code{rhc_prepare()}.
#'
#' @return \code{data}, invisibly.
rhc_report <- function(data) {
  candidates <- setdiff(names(data), c("rhc", "death30"))
  columns    <- ncol(stats::model.matrix(~ ., data[, candidates])) - 1L

  treated <- data$rhc == 1
  message("Patients:            ", nrow(data),
          " (", sum(treated), " catheterised, ", sum(!treated), " not)")
  message("30-day mortality:    ",
          format(mean(data$death30[treated]), digits = 3), " catheterised, ",
          format(mean(data$death30[!treated]), digits = 3), " not, ",
          "crude risk difference ",
          format(mean(data$death30[treated]) - mean(data$death30[!treated]),
                 digits = 3))
  message("Candidates:          ", length(candidates),
          ", expanding to ", columns, " design matrix columns")
  message("Paper section 7.3:   48 candidates, 58 columns")

  invisible(data)
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

rhc_output <- file.path(rhc_cache_dir(), "rhc.rds")

rhc <- rhc_prepare(rhc_download())
rhc_report(rhc)

saveRDS(rhc, rhc_output)
message("Wrote ", rhc_output, ".")
message("scripts/application_rhc.R reads this file.")
