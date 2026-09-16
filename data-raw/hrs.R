# Health and Retirement Study data: preparation recipe
#
# Prepares the cohort analysed in section 6 of Orihara, Momozaki and Sugasawa
# (2025), "Effects of Antihypertensive Treatment on Dementia and Sensitivity
# Analysis".
#
# UNLIKE data-raw/rhc.R, THIS SCRIPT DOWNLOADS NOTHING. The Health and
# Retirement Study is free but not public: every file needs a registered
# account, and the genetic file needs a further agreement. There is no address
# a script can fetch from, so what follows is the derivation written out as
# runnable code, to be pointed at files you have already obtained.
#
# What to obtain, from https://hrsdata.isr.umich.edu/data-products/ :
#
#   1. RAND HRS Longitudinal File. One row per respondent, one set of columns
#      per wave, cleaned and named consistently across waves. Nearly every
#      variable below comes from here. Register, accept the conditions of use,
#      and download the R or Stata version.
#   2. Core interview physical measures, or the RAND HRS Fat File for the same
#      waves. Blood pressure is measured in the enhanced face-to-face module,
#      which the RAND Longitudinal File does not carry. Only a rotating half
#      sample is measured at each wave, which is why the cohort below is much
#      smaller than the number of respondents aged 65 and over.
#   3. APOE genotype. Distributed with the HRS genetic data under a separate
#      agreement, not with the core files. Section 6 lists the apolipoprotein
#      E epsilon-4 allele as a confounder, so the analysis cannot be run
#      without it; if you cannot obtain it, drop it from `hrs_confounders`
#      below and say so when reporting the result.
#
# Cite the study as the conditions of use require: the HRS is sponsored by the
# National Institute on Aging (grant NIA U01AG009740) and conducted by the
# University of Michigan.
#
# Run it with
#   HRS_DATA_DIR=/path/to/your/hrs/files Rscript data-raw/hrs.R
# from the top of the package. Without HRS_DATA_DIR the script defines its
# functions, prints the recipe and stops, which is all it can do.
#
# References
#   Bennett, E. E., et al. (2025), the analysis section 6 takes its confounder
#     list from.
#   Orihara, S., Momozaki, T., & Sugasawa, S. (2025). Bayesian Doubly Robust
#     Causal Inference via Posterior Coupling. arXiv:2506.04868.


# ---------------------------------------------------------------------------
# Variable names
# ---------------------------------------------------------------------------
#
# RAND HRS names a respondent-level variable r<wave><concept>, so age at wave
# 10 is r10agey_e, and a time-invariant one ra<concept>, so ragender. Waves are
# biennial from 1992: wave 10 is 2010, 11 is 2012, 12 is 2014, 13 is 2016.

#' Wave number of a survey year
#'
#' @param year Survey year, for example 2010.
#'
#' @return The wave number RAND HRS uses in its variable names, so 10 for 2010.
hrs_wave <- function(year) {
  if (!year %in% seq(1992, 2020, by = 2)) {
    stop("HRS waves are biennial from 1992; ", year, " is not one of them.",
         call. = FALSE)
  }
  (year - 1992) / 2 + 1
}

#' RAND HRS variable names for one wave
#'
#' @param year Survey year, for example 2010.
#'
#' @return A named character vector mapping the concepts used below to the
#'   RAND HRS column names for that wave.
hrs_rand_names <- function(year) {
  w <- hrs_wave(year)
  r <- function(concept) paste0("r", w, concept)

  c(
    interviewed  = r("iwstat"),   # 1 alive and interviewed, 5 died
    age          = r("agey_e"),
    hypertension = r("hibpe"),    # ever told had high blood pressure
    bp_medication = r("rxhibp"),  # currently takes medication for it
    dementia     = r("demene"),
    alzheimers   = r("alzhe"),
    memory       = r("memrye"),   # any memory-related disease
    bmi          = r("bmi"),
    smoking      = r("smoken"),   # smokes now
    depression   = r("cesd"),     # CES-D score, 0 to 8
    heart        = r("hearte"),   # ever had a heart condition
    stroke       = r("stroke"),
    diabetes     = r("diabe")
  )
}

# Time-invariant, so no wave in the name.
hrs_fixed_names <- c(sex = "ragender", race = "raracem", education = "raedyrs")

# Systolic blood pressure and APOE do not come from the RAND Longitudinal
# File, so their column names depend on which product you downloaded and on
# the wave. Look them up in the codebook of that product and put them here:
# the physical measures section of the core interview numbers its blood
# pressure items differently in every wave, and the genetic file is
# distributed in more than one layout.
hrs_external_names <- c(
  systolic = "SYSTOLIC_BP",   # first systolic reading at the baseline wave
  apoe4    = "APOE4_ALLELES"  # count of epsilon-4 alleles, 0, 1 or 2
)

# Section 6 takes its confounders from Bennett et al. (2025).
hrs_confounders <- c("systolic", "age", "sex", "education", "race", "bmi",
                     "apoe4", "smoking", "depression", "heart", "stroke",
                     "diabetes")


# ---------------------------------------------------------------------------
# Cohort derivation
# ---------------------------------------------------------------------------

#' Whether a respondent is on an antihypertensive at one wave
#'
#' RAND HRS asks about blood pressure medication only of respondents who report
#' ever having been told they have high blood pressure, so r<wave>rxhibp is NA
#' under that skip pattern rather than 0. A respondent who has never been
#' diagnosed is by construction not taking one, and reading their NA as missing
#' would quietly delete the whole never-hypertensive stratum. A respondent who
#' was asked and did not answer is genuinely unknown and stays NA.
#'
#' This is also the exposure of section 6.1, "using any antihypertensive
#' treatment and having hypertension". The second half matters: someone taking
#' one of these drugs for heart failure or migraine is not the patient the
#' question is about, and counting them as treated would dilute the effect.
#' In this variable the two conditions coincide, because rxhibp records
#' medication taken for high blood pressure specifically.
#'
#' @param hrs Data frame with the RAND HRS columns.
#' @param names Output of \code{hrs_rand_names()} for the wave wanted.
#'
#' @return A logical vector, NA only where the answer is missing rather than
#'   skipped.
on_antihypertensive <- function(hrs, names) {
  hypertensive <- hrs[[names["hypertension"]]]
  medication   <- hrs[[names["bp_medication"]]]
  ifelse(!is.na(hypertensive) & hypertensive == 0, FALSE, medication == 1)
}


#' Build the analysis cohort from a merged HRS data frame
#'
#' Section 6 describes a new-user design over three waves: baseline
#' characteristics and the exclusions come from the baseline survey, treatment
#' is ascertained at the next survey two years later, and the outcome is read
#' four years after that. Excluding anyone already on an antihypertensive at
#' baseline is what makes the treated group new users rather than prevalent
#' ones, so that the comparison is not distorted by patients who had already
#' tolerated the drug for years.
#'
#' The paper says only that the outcome is taken from "the survey conducted
#' four years later" without naming the anchor. This function reads that as
#' four years after treatment ascertainment, so 2010, 2012 and 2016, which
#' spans the 2010 to 2016 range the paper mentions. Set
#' \code{outcome_year = 2014} for the other reading, four years after
#' baseline. That one reads waves 10, 11 and 12 and never touches wave 13, so
#' the 2016 survey plays no part in it.
#'
#' @param hrs Data frame with one row per respondent, holding the RAND HRS
#'   columns for every wave used together with the blood pressure and APOE
#'   columns named in \code{hrs_external_names}.
#' @param baseline_year Survey the exclusions and confounders are read from.
#' @param treatment_year Survey treatment is ascertained at.
#' @param outcome_year Survey the outcome is read from.
#'
#' @return A data frame with the treatment \code{antihypertensive}, the
#'   outcome \code{dementia} and the confounders of \code{hrs_confounders},
#'   and nothing else.
hrs_prepare <- function(hrs, baseline_year = 2010, treatment_year = 2012,
                        outcome_year = 2016) {

  base <- hrs_rand_names(baseline_year)
  trt  <- hrs_rand_names(treatment_year)
  out  <- hrs_rand_names(outcome_year)

  required <- c(base, trt[c("interviewed", "hypertension", "bp_medication")],
                out[c("interviewed", "dementia", "alzheimers", "memory")],
                hrs_fixed_names, hrs_external_names)
  absent <- setdiff(required, names(hrs))
  if (length(absent) > 0) {
    stop("These columns are not in the merged data frame: ", toString(absent),
         ". Check the wave numbering against your codebook, and set ",
         "hrs_external_names to the blood pressure and APOE column names ",
         "that your download actually uses.", call. = FALSE)
  }

  # Cognitive impairment at any of the three waves is recorded as ever having
  # been told, so a respondent who reports it at baseline is prevalent rather
  # than incident and cannot contribute an incident outcome.
  impaired <- function(names) {
    hrs[[names["dementia"]]] == 1 | hrs[[names["alzheimers"]]] == 1 |
      hrs[[names["memory"]]] == 1
  }

  # Each step is reported, because the size of the analysed cohort relative to
  # the study is the first thing a reader of an HRS analysis asks about.
  keep <- rep(TRUE, nrow(hrs))
  drop_to <- function(keep, condition, why) {
    kept <- keep & !is.na(condition) & condition
    message("  ", format(sum(kept), width = 6), "  ", why)
    kept
  }

  message("Cohort derivation (section 6.1):")
  message("  ", format(sum(keep), width = 6), "  respondents in the file")
  keep <- drop_to(keep, hrs[[base["interviewed"]]] == 1,
                  paste("interviewed at", baseline_year))
  keep <- drop_to(keep, hrs[[base["age"]]] >= 65, "aged 65 or over")
  keep <- drop_to(keep, hrs[[hrs_external_names["systolic"]]] > 130,
                  "systolic blood pressure above 130")
  keep <- drop_to(keep, !impaired(base),
                  "no dementia, Alzheimer's or memory problem at baseline")
  keep <- drop_to(keep, !on_antihypertensive(hrs, base),
                  "not already on an antihypertensive")
  keep <- drop_to(keep, hrs[[trt["interviewed"]]] == 1,
                  paste("interviewed at", treatment_year))
  keep <- drop_to(keep, hrs[[out["interviewed"]]] == 1,
                  paste("interviewed at", outcome_year))

  cohort <- hrs[keep, , drop = FALSE]

  # The exposure of section 6.1, read at the treatment wave. Everyone here was
  # free of antihypertensives at baseline, so this is a new-user indicator.
  antihypertensive <- as.integer(on_antihypertensive(cohort, trt))

  prepared <- data.frame(
    antihypertensive = antihypertensive,
    dementia         = as.integer(impaired(out)[keep]),
    systolic         = cohort[[hrs_external_names["systolic"]]],
    age              = cohort[[base["age"]]],
    sex              = factor(cohort[[hrs_fixed_names["sex"]]],
                              levels = 1:2, labels = c("male", "female")),
    education        = cohort[[hrs_fixed_names["education"]]],
    race             = factor(cohort[[hrs_fixed_names["race"]]],
                              levels = 1:3,
                              labels = c("white", "black", "other")),
    bmi              = cohort[[base["bmi"]]],
    apoe4            = cohort[[hrs_external_names["apoe4"]]],
    smoking          = as.integer(cohort[[base["smoking"]]] == 1),
    depression       = cohort[[base["depression"]]],
    heart            = as.integer(cohort[[base["heart"]]] == 1),
    stroke           = as.integer(cohort[[base["stroke"]]] == 1),
    diabetes         = as.integer(cohort[[base["diabetes"]]] == 1)
  )

  incomplete <- !stats::complete.cases(prepared)
  if (any(incomplete)) {
    message("  ", format(sum(!incomplete), width = 6),
            "  complete on every analysis variable")
    prepared <- prepared[!incomplete, , drop = FALSE]
  }

  rownames(prepared) <- NULL
  prepared
}


#' Formulas for the section 6 analysis
#'
#' @return A list with the outcome and propensity score formulas, ready for
#'   \code{DRBayes::drbayes_pc()}.
hrs_formulas <- function() {
  covariates <- paste(hrs_confounders, collapse = " + ")
  list(
    outcome = stats::as.formula(
      paste("dementia ~ antihypertensive +", covariates)),
    ps = stats::as.formula(paste("antihypertensive ~", covariates))
  )
}


# ---------------------------------------------------------------------------
# Run, if the files are there
# ---------------------------------------------------------------------------

hrs_dir <- Sys.getenv("HRS_DATA_DIR", unset = "")

if (!nzchar(hrs_dir)) {
  writeLines(c(
    "HRS_DATA_DIR is not set, so there is nothing to prepare.",
    "",
    "Obtain the three products listed at the top of this script from",
    "https://hrsdata.isr.umich.edu/data-products/, merge them into one data",
    "frame with one row per respondent, save it as hrs_merged.rds in a",
    "directory of your choice, and rerun with HRS_DATA_DIR set to that",
    "directory. hrs_prepare() then derives the exposure, the outcome and the",
    "confounders.",
    "",
    paste("Confounders (section 6.1):", toString(hrs_confounders)),
    "Formulas:"))
  for (formula in hrs_formulas()) {
    writeLines(paste("  ", paste(trimws(deparse(formula)), collapse = " ")))
  }
} else {
  merged_file <- file.path(hrs_dir, "hrs_merged.rds")
  if (!file.exists(merged_file)) {
    stop("HRS_DATA_DIR is ", hrs_dir, " but there is no hrs_merged.rds in ",
         "it. Merge the RAND HRS Longitudinal File, the physical measures ",
         "and the APOE file into one data frame with one row per respondent ",
         "and save it there.", call. = FALSE)
  }

  hrs <- hrs_prepare(readRDS(merged_file))

  output <- file.path(hrs_dir, "hrs_cohort.rds")
  saveRDS(hrs, output)
  message("Wrote ", output, " with ", nrow(hrs), " respondents, ",
          sum(hrs$antihypertensive), " treated.")
  message("The HRS conditions of use do not allow this file to be ",
          "redistributed, so keep it out of the repository.")
}
