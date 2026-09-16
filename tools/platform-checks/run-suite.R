cat(sprintf("  %s   long.double %s   BLAS %s\n", R.version$platform,
            capabilities("long.double"), basename(extSoftVersion()[["BLAS"]])))
suppressPackageStartupMessages(library(testthat))
r <- test_dir("tests/testthat", package = "DRBayes", load_package = "installed",
              reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(r)
cat(sprintf("  passed %d  failed %d  errors %d  skipped %d\n",
            sum(d$passed), sum(d$failed), sum(d$error), sum(d$skipped)))
bad <- d[d$failed > 0 | d$error, c("file", "test")]
if (nrow(bad)) print(bad, row.names = FALSE) else cat("  no failures\n")
