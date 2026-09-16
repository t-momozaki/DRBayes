# These pin the build rather than the arithmetic. A package that failed to
# compile, or that lost its useDynLib entry when the documentation was
# regenerated, would otherwise load and run correctly on everything except the
# compiled paths, and the failure would surface as a wrong answer somewhere far
# away rather than here.

test_that("the compiled routines are registered and reachable", {
  info <- getLoadedDLLs()[["DRBayes"]]
  skip_if(is.null(info), "the source tree is loaded, not the installed package")

  registered <- getDLLRegisteredRoutines(info)$.Call
  expect_true(length(registered) > 0)

  # Rcpp prefixes each entry point with the package name, so the compiled
  # routines are the ones whose registered name carries it.
  expect_true(any(grepl("^_DRBayes_", names(registered))))
})


test_that("Armadillo reaches the same BLAS R does", {
  skip_if_not(exists("drb_blas_check", envir = asNamespace("DRBayes")),
              "no compiled kernel in this build")

  set.seed(1)
  a <- matrix(stats::rnorm(60), 6, 10)
  b <- matrix(stats::rnorm(70), 10, 7)

  # A build that linked against no BLAS, or against one that disagrees with R's,
  # is the failure this catches: the kernels are written on the assumption that
  # their matrix products are the products R would have made.
  expect_equal(DRBayes:::drb_blas_check(a, b), sum(a %*% b))
})


test_that("the kernels are told the width R actually sums in", {
  skip_if_not(exists("long_double_sums", envir = asNamespace("DRBayes")),
              "no compiled kernel in this build")

  # Not a constant. R's own answer has to be forwarded, because the macro that
  # records it is private to R's build and no #ifdef in the package can see it.
  # A kernel wired to either answer would agree with R on about half the
  # platforms CRAN checks and disagree silently on the rest.
  expect_identical(DRBayes:::long_double_sums(),
                   isTRUE(capabilities("long.double")))

  wide <- DRBayes:::long_double_sums()
  x <- c(1e16, 1, 1, 1, -1e16)
  expect_identical(DRBayes:::drb_mean_one_pass(x, wide),
                   .colMeans(x, length(x), 1L))

  if (wide) {
    # Where long double is wider than double the two instantiations must part
    # company on a sum that cancels, or the argument is being ignored. Where it
    # is not wider, as on aarch64 macOS, they are the same type and agreeing is
    # the correct answer, so there is nothing to assert.
    expect_false(identical(DRBayes:::drb_mean_one_pass(x, TRUE),
                           DRBayes:::drb_mean_one_pass(x, FALSE)))
  }
})


test_that("the Stan backend still finds its models now that src exists", {
  skip_on_cran()
  skip_if_not_installed("instantiate")
  skip_if_not(nzchar(system.file(package = "DRBayes")),
              "installed package only")

  # instantiate would normally put the models under src/stan and compile them at
  # install time; this package ships them in inst/stan and compiles on first use
  # instead, so a src/ directory of our own has nothing to collide with. The
  # test exists because that separation is an assumption, not a guarantee.
  for (model in DRBayes:::stan_model_names()) {
    file <- system.file("stan", paste0("bayes_", model, ".stan"),
                        package = "DRBayes")
    expect_true(nzchar(file), info = model)
  }
})
