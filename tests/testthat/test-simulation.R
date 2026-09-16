# Reproducibility and global state. Every one of these fails against the
# pre-fix code, which set a seed inside each worker, switched RNGkind without
# restoring the seed, and reset the future plan to sequential.

test_that("a seed given to generate_dataset does not disturb the caller", {
  set.seed(7)
  reference <- rnorm(3)

  set.seed(7)
  invisible(generate_dataset(nn = 10, pp = 0, seed = 123))
  expect_equal(rnorm(3), reference)

  # The same seed still gives the same dataset.
  a <- generate_dataset(nn = 10, pp = 0, seed = 123)
  b <- generate_dataset(nn = 10, pp = 0, seed = 123)
  expect_identical(a, b)
})

test_that("generate_dataset implements the data generating process of the paper", {
  d <- generate_dataset(nn = 5000, pp = 0, seed = 1)

  expect_named(d, c("AA", "YY", "WW1", "WW2", "WW3", "WW4"))
  expect_true(all(d$AA %in% c(0, 1)))

  # No intercept in the propensity score, so treatment prevalence is one half.
  expect_lt(abs(mean(d$AA) - 0.5), 0.03)

  # The coefficients the paper specifies, recovered from the generated data.
  ps <- stats::glm(AA ~ WW1 + WW2 + WW3 + WW4, binomial, d)
  expect_equal(unname(coef(ps)), c(0, 1, -0.5, 0.25, 0.1), tolerance = 0.15)

  outcome <- stats::lm(YY ~ AA + WW1 + WW2 + WW3 + WW4, d)
  expect_equal(unname(coef(outcome)),
               c(100, 110, 27.4, 13.7, 13.7, 13.7), tolerance = 0.15)
})

test_that("the extra covariates are irrelevant, not confounders", {
  d <- generate_dataset(nn = 20000, pp = 5, seed = 2)
  extra <- d[, paste0("WW", 5:9)]

  # They enter neither model, so they are uncorrelated with both.
  expect_true(all(abs(cor(extra, d$AA)) < 0.05))
  expect_true(all(abs(cor(extra, d$YY)) < 0.05))

  # Each column keeps its own mean, which is what the byrow argument buys.
  expect_gt(diff(range(colMeans(extra))), 0.5)
})

test_that("a study is reproducible and replayable", {
  a <- run_parallel_simulation(4, nn = 20, pp = 0, seed = 2024, verbose = FALSE)
  b <- run_parallel_simulation(4, nn = 20, pp = 0, seed = 2024, verbose = FALSE)
  expect_identical(a, b)

  set.seed(11)
  first  <- run_parallel_simulation(4, nn = 20, pp = 0, verbose = FALSE)
  second <- run_parallel_simulation(4, nn = 20, pp = 0, verbose = FALSE)
  # Two successive calls must give different studies, otherwise a user
  # enlarging a simulation gets the same datasets back and reports a Monte
  # Carlo error near zero.
  expect_false(identical(first, second))

  # The seed actually used is recorded, so the run can be replayed.
  replay <- run_parallel_simulation(4, nn = 20, pp = 0,
                                    seed = attr(first, "seed"), verbose = FALSE)
  expect_identical(unclass(first), unclass(replay))
})

test_that("the datasets do not depend on how many workers run", {
  skip_on_cran()
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  sequential_run <- run_parallel_simulation(6, nn = 20, pp = 0, seed = 5,
                                            verbose = FALSE)
  parallel_run <- run_parallel_simulation(6, nn = 20, pp = 0, seed = 5,
                                          n_cores = 2, verbose = FALSE)
  expect_identical(sequential_run, parallel_run)
})

test_that("the caller's plan, RNG state and options survive the call", {
  skip_on_cran()

  old <- future::plan(future::multisession, workers = 2)
  on.exit({
    future::plan(old)
  }, add = TRUE)
  expect_equal(future::nbrOfWorkers(), 2)

  set.seed(5)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  option_set_before <- "future.rng.onMisuse" %in% names(options())

  invisible(run_parallel_simulation(3, nn = 20, pp = 0, seed = 1,
                                    n_cores = 2, verbose = FALSE))

  # The pre-fix code reset the plan to sequential, silently shutting down the
  # caller's workers for the rest of their session.
  expect_equal(future::nbrOfWorkers(), 2)
  # Switching RNGkind re-seeds, and the pre-fix on.exit restored only the kind.
  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)
  expect_identical("future.rng.onMisuse" %in% names(options()),
                   option_set_before)
})

test_that("a session that has never drawn gets its generator back", {
  # CRAN policy forbids a package changing the kind or the seed of the user's
  # generator. furrr switches to L'Ecuyer-CMRG to derive its substreams and
  # puts back whatever .Random.seed it found, so in a session that has never
  # drawn a random number there was nothing to put back and the switch stuck:
  # every set.seed() the caller made afterwards produced a different stream
  # than it would have done. The other tests here all run after something has
  # already seeded, which is exactly why they miss it.
  set.seed(7)
  reference <- rnorm(2)

  saved <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
  kind_before <- RNGkind()
  rm(".Random.seed", envir = globalenv())

  invisible(run_parallel_simulation(3, nn = 20, pp = 0, seed = 1,
                                    verbose = FALSE))

  expect_identical(RNGkind(), kind_before)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))

  # What the leak actually cost, seen from the caller's side.
  set.seed(7)
  expect_identical(rnorm(2), reference)
})

test_that("the call takes one draw without a seed and none with one", {
  set.seed(11)
  untouched <- .Random.seed

  invisible(run_parallel_simulation(3, nn = 20, pp = 0, seed = 99,
                                    verbose = FALSE))
  # A seeded study is a pure function of its arguments, so the caller's stream
  # has to come back exactly as it was.
  expect_identical(.Random.seed, untouched)

  set.seed(11)
  one_draw <- sample.int(.Machine$integer.max, 1L)
  after_one_draw <- .Random.seed

  set.seed(11)
  study <- run_parallel_simulation(3, nn = 20, pp = 0, verbose = FALSE)

  # Without a seed the documented effect on the caller's stream is that single
  # draw, the one that picks the seed. Leaving it anywhere else means a script
  # that seeds itself cannot say what it draws next.
  expect_identical(attr(study, "seed"), one_draw)
  expect_identical(.Random.seed, after_one_draw)
})

test_that("arguments are validated and progress can be silenced", {
  expect_error(run_parallel_simulation(0), "num_simulations must be")
  expect_error(run_parallel_simulation(2, nn = 0), "nn must be")
  expect_error(run_parallel_simulation(2, nn = 10, pp = 0, n_cores = 0),
               "n_cores must be")
  expect_error(run_parallel_simulation(2, nn = 10, pp = 0, seed = NA),
               "seed must be")
  expect_error(run_parallel_simulation(2, nn = 10, pp = 0, verbose = NA),
               "verbose must be")

  expect_silent(invisible(run_parallel_simulation(2, nn = 10, pp = 0, seed = 1,
                                                  verbose = FALSE)))
  expect_message(invisible(run_parallel_simulation(2, nn = 10, pp = 0, seed = 1,
                                                   verbose = TRUE)),
                 "Generating 2 datasets")
})

test_that("a rejected argument is named, and the right one", {
  # The validators of generate_dataset, which nothing else here reaches, and
  # the covariate count of run_parallel_simulation, which the test above
  # leaves out. What makes them worth a test is not that they fire but what
  # they say: a message naming the wrong argument sends a reader off to look
  # at data that are fine.
  expect_error(generate_dataset(nn = 0), "nn must be a single positive integer")
  expect_error(generate_dataset(nn = 10, pp = -1),
               "pp must be a single non-negative integer")
  expect_error(generate_dataset(nn = 10, pp = 2.5),
               "pp must be a single non-negative integer")
  expect_error(generate_dataset(nn = 10, pp = 0, seed = NA),
               "seed must be NULL or a single non-missing number")
  expect_error(run_parallel_simulation(2, nn = 10, pp = -1),
               "pp must be a single non-negative integer")
})

test_that("convert_to_array is reachable and returns a numeric array", {
  sims <- run_parallel_simulation(3, nn = 15, pp = 2, seed = 9, verbose = FALSE)
  arr <- convert_to_array(sims)

  expect_identical(dim(arr), c(15L, 8L, 3L))
  expect_identical(typeof(arr), "double")
  expect_identical(dimnames(arr)[[2]], names(sims[[1]]))
  expect_equal(arr[, , 1], as.matrix(sims[[1]]), ignore_attr = TRUE)

  expect_error(convert_to_array(list()), "non-empty list")
  expect_error(convert_to_array(list(1, 2)), "must be a data.frame")
  expect_error(convert_to_array(list(sims[[1]], sims[[1]][1:5, ])),
               "same dimensions")
})
