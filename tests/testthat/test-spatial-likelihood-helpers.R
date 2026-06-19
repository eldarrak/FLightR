fake_likelihood_calibration <- function() {
  list(
    Parameters = list(
      calibration.type = "parametric.slope",
      LogSlope = c(0.5, 1),
      log.light.borders = c(-Inf, Inf),
      log.irrad.borders = c(-Inf, Inf),
      impute.on.boundaries = FALSE,
      mean.of.individual.slope.sigma = 0.1
    ),
    time_correction_fun = function(declination, dusk) 0.5,
    lat_correction_fun = function(lat) 0,
    c_fun = NULL
  )
}

fake_twilight_matrices <- function(n_dusk = 2, n_dawn = 2) {
  base <- as.POSIXct("2020-03-20 04:00:00", tz = "GMT")
  time_col <- as.numeric(base + seq_len(49) * 300)
  list(
    dusk_time = matrix(rep(time_col, n_dusk), nrow = 49),
    dawn_time = matrix(rep(time_col, n_dawn), nrow = 49),
    dusk_light = matrix(rep(seq(2, 4, length.out = 49), n_dusk), nrow = 49),
    dawn_light = matrix(rep(seq(4, 2, length.out = 49), n_dawn), nrow = 49)
  )
}

test_that("get.probs.lm returns finite probability and optional slope details", {
  model <- list(coefficients = c(0, 0.5), stderr = c(0, 0))

  probability <- get.probs.lm(
    model,
    calibration = NULL,
    time_correction = 0.5,
    Calib.param = c(0.5, 1),
    delta = 0,
    likelihood.correction = FALSE
  )
  expect_type(probability, "double")
  expect_length(probability, 1)
  expect_true(is.finite(probability))
  expect_gt(probability, 0)

  with_slopes <- get.probs.lm(
    model,
    calibration = NULL,
    time_correction = 0.5,
    Calib.param = c(0.5, 1),
    delta = 0,
    return.slopes = TRUE,
    likelihood.correction = FALSE
  )
  expect_equal(length(with_slopes), 3)
  expect_equal(with_slopes[2], 0.5)
})

test_that("get.probs.lm falls back for non-finite slope stderr", {
  model <- list(coefficients = c(0, 0.5), stderr = c(0, Inf))
  calibration <- fake_likelihood_calibration()

  set.seed(1)
  probability <- get.probs.lm(
    model,
    calibration = calibration,
    Twilight.time.vector = as.numeric(as.POSIXct("2020-03-20 00:00:00", tz = "GMT") + seq_len(49) * 300),
    Calib.param = c(0.5, 1),
    x = c(0, 0),
    likelihood.correction = FALSE
  )
  expect_type(probability, "double")
  expect_true(is.finite(probability))
})

test_that("get.prob.surface returns a probability per grid cell", {
  mats <- fake_twilight_matrices(n_dusk = 1, n_dawn = 1)
  calibration <- fake_likelihood_calibration()
  grid <- matrix(c(0, 0, 1, 1, 1, 1), ncol = 3, byrow = TRUE)
  colnames(grid) <- c("lon", "lat", "mask")

  set.seed(1)
  probs <- suppressMessages(get.prob.surface(
    Twilight.ID = 1,
    dusk = TRUE,
    Twilight.time.mat = mats$dusk_time,
    Twilight.log.light.mat = mats$dusk_light,
    Calib.param = calibration$Parameters$LogSlope,
    log.irrad.borders = calibration$Parameters$log.irrad.borders,
    Grid = grid,
    log.light.borders = calibration$Parameters$log.light.borders,
    calibration = calibration,
    impute.on.boundaries = FALSE,
    likelihood.correction = FALSE
  ))

  expect_type(probs, "double")
  expect_length(probs, nrow(grid))
  expect_true(all(is.finite(probs)))
  expect_true(all(probs >= 0))
})

test_that("nonparametric and current slope probability helpers return finite values", {
  calibration <- fake_likelihood_calibration()
  times <- as.numeric(as.POSIXct("2020-03-20 04:00:00", tz = "GMT") + seq_len(49) * 300)
  solar <- solar.FLightR(as.POSIXct(times, origin = "1970-01-01", tz = "GMT"))
  log_light <- seq(2, 4, length.out = 49)

  nonparam <- get.probs.nonparam.slope(
    Slopes = c(0.4, 0.5, 0.6),
    calibration = calibration,
    time_correction = 0.5,
    Calib.param = c(0.5, 1),
    delta = 0,
    return.slopes = TRUE
  )
  expect_equal(length(nonparam), 3)
  expect_true(is.finite(nonparam[1]))

  boundary_data <- check.boundaries(
    x = c(0, 0),
    Twilight.solar.vector = solar,
    Twilight.log.light.vector = log_light,
    log.light.borders = c(-Inf, Inf),
    log.irrad.borders = c(-Inf, Inf),
    dusk = TRUE,
    Twilight.time.vector = times
  )
  expect_true(is.matrix(boundary_data))
  expect_true(ncol(boundary_data) >= 2)

  set.seed(1)
  current_parametric <- get.current.slope.prob(
    x = c(0, 0),
    calibration = calibration,
    Twilight.solar.vector = solar,
    Twilight.log.light.vector = log_light,
    log.light.borders = c(-Inf, Inf),
    log.irrad.borders = c(-Inf, Inf),
    dusk = TRUE,
    Twilight.time.vector = times,
    time_correction = 0.5,
    Calib.param = calibration$Parameters$LogSlope,
    delta = 0,
    likelihood.correction = FALSE
  )
  expect_type(current_parametric, "double")
  expect_true(is.finite(current_parametric))

  calibration$Parameters$calibration.type <- "nonparametric.slope"
  current_nonparametric <- get.current.slope.prob(
    x = c(0, 0),
    calibration = calibration,
    Twilight.solar.vector = solar,
    Twilight.log.light.vector = log_light,
    log.light.borders = c(-Inf, Inf),
    log.irrad.borders = c(-Inf, Inf),
    dusk = TRUE,
    Twilight.time.vector = times,
    time_correction = 0.5,
    Calib.param = calibration$Parameters$LogSlope,
    delta = 0,
    likelihood.correction = FALSE
  )
  expect_type(current_nonparametric, "double")
  expect_true(is.finite(current_nonparametric))
})

test_that("get.Phys.Mat.parallel orders dusk and dawn columns and floors tiny probabilities", {
  mats <- fake_twilight_matrices(n_dusk = 2, n_dawn = 2)
  calibration <- fake_likelihood_calibration()
  all.out <- list(
    Spatial = list(Grid = matrix(0, nrow = 2, ncol = 3)),
    Indices = list(Matrix.Index.Table = data.frame(Dusk = c(TRUE, FALSE, TRUE, FALSE)))
  )

  testthat::local_mocked_bindings(
    get.prob.surface = function(Twilight.ID, dusk, ...) {
      if (dusk) {
        c(Twilight.ID, 0)
      } else {
        c(Twilight.ID + 100, 1e-100)
      }
    }
  )

  phys <- suppressMessages(get.Phys.Mat.parallel(
    all.out = all.out,
    Twilight.time.mat.dusk = mats$dusk_time,
    Twilight.log.light.mat.dusk = mats$dusk_light,
    Twilight.time.mat.dawn = mats$dawn_time,
    Twilight.log.light.mat.dawn = mats$dawn_light,
    threads = 1,
    calibration = calibration,
    likelihood.correction = FALSE
  ))

  expect_equal(dim(phys), c(2, 4))
  expect_equal(phys[1, ], c(1, 101, 2, 102))
  expect_equal(phys[2, ], rep(1e-70, 4))
})
