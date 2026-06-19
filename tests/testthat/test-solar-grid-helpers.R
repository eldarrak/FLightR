test_that("solar helper functions return finite, plausible values", {
  times <- as.POSIXct(c("2020-03-20 12:00:00", "2020-06-21 12:00:00"),
    tz = "GMT")

  trip <- solar.FLightR(times, mode = "tripEstimation")
  expect_named(trip, c("solarTime", "sinSolarDec", "cosSolarDec"))
  expect_length(trip$solarTime, 2)
  expect_true(all(is.finite(trip$solarTime)))
  expect_true(all(abs(trip$sinSolarDec) <= 1))
  expect_true(all(abs(trip$cosSolarDec) <= 1))

  geo <- solar.FLightR(times, mode = "GeoLight")
  expect_named(geo, c("theta.G", "alpha", "sinSolarDec", "cosSolarDec"))
  expect_length(geo$theta.G, 2)
  expect_true(all(is.finite(geo$theta.G)))

  equinox_sun <- list(
    solarTime = trip$solarTime[1],
    sinSolarDec = trip$sinSolarDec[1],
    cosSolarDec = trip$cosSolarDec[1]
  )
  equinox_elevation <- elevation.FLightR(0, 0, equinox_sun, mode = "tripEstimation")
  expect_gt(equinox_elevation, 80)

  declination <- get.declination(times)
  expect_length(declination, 2)
  expect_true(all(is.finite(declination)))
  expect_gt(declination[2], 20)
})

test_that("declination accepts numeric epoch seconds", {
  times <- as.POSIXct(c("2020-01-01 00:00:00", "2020-07-01 00:00:00"),
    tz = "GMT")
  expect_equal(
    get.declination(as.numeric(times)),
    get.declination(times),
    tolerance = 1e-10
  )
})

test_that("make.grid records bounds and supports dateline-crossing extents", {
  grid <- make.grid(
    left = 170, right = -170, bottom = -10, top = 10,
    distance.from.land.allowed.to.use = c(-Inf, Inf),
    distance.from.land.allowed.to.stay = c(-Inf, Inf),
    plot = FALSE
  )

  expect_true(nrow(grid) > 0)
  expect_equal(attr(grid, "left"), 170)
  expect_equal(attr(grid, "right"), -170)
  expect_true(all(grid[, 1] >= 170 | grid[, 1] <= -170))
  expect_true(all(grid[, 2] >= -10 & grid[, 2] <= 10))
})

test_that("make.grid can return distance-to-land column", {
  grid <- make.grid(
    left = 0, right = 5, bottom = 50, top = 52,
    distance.from.land.allowed.to.use = c(-Inf, Inf),
    distance.from.land.allowed.to.stay = c(-Inf, Inf),
    return.distances = TRUE,
    plot = FALSE
  )

  expect_equal(ncol(grid), 4)
  expect_true(all(is.finite(grid[, 4])))
})

test_that("make.grid applies probability of staying and adjusts stay bounds", {
  grid <- make.grid(
    left = 0, right = 10, bottom = 50, top = 56,
    distance.from.land.allowed.to.use = c(-100, Inf),
    distance.from.land.allowed.to.stay = c(0, 10),
    probability.of.staying = 0.25,
    plot = FALSE
  )

  expect_true(any(grid[, 3] == 0.25))
  expect_true(any(grid[, 3] == 1))
  expect_equal(attr(grid, "distance.from.land.allowed.to.stay"), c(-25, 25))
  expect_equal(attr(grid, "distance.from.land.allowed.to.use"), c(-100, Inf))
})
