fake_summary_result <- function() {
  times <- as.POSIXct("2020-01-01 00:00:00", tz = "GMT") + seq(0, 4) * 3600
  grid <- matrix(c(
    0, 0, 1,
    1, 0, 1,
    0, 1, 1
  ), ncol = 3, byrow = TRUE)
  colnames(grid) <- c("lon", "lat", "mask")

  list(
    Indices = list(Matrix.Index.Table = data.frame(time = times)),
    Spatial = list(Grid = grid),
    Results = list(
      Points.rle = list(
        list(values = c(1, 2), lengths = c(3, 1)),
        list(values = c(1, 2), lengths = c(2, 2)),
        list(values = c(2, 3), lengths = c(3, 1)),
        list(values = c(3, 2), lengths = c(3, 1)),
        list(values = c(3, 1), lengths = c(2, 2))
      ),
      Transitions.rle = list(
        get.transition.rle(c(1, 1, 2, 2), c(1, 2, 2, 3), n_grid = 3),
        get.transition.rle(c(2, 2, 3, 3), c(2, 3, 3, 1), n_grid = 3)
      ),
      Movement.results = data.frame(
        time = times[1:2],
        Decision = c(0.1, 0.9)
      ),
      Quantiles = data.frame(
        time = times,
        Meanlon = c(0, 0.2, 0.8, 0.2, 0),
        Meanlat = c(0, 0.1, 0.1, 0.8, 1),
        Medianlon = c(0, 0, 1, 0, 0),
        Medianlat = c(0, 0, 0, 1, 1),
        FstQu.lon = c(0, 0, 0, 0, 0),
        TrdQu.lon = c(1, 1, 1, 1, 1),
        FstQu.lat = c(0, 0, 0, 0, 0),
        TrdQu.lat = c(1, 1, 1, 1, 1),
        LCI.lon = c(0, 0, 0, 0, 0),
        UCI.lon = c(1, 1, 1, 1, 1),
        LCI.lat = c(0, 0, 0, 0, 0),
        UCI.lat = c(1, 1, 1, 1, 1)
      )
    )
  )
}

test_that("get_ZI_distances summarizes transition RLE distances", {
  result <- fake_summary_result()
  distances <- get_ZI_distances(result)

  expect_s3_class(distances, "data.frame")
  expect_equal(nrow(distances), 2)
  expect_named(distances, c("Departure", "Arrival", "25%", "50%", "75%", "Mean"))
  expect_true(all(distances$Mean >= 0))
})

test_that("stationary stats summarize point distributions", {
  set.seed(1)
  result <- fake_summary_result()
  stats <- get_stationary_stats(result, 1:2)

  expect_s3_class(stats, "data.frame")
  expect_true(all(c("Medianlat", "Meanlon", "LCI.lat", "UCI.lon") %in% names(stats)))
  expect_equal(nrow(stats), 1)
})

test_that("time boundary helper returns arrival and departure columns", {
  result <- fake_summary_result()
  bounds <- get_time_boundaries(result, 2:4, utilisation.distribution.percentile = 0.5)

  expect_s3_class(bounds, "data.frame")
  expect_true(any(grepl("^Arrival", names(bounds))))
  expect_true(any(grepl("^Departure", names(bounds))))
})

test_that("buffer helper and time-spent buffer return sf output", {
  result <- fake_summary_result()

  buffer <- get_buffer(c(0, 0), r = 1000)
  expect_s3_class(buffer, "sf")

  spent <- get_time_spent_buffer(result, dates = 1, percentile = 0.5, r = 1000)
  expect_type(spent, "list")
  expect_s3_class(spent$Buffer, "sf")
  expect_true(spent$nPoints >= 1)
})
