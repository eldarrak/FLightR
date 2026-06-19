fake_schedule_result <- function() {
  times <- as.POSIXct("2020-01-01 00:00:00", tz = "GMT") + seq(0, 4) * 3600
  list(
    Indices = list(Matrix.Index.Table = data.frame(time = times)),
    Spatial = list(Grid = matrix(c(0, 0, 1, 1, 0, 1), ncol = 3, byrow = TRUE)),
    Results = list(
      Points.rle = list(
        list(values = 2, lengths = 4),
        list(values = c(1, 2), lengths = c(1, 3)),
        list(values = c(1, 2), lengths = c(3, 1)),
        list(values = 1, lengths = 4),
        list(values = 2, lengths = 4)
      ),
      Quantiles = data.frame(
        time = times[1:2],
        Medianlon = c(5, 6),
        Medianlat = c(52, 53),
        TrdQu.lon = c(5.5, 6.5),
        FstQu.lon = c(4.5, 5.5),
        TrdQu.lat = c(52.5, 53.5),
        FstQu.lat = c(51.5, 52.5),
        UCI.lon = c(6, 7),
        LCI.lon = c(4, 5),
        UCI.lat = c(53, 54),
        LCI.lat = c(51, 52)
      )
    )
  )
}

test_that("get.prob.of.being.in calculates marginal probabilities from RLE points", {
  result <- fake_schedule_result()

  expect_equal(
    get.prob.of.being.in(result, 1),
    c(0, 0.25, 0.75, 1, 0)
  )
  expect_equal(
    get.prob.of.being.in(result, 2),
    c(1, 0.75, 0.25, 0, 1)
  )
})

test_that("find.time interpolates threshold crossings and returns NA when absent", {
  times <- as.POSIXct("2020-01-01 00:00:00", tz = "GMT") + seq(0, 4) * 3600
  probabilities <- c(0, 0.25, 0.75, 1, 0)

  crossings <- find.time(probabilities, times, prob.cutoff = 0.5)
  expect_s3_class(crossings, "POSIXct")
  expect_equal(length(crossings), 2)
  expect_equal(
    as.numeric(crossings[1] - times[2], units = "mins"),
    30,
    tolerance = 1e-8
  )

  expect_true(is.na(find.time(rep(0.1, 5), times, prob.cutoff = 0.5)))
})

test_that("find.times.distribution returns crossing quantiles", {
  result <- fake_schedule_result()
  schedule <- find.times.distribution(result, 1)

  expect_s3_class(schedule, "data.frame")
  expect_named(schedule, c("Q.025", "Q.25", "Q.50", "Q.75", "Q.975"))
  expect_s3_class(schedule$Q.50, "POSIXct")
  expect_true(all(!is.na(schedule$Q.50)))
})

test_that("point distribution helpers summarize selected twilights", {
  result <- fake_schedule_result()

  distribution <- get_points_distribution(result, c(2, 3))
  expect_equal(distribution, c(4, 4))

  most_likely <- get_utilisation_points(result, c(2, 3), percentile = 0.4)
  expect_equal(most_likely, 1:2)

  top_point <- get_utilisation_points(result, c(4), percentile = 0.01)
  expect_equal(top_point, 1)
})

test_that("FLightR2Movebank formats median and 95 percent intervals", {
  result <- fake_schedule_result()

  movebank_50 <- FLightR2Movebank(result, alpha = 0.5)
  expect_s3_class(movebank_50, "data.frame")
  expect_equal(nrow(movebank_50), 2)
  expect_match(movebank_50$timestamp[1], "2020-01-01T00:00:00.000Z")
  expect_equal(movebank_50$location.long, c(5, 6))
  expect_equal(movebank_50$long.upper, c(5.5, 6.5))

  movebank_95 <- FLightR2Movebank(result, alpha = 0.95)
  expect_equal(movebank_95$long.upper, c(6, 7))
  expect_equal(movebank_95$lat.lower, c(51, 52))

  expect_error(FLightR2Movebank(result, alpha = 0.75), "0.5 or 0.95")
})
