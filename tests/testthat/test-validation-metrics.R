source(testthat::test_path("../../R/validation_metrics.R"))

fake_result <- function() {
  list(
    twilight = as.POSIXct(
      c("2020-01-01 06:00:00", "2020-01-01 18:00:00", "2020-01-02 06:00:00"),
      tz = "UTC"
    ),
    posterior_points = list(
      data.frame(lon = c(4.9, 5.0, 5.1), lat = c(51.9, 52.0, 52.1), weight = c(0.2, 0.7, 0.1)),
      data.frame(lon = c(5.4, 5.5, 5.6), lat = c(52.4, 52.5, 52.6), weight = c(0.2, 0.6, 0.2)),
      data.frame(lon = c(5.9, 6.0, 6.1), lat = c(52.9, 53.0, 53.1), weight = c(0.2, 0.5, 0.3))
    )
  )
}

fake_gps <- function() {
  data.frame(
    time = as.POSIXct(
      c("2020-01-01 06:15:00", "2020-01-01 17:45:00", "2020-01-03 06:00:00"),
      tz = "UTC"
    ),
    lon = c(5.0, 5.5, 9.0),
    lat = c(52.0, 52.5, 55.0)
  )
}

test_that("GPS fixes are matched to nearby twilight times", {
  matched <- match_gps_to_twilights(fake_result(), fake_gps(), max_time_diff_hours = 1)

  expect_equal(nrow(matched), 2)
  expect_named(
    matched,
    c("twilight_index", "twilight_time", "gps_time", "time_diff_hours", "gps_lon", "gps_lat")
  )
  expect_equal(matched$twilight_index, c(1, 2))
})

test_that("posterior point distribution is normalized", {
  points <- posterior_point_distribution(fake_result(), 1)

  expect_named(points, c("lon", "lat", "weight"))
  expect_equal(sum(points$weight), 1)
})

test_that("validation summary returns expected columns", {
  result <- fake_result()
  matched <- match_gps_to_twilights(result, fake_gps(), max_time_diff_hours = 1)
  summary <- validation_summary(result, matched)

  expect_equal(nrow(summary), 1)
  expect_named(
    summary,
    c(
      "n_matched",
      "median_error_km",
      "mean_error_km",
      "mode_error_km",
      "ci95_coverage",
      "notes"
    )
  )
  expect_equal(summary$n_matched, 2)
})
