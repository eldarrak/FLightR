test_that("tsoutliers.test returns a numeric score for each observation", {
  series <- c(1, 1.1, 0.9, 1.05, 8, 1, 0.95, 1.1, 1)

  score <- tsoutliers.test(series, plot = FALSE)
  expect_type(score, "double")
  expect_length(score, length(series))
  expect_true(all(score >= 0))
  expect_true(score[5] > 0)
})

test_that("tsoutliers.test gives low scores for a smooth series", {
  smooth <- sin(seq(0, 2 * pi, length.out = 30))

  score <- tsoutliers.test(smooth, plot = FALSE)
  expect_length(score, length(smooth))
  expect_true(all(score >= 0))
  expect_lt(max(score), 1)
})

test_that("tsoutliers.test plot path returns invisibly", {
  series <- c(rep(1, 10), 5, rep(1, 10))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- tsoutliers.test(series, plot = TRUE)

  expect_type(plotted, "double")
  expect_length(plotted, length(series))
})
