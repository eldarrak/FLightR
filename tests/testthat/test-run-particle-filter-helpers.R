test_that("pf message helpers respect quiet, normal, and debug modes", {
  expect_silent(pf_message("quiet", verbose = FALSE))
  expect_message(pf_message("normal", verbose = TRUE), "normal")
  expect_message(pf_message("debug-normal", verbose = "debug"), "debug-normal")

  expect_silent(pf_debug("quiet-debug", verbose = TRUE))
  expect_message(pf_debug("debug", verbose = "debug"), "debug")
})

test_that("small particle-filter math helpers return expected structures", {
  expect_equal(mu.sigma.truncnorm(rep(10, 4), a = 0, b = 20), c(10, 0))

  vm <- flightr.dvonmises(list(0, 0), mykap = 1)
  expect_type(vm, "double")
  expect_gt(vm, 0)

  grid <- matrix(c(
    0, 0, 1,
    1, 0, 1,
    0, 1, 1
  ), ncol = 3, byrow = TRUE)
  colnames(grid) <- c("lon", "lat", "mask")
  in_data <- list(Spatial = list(Grid = grid))

  expect_equal(round(dir_fun(c(1, 2), in_data)), 90)

  transition <- get.transition.rle(c(1, 1, 2), c(1, 2, 3), n_grid = 3)
  result <- list(Spatial = list(Grid = grid))
  distances <- dist.fun(transition$values, result, transition)
  expect_length(distances, length(transition$values))
  expect_true(all(distances >= 0))
})

test_that("movement candidate builders produce bounded candidate lists", {
  grid <- matrix(c(
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    10, 10, 1
  ), ncol = 3, byrow = TRUE)
  colnames(grid) <- c("lon", "lat", "mask")

  candidates <- build.grid.movement.candidates(grid, a = 0, b = 200)
  expect_equal(candidates$grid_n, nrow(grid))
  expect_length(candidates$to, nrow(grid))
  expect_true(all(candidates$distance[[1]] <= 200))
  expect_true(1 %in% candidates$to[[1]])

  cache <- build.partial.movement.cache(grid, a = 0, b = 200, max_cache_mb = 1)
  entry1 <- get.partial.movement.entry(cache, 1)
  entry2 <- get.partial.movement.entry(cache, 1)
  expect_equal(entry1$to, entry2$to)
  expect_equal(cache$misses, 1L)
  expect_equal(cache$hits, 1L)

  diagnostics <- partial.movement.cache.diagnostics(cache)
  expect_s3_class(diagnostics, "data.frame")
  expect_equal(diagnostics$cached_origins, 1)
  expect_equal(diagnostics$cache_hits, 1)
})

test_that("get.LL.PF handles matrix and RLE particle histories", {
  phys <- matrix(c(
    0.5, 0.25, 0.2,
    0.5, 0.75, 0.8
  ), nrow = 2, byrow = TRUE)
  in_data <- list(Spatial = list(Phys.Mat = phys))

  matrix_history <- matrix(c(1, 1, 1, 2, 2, 2), nrow = 2, byrow = TRUE)
  ll_matrix <- get.LL.PF(in_data, matrix_history)
  expect_type(ll_matrix, "double")
  expect_true(is.finite(ll_matrix))

  rle_history <- list(
    list(values = c(1, 2), lengths = c(1, 1)),
    list(values = c(1, 2), lengths = c(1, 1)),
    list(values = c(1, 2), lengths = c(1, 1))
  )
  ll_rle <- get.LL.PF(in_data, rle_history)
  expect_type(ll_rle, "double")
  expect_true(is.finite(ll_rle))
})
