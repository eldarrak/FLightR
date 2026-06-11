context("transition encoding")

env<-new.env(parent=globalenv())
source_checkout_r("run_particle_filter.R", env=env)

test_that("transition encoding round-trips on small grids", {
  from<-c(1L, 12L, 99L)
  to<-c(2L, 25L, 100L)

  key<-env$encode_transition(from, to, n_grid=100L)
  decoded<-env$decode_transition(key, transition_base=env$transition_base(100L))

  expect_equal(unname(decoded[, "from"]), from)
  expect_equal(unname(decoded[, "to"]), to)
})

test_that("transition encoding round-trips above the legacy base", {
  from<-c(1L, 100000L, 120001L)
  to<-c(120001L, 100001L, 2L)

  key<-env$encode_transition(from, to, n_grid=120001L)
  decoded<-env$decode_transition(key, transition_base=env$transition_base(120001L))

  expect_equal(unname(decoded[, "from"]), from)
  expect_equal(unname(decoded[, "to"]), to)
})

test_that("transition RLE stores enough metadata for large-grid decoding", {
  from<-c(100000L, 2L, 100000L)
  to<-c(100001L, 3L, 100001L)

  transition_rle<-env$get.transition.rle(from, to, n_grid=150000L)
  decoded<-env$decode_transition_rle_values(transition_rle)
  expected_order<-order(env$encode_transition(from, to, n_grid=150000L))

  expect_equal(attr(transition_rle, "transition_base"), env$transition_base(150000L))
  expect_equal(attr(transition_rle, "n_grid"), 150000L)
  expect_equal(unname(decoded[, "from"]), from[expected_order][c(TRUE, diff(env$encode_transition(from, to, n_grid=150000L)[expected_order])!=0)])
  expect_equal(unname(decoded[, "to"]), to[expected_order][c(TRUE, diff(env$encode_transition(from, to, n_grid=150000L)[expected_order])!=0)])
})

test_that("legacy transition keys still decode with the old base", {
  legacy_key<-1*1e5+2

  decoded<-env$decode_transition(legacy_key)

  expect_equal(unname(decoded[, "from"]), 1L)
  expect_equal(unname(decoded[, "to"]), 2L)
})
