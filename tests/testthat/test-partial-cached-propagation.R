test_that("partial cached candidates match full cached candidates on a small grid", {
  testthat::skip_if_not_installed("sf")
  testthat::skip_if_not_installed("geosphere")
  env<-new.env(parent=globalenv())
  source_checkout_r("run_particle_filter.R", env=env)

  Grid<-matrix(c(0, 0, 1,
                 1, 0, 1,
                 0, 1, 1,
                 3, 0, 1,
                 10, 0, 1),
               ncol=3, byrow=TRUE)
  colnames(Grid)<-c("lon", "lat", "mask")

  full<-env$build.grid.movement.candidates(Grid, a=50, b=250)
  cache<-env$build.partial.movement.cache(Grid, a=50, b=250)
  entry<-env$get.partial.movement.entry(cache, 1L)

  expect_equal(sort(entry$to), sort(full$to[[1]]))
  expect_equal(length(entry$distance), length(full$distance[[1]]))
  expect_equal(length(entry$bearing), length(full$bearing[[1]]))
  expect_equal(cache$misses, 1L)
  expect_equal(cache$hits, 0L)

  entry2<-env$get.partial.movement.entry(cache, 1L)
  expect_equal(entry2$to, entry$to)
  expect_equal(cache$hits, 1L)
})

test_that("partial cached propagation matches full cached sampling on a small grid", {
  testthat::skip_if_not_installed("sf")
  testthat::skip_if_not_installed("geosphere")
  testthat::skip_if_not_installed("truncnorm")
  testthat::skip_if_not_installed("circular")
  env<-new.env(parent=globalenv())
  source_checkout_r("run_particle_filter.R", env=env)

  Grid<-matrix(c(0, 0, 1,
                 1, 0, 1,
                 0, 1, 1,
                 10, 0, 1),
               ncol=3, byrow=TRUE)
  colnames(Grid)<-c("lon", "lat", "mask")
  Current.Proposal<-list(M.mean=111, M.sd=20, Direction=90, Kappa=10)

  full.Data<-list(Spatial=list(Grid=Grid, tmp=list()))
  full.Data$Spatial$tmp$movement_candidates<-env$build.grid.movement.candidates(Grid, a=50, b=200)
  partial.Data<-list(Spatial=list(Grid=Grid, tmp=list()))
  partial.Data$Spatial$tmp$partial_movement_cache<-env$build.partial.movement.cache(Grid, a=50, b=200)

  set.seed(42)
  full<-env$generate.points.dirs.cached(c(1, 500, 425), full.Data, Current.Proposal, a=50, b=200)
  set.seed(42)
  partial<-env$generate.points.dirs.partial(c(1, 500, 425), partial.Data, Current.Proposal, a=50, b=200)

  expect_equal(partial, full)
  expect_equal(sum(partial == 1), 75)
  expect_true(all(partial %in% seq_len(nrow(Grid))))
})

test_that("partial cached diagnostics show lazy origin caching", {
  testthat::skip_if_not_installed("sf")
  testthat::skip_if_not_installed("geosphere")
  env<-new.env(parent=globalenv())
  source_checkout_r("run_particle_filter.R", env=env)

  Grid<-matrix(c(0, 0, 1,
                 1, 0, 1,
                 0, 1, 1,
                 2, 0, 1,
                 10, 0, 1),
               ncol=3, byrow=TRUE)
  colnames(Grid)<-c("lon", "lat", "mask")
  cache<-env$build.partial.movement.cache(Grid, a=50, b=250)

  env$get.partial.movement.entry(cache, 1L)
  env$get.partial.movement.entry(cache, 1L)
  diag<-env$partial.movement.cache.diagnostics(cache)

  expect_equal(diag$propagation_backend, "partial_cached")
  expect_equal(diag$grid_size, nrow(Grid))
  expect_equal(diag$cached_origins, 1L)
  expect_lt(diag$cached_origins, nrow(Grid))
  expect_equal(diag$cache_misses, 1L)
  expect_equal(diag$cache_hits, 1L)
  expect_gt(diag$total_candidate_entries, 0)
  expect_gt(diag$approximate_cache_object_bytes, 0)
})

test_that("partial cached backend argument is accepted by run.particle.filter", {
  env<-new.env(parent=globalenv())
  source_checkout_r("run_particle_filter.R", env=env)

  call.args<-NULL
  env$pf.run.parallel.SO.resample<-function(...) {
    call.args<<-list(...)
    stop("stubbed particle filter", call.=FALSE)
  }

  expect_error(
    env$run.particle.filter(list(), threads=1, nParticles=1, plot=FALSE, propagation.backend="partial_cached"),
    "stubbed particle filter"
  )
  expect_equal(call.args$propagation.backend, "partial_cached")
})