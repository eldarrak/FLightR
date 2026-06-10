test_that("cached movement candidates respect distance bounds", {
  testthat::skip_if_not_installed("geosphere")
  env<-new.env(parent=globalenv())
  sys.source(testthat::test_path("../../R/run_particle_filter.R"), envir=env)

  Grid<-matrix(c(0, 0, 1,
                 1, 0, 1,
                 0, 1, 1,
                 10, 0, 1),
               ncol=3, byrow=TRUE)
  colnames(Grid)<-c("lon", "lat", "mask")

  candidates<-env$build.grid.movement.candidates(Grid, a=50, b=200)
  expect_equal(candidates$grid_n, 4)
  expect_true(all(candidates$to[[1]] %in% c(2L, 3L)))
  expect_true(all(candidates$distance[[1]]>=50))
  expect_true(all(candidates$distance[[1]]<=200))
})

test_that("cached propagation returns valid indices and preserves non-moving particles", {
  testthat::skip_if_not_installed("geosphere")
  env<-new.env(parent=globalenv())
  sys.source(testthat::test_path("../../R/run_particle_filter.R"), envir=env)

  Grid<-matrix(c(0, 0, 1,
                 1, 0, 1,
                 0, 1, 1),
               ncol=3, byrow=TRUE)
  colnames(Grid)<-c("lon", "lat", "mask")
  in.Data<-list(Spatial=list(Grid=Grid, tmp=list()))
  in.Data$Spatial$tmp$movement_candidates<-env$build.grid.movement.candidates(Grid, a=50, b=200)
  Current.Proposal<-list(M.mean=111, M.sd=20, Direction=90, Kappa=0)

  set.seed(10)
  proposed<-env$generate.points.dirs(
    x=c(1, 100, 25),
    in.Data=in.Data,
    Current.Proposal=Current.Proposal,
    a=50,
    b=200
  )

  expect_length(proposed, 100)
  expect_true(all(proposed %in% seq_len(nrow(Grid))))
  expect_equal(sum(proposed == 1), 75)
  expect_true(all(proposed[proposed != 1] %in% c(2L, 3L)))
})

test_that("cached distance-only proposal is close to legacy on a simple grid", {
  testthat::skip_if_not_installed("geosphere")
  testthat::skip_if_not_installed("sf")
  testthat::skip_if_not_installed("truncnorm")
  testthat::skip_if_not_installed("circular")
  env<-new.env(parent=globalenv())
  sys.source(testthat::test_path("../../R/run_particle_filter.R"), envir=env)

  Grid<-matrix(c(0, 0, 1,
                 1, 0, 1,
                 0, 1, 1),
               ncol=3, byrow=TRUE)
  colnames(Grid)<-c("lon", "lat", "mask")
  Current.Proposal<-list(M.mean=111, M.sd=25, Direction=0, Kappa=0)

  legacy.Data<-list(Spatial=list(Grid=Grid, tmp=list()))
  cached.Data<-legacy.Data
  cached.Data$Spatial$tmp$movement_candidates<-env$build.grid.movement.candidates(Grid, a=50, b=200)

  set.seed(123)
  legacy<-env$generate.points.dirs.legacy(c(1, 5000, 5000), legacy.Data, Current.Proposal, a=50, b=200)
  set.seed(123)
  cached<-env$generate.points.dirs.cached(c(1, 5000, 5000), cached.Data, Current.Proposal, a=50, b=200)

  legacy.freq<-tabulate(legacy, nbins=3)/length(legacy)
  cached.freq<-tabulate(cached, nbins=3)/length(cached)
  expect_lt(max(abs(legacy.freq-cached.freq)), 0.05)
})

test_that("cached directional proposal favors intended direction", {
  testthat::skip_if_not_installed("geosphere")
  env<-new.env(parent=globalenv())
  sys.source(testthat::test_path("../../R/run_particle_filter.R"), envir=env)

  Grid<-matrix(c(0, 0, 1,
                 1, 0, 1,
                 0, 1, 1),
               ncol=3, byrow=TRUE)
  colnames(Grid)<-c("lon", "lat", "mask")
  in.Data<-list(Spatial=list(Grid=Grid, tmp=list()))
  in.Data$Spatial$tmp$movement_candidates<-env$build.grid.movement.candidates(Grid, a=50, b=200)
  Current.Proposal<-list(M.mean=111, M.sd=20, Direction=90, Kappa=10)

  set.seed(1)
  proposed<-env$generate.points.dirs.cached(
    x=c(1, 500, 500),
    in.Data=in.Data,
    Current.Proposal=Current.Proposal,
    a=50,
    b=200
  )

  expect_gt(mean(proposed == 2), 0.95)
  expect_lt(mean(proposed == 3), 0.05)
})

test_that("legacy backend and run.particle.filter backend argument remain available", {
  testthat::skip_if_not_installed("sf")
  testthat::skip_if_not_installed("geosphere")
  testthat::skip_if_not_installed("truncnorm")
  testthat::skip_if_not_installed("circular")
  env<-new.env(parent=globalenv())
  sys.source(testthat::test_path("../../R/run_particle_filter.R"), envir=env)

  Grid<-matrix(c(0, 0, 1,
                 1, 0, 1,
                 0, 1, 1),
               ncol=3, byrow=TRUE)
  colnames(Grid)<-c("lon", "lat", "mask")
  in.Data<-list(Spatial=list(Grid=Grid, tmp=list()))
  Current.Proposal<-list(M.mean=111, M.sd=20, Direction=90, Kappa=10)
  set.seed(1)
  proposed<-env$generate.points.dirs.legacy(c(1, 20, 20), in.Data, Current.Proposal, a=50, b=200)
  expect_true(all(proposed %in% seq_len(nrow(Grid))))

  call.args<-NULL
  env$pf.run.parallel.SO.resample<-function(...) {
    call.args<<-list(...)
    stop("stubbed particle filter", call.=FALSE)
  }
  for (backend in c("legacy", "cached", "auto")) {
    expect_error(
      env$run.particle.filter(list(), threads=1, nParticles=1, plot=FALSE, propagation.backend=backend),
      "stubbed particle filter"
    )
    expect_equal(call.args$propagation.backend, backend)
  }
})
