load_run_pf_env <- function() {
  env<-new.env(parent=globalenv())
  source_checkout_r("run_particle_filter.R", env=env)
  env
}

stub_run_input <- function() {
  list(
    Indices=list(
      Main.Index=data.frame(x=1),
      Matrix.Index.Table=data.frame(time=as.POSIXct("2020-01-01", tz="UTC"))
    ),
    Spatial=list(tmp=list()),
    Results=list()
  )
}

install_deterministic_pf_stubs <- function(env) {
  env$pf.run.parallel.SO.resample<-function(..., verbose=FALSE) {
    env$pf_message("normal progress", verbose=verbose)
    env$pf_debug("debug detail", verbose=verbose)
    list(
      Points=list(
        structure(list(lengths=1L, values=1L), class="rle"),
        structure(list(lengths=1L, values=1L), class="rle")
      ),
      Trans=list(),
      Results=list(
        outliers=integer(),
        tmp.results=list(),
        partial_cache_diagnostics=NULL
      )
    )
  }
  env$get.LL.PF<-function(...) -12.34
  env$get.coordinates.PF<-function(Points, in.Data, add.jitter=FALSE, verbose=FALSE) {
    env$pf_message("coordinates progress", verbose=verbose)
    in.Data$Results$Quantiles<-data.frame(time=as.POSIXct("2020-01-01", tz="UTC"))
    in.Data$Results$Points.rle<-Points
    in.Data
  }
  env$estimate.movement.parameters<-function(..., verbose=FALSE) {
    env$pf_message("movement progress", verbose=verbose)
    list(Movement.results=data.frame(value=1), Transitions.rle=list())
  }
}

test_that("run.particle.filter default is quiet and passes verbose FALSE to PF core", {
  env<-load_run_pf_env()
  call.args<-NULL
  env$pf.run.parallel.SO.resample<-function(...) {
    call.args<<-list(...)
    stop("stubbed particle filter", call.=FALSE)
  }

  expect_error(
    expect_no_message(
      env$run.particle.filter(list(), threads=1, nParticles=1, plot=FALSE)
    ),
    "stubbed particle filter"
  )
  expect_false(call.args$verbose)
})

test_that("verbose TRUE prints normal progress", {
  env<-load_run_pf_env()
  env$pf.run.parallel.SO.resample<-function(..., verbose=FALSE) {
    env$pf_message("normal progress", verbose=verbose)
    stop("stubbed particle filter", call.=FALSE)
  }

  expect_message(
    expect_error(
      env$run.particle.filter(list(), threads=1, nParticles=1, plot=FALSE, verbose=TRUE),
      "stubbed particle filter"
    ),
    "normal progress"
  )
})

test_that("verbose debug prints debug messages", {
  env<-load_run_pf_env()
  expect_no_message(env$pf_debug("debug detail"))
  expect_no_message(env$pf_debug("debug detail", verbose=TRUE))
  expect_no_message(env$pf_debug("debug detail", verbose=FALSE))
  expect_message(env$pf_debug("debug detail", verbose="debug"), "debug detail")
})

test_that("warnings are not suppressed by quiet default", {
  env<-load_run_pf_env()
  env$pf.run.parallel.SO.resample<-function(...) {
    stop("stubbed particle filter", call.=FALSE)
  }

  expect_warning(
    expect_error(
      env$run.particle.filter(list(), cpus=1, nParticles=1, plot=FALSE),
      "stubbed particle filter"
    ),
    "use threads instead of cpus"
  )
})

test_that("LL and returned results are unchanged between verbose modes", {
  testthat::skip_if_not_installed("FLightR")
  env<-load_run_pf_env()
  install_deterministic_pf_stubs(env)

  run_one<-function(verbose) {
    out<-testthat::capture_messages(
      res<-env$run.particle.filter(
        stub_run_input(),
        threads=1,
        nParticles=1,
        plot=FALSE,
        verbose=verbose
      )
    )
    list(result=res, output=out)
  }

  quiet<-run_one(FALSE)
  normal<-run_one(TRUE)
  debug<-run_one("debug")

  expect_identical(quiet$result$Results$LL, normal$result$Results$LL)
  expect_identical(quiet$result$Results$LL, debug$result$Results$LL)
  expect_equal(quiet$result$Results$Quantiles, normal$result$Results$Quantiles)
  expect_equal(quiet$result$Results$Quantiles, debug$result$Results$Quantiles)

  expect_length(quiet$output, 0)
  expect_true(any(grepl("normal progress", normal$output)))
  expect_true(any(grepl("debug detail", debug$output)))
})
