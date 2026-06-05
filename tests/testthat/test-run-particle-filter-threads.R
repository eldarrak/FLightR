test_that("run.particle.filter defines Threads for threads=1", {
   env<-new.env(parent=globalenv())
   sys.source(testthat::test_path("../../R/run_particle_filter.R"), envir=env)

   call.args<-NULL
   env$pf.run.parallel.SO.resample<-function(...) {
      call.args<<-list(...)
      stop("stubbed particle filter", call.=FALSE)
   }

   expect_error(
      env$run.particle.filter(list(), threads=1, nParticles=1, plot=FALSE),
      "stubbed particle filter"
   )
   expect_equal(call.args$threads, 1)
   expect_false(call.args$parallel)
})
