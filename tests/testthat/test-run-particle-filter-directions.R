test_that("directional proposal uses bearings from current point to candidates", {
   env<-new.env(parent=globalenv())
   sys.source(testthat::test_path("../../R/run_particle_filter.R"), envir=env)

   Grid<-matrix(c(0, 0, 1,
                  1, 0, 1,
                  0, 1, 1),
                ncol=3, byrow=TRUE)
   colnames(Grid)<-c("lon", "lat", "mask")
   in.Data<-list(Spatial=list(Grid=Grid))

   expect_equal(as.numeric(env$dir_fun(c(1, 2), in.Data)), 90, tolerance=1)
   expect_equal(as.numeric(env$dir_fun(c(1, 3), in.Data)), 0, tolerance=1)

   Current.Proposal<-list(
      M.mean=111,
      M.sd=20,
      Direction=90,
      Kappa=10
   )

   set.seed(1)
   proposed<-env$generate.points.dirs(
      x=c(1, 200, 200),
      in.Data=in.Data,
      Current.Proposal=Current.Proposal,
      a=0,
      b=500
   )

   expect_gt(mean(proposed == 2), 0.95)
   expect_lt(mean(proposed == 3), 0.05)
})
