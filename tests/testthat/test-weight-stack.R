test_that("preallocated weight stack matches repeated cbind", {
   env<-new.env(parent=globalenv())
   source_checkout_r("run_particle_filter.R", env=env)

   columns<-list(
      c(0.2, 0.3, 0.5),
      c(1, 0, 2),
      c(1e-300, 2e-300, 3e-300)
   )

   old<-as.matrix(columns[[1]])
   values<-matrix(NA_real_, nrow=length(columns[[1]]), ncol=length(columns)+1L)
   values[,1]<-columns[[1]]
   active<-1L
   start<-1L

   expect_equal(env$weight_stack_matrix(values, start, active), old)
   expect_equal(env$weight_stack_last(values, start, active), old[,1])

   for (i in 2:length(columns)) {
      old<-cbind(old, columns[[i]])
      active<-active+1L
      values[, env$weight_stack_last_column(start, active, ncol(values))]<-columns[[i]]
      expect_equal(unname(env$weight_stack_matrix(values, start, active)), unname(old), tolerance=0)
      expect_equal(unname(env$weight_stack_last(values, start, active)), unname(old[,ncol(old)]), tolerance=0)
   }
})

test_that("preallocated weight stack row products match old matrix path", {
   env<-new.env(parent=globalenv())
   source_checkout_r("run_particle_filter.R", env=env)

   columns<-list(
      c(0.2, 0.3, 0.5),
      c(0.7, 0.8, 0.9),
      c(1e-100, 2e-100, 3e-100)
   )
   old<-as.matrix(columns[[1]])
   values<-matrix(NA_real_, nrow=3, ncol=4)
   values[,1]<-columns[[1]]
   active<-1L
   start<-1L
   for (i in 2:3) {
      old<-cbind(old, columns[[i]])
      active<-active+1L
      values[, env$weight_stack_last_column(start, active, ncol(values))]<-columns[[i]]
   }

   rowProds<-function(a) exp(rowSums(log(a)))
   expect_equal(rowProds(env$weight_stack_matrix(values, start, active)), rowProds(old), tolerance=0)

   current<-c(0.4, 0.5, 0.6)
   expect_equal(
      rowProds(cbind(env$weight_stack_matrix(values, start, active), current)),
      rowProds(cbind(old, current)),
      tolerance=0
   )
})

test_that("preallocated weight stack handles reset, set last, and drop first", {
   env<-new.env(parent=globalenv())
   source_checkout_r("run_particle_filter.R", env=env)

   values<-matrix(NA_real_, nrow=4, ncol=4)
   values[,1]<-rep(0.25, 4)
   active<-1L
   start<-1L

   for (col in list(c(1, 2, 3, 4), c(5, 6, 7, 8))) {
      active<-active+1L
      values[, env$weight_stack_last_column(start, active, ncol(values))]<-col
   }
   old<-cbind(rep(0.25, 4), c(1, 2, 3, 4), c(5, 6, 7, 8))
   expect_equal(unname(env$weight_stack_matrix(values, start, active)), unname(old), tolerance=0)

   values[, env$weight_stack_last_column(start, active, ncol(values))]<-rep(0.25, 4)
   old[, ncol(old)]<-rep(0.25, 4)
   expect_equal(unname(env$weight_stack_matrix(values, start, active)), unname(old), tolerance=0)

   start<-(start %% ncol(values))+1L
   active<-active-1L
   old<-old[, -1, drop=FALSE]
   expect_equal(unname(env$weight_stack_matrix(values, start, active)), unname(old), tolerance=0)

   values[,1]<-rep(0.25, 4)
   active<-1L
   start<-1L
   expect_equal(env$weight_stack_matrix(values, start, active), as.matrix(rep(0.25, 4)), tolerance=0)
})

test_that("preallocated weight stack preserves logical order after circular wrap", {
   env<-new.env(parent=globalenv())
   source_checkout_r("run_particle_filter.R", env=env)

   values<-matrix(NA_real_, nrow=3, ncol=4)
   values[,1]<-c(1, 1, 1)
   active<-1L
   start<-1L
   old<-as.matrix(c(1, 1, 1))

   for (i in 2:4) {
      col<-rep(i, 3)
      old<-cbind(old, col)
      active<-active+1L
      values[, env$weight_stack_last_column(start, active, ncol(values))]<-col
   }
   expect_equal(unname(env$weight_stack_matrix(values, start, active)), unname(old), tolerance=0)

   start<-(start %% ncol(values))+1L
   active<-active-1L
   old<-old[, -1, drop=FALSE]
   expect_equal(unname(env$weight_stack_matrix(values, start, active)), unname(old), tolerance=0)

   col<-c(5, 5, 5)
   old<-cbind(old, col)
   active<-active+1L
   values[, env$weight_stack_last_column(start, active, ncol(values))]<-col
   expect_equal(unname(env$weight_stack_matrix(values, start, active)), unname(old), tolerance=0)
   expect_equal(unname(env$weight_stack_last(values, start, active)), unname(old[, ncol(old)]), tolerance=0)

   start<-(start %% ncol(values))+1L
   active<-active-1L
   start<-(start %% ncol(values))+1L
   active<-active-1L
   old<-old[, -c(1, 2), drop=FALSE]
   expect_equal(unname(env$weight_stack_matrix(values, start, active)), unname(old), tolerance=0)
})
