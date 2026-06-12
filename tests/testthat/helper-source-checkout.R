source_checkout_r<-function(file, env=parent.frame()) {
  path<-testthat::test_path("..", "..", "R", file)
  testthat::skip_if_not(
    file.exists(path),
    paste("Source-checkout-only internal helper test:", file)
  )
  sys.source(path, envir=env)
}
