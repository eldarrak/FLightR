my_google_api_key <- function() {
  val <- Sys.getenv("google_api_key")
  if (identical(val, "")) {
    cat("`google_api_key` env var has not been set")
  } else {
    return(val)
  }
}

google_api_key<-my_google_api_key()

if (!is.null(a)) {
   ggmap::register_google(google_api_key)
}