
get.transition.rle=function(From, To) {
  rle(sort.int(From*1e5+To, method = "quick"))
}