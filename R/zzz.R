.First.lib <- function(lib, pkg) {
  library.dynam( "som", pkg, lib )
  require(mva)
}
