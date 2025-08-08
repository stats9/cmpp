.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  packageStartupMessage(
    paste("Welcome to the cmpp package | Version:", v)
  )
}