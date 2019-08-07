.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to receptoR!\nTo launch the shiny application, use 'launchApp()'")
  shiny::addResourcePath('www', system.file("shiny/www", package = "receptoR"))
}