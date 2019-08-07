#' receptoR
#' This stand-alone package allows for analysis of high-throughput transcriptome data, specifically identifying molecular receptors present on specific cell types of interest.
#'
#' @export launchApp
#'
#' @return shiny application object
#'
#' @example \dontrun {launchApp()}
#'
#' @import shiny
#'


# wrapper for shiny::shinyApp()
launchApp <- function() {
    appDir <- system.file("shiny", package="receptoR")
  shiny::runApp(appDir)
}
