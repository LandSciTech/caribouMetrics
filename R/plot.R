#' @include AAAClassDefinitions.R
NULL

# Overwrite default plot methods to improve display of CaribouHabitat objects.

#' Plot method for CaribouHabitat Objects
#'
#' This will plot the predicted habitat use from a CaribouHabitat object.
#'
#' @param x A CaribouHabitat object
#' @param season character. By default "all" or supply a vector of seasons ie
#'   \code{c("Spring", "Summer", "Fall", "Winter")}
#' @param raster.title character. Title to give map
#' @param tmap logical. Should tmap be used for plotting by default it is used
#'   if installed
#' @param ... Other agruments passed to \code{tmap::qtm} or \code{raster::plot}
#'
#' @return If tmap is TRUE a tmap object if FALSE a plot is created in the viewer
#' @export
#'
setMethod("plot", "CaribouHabitat", 
          function(x, season = "all", raster.title = "Probability\nof use",
                   tmap = requireNamespace("tmap", quietly = TRUE), ...) {
  if (season[1] == "all"){
    season <- c("Spring", "Summer", "Fall", "Winter")
  }
  if(tmap){
    tmap::qtm(x@habitatUse[[season]], raster.title = raster.title, ...)+
      tmap::qtm(x@projectPoly, fill = NULL)
  }else {
    plot(x@habitatUse[[season]], ...)
  }
  
})
