#' @include AAAClassDefinitions.R
NULL

# Overwrite default plot methods to improve display of CaribouHabitat objects.

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
