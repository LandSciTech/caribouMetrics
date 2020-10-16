#' @include AAAClassDefinitions.R
NULL

# Overwrite default plot methods to improve display of CaribouHabitat objects.

setMethod("plot", "CaribouHabitat", 
          function(x, season = "all", 
                   tmap = requireNamespace("tmap", quietly = TRUE), ...) {
  if (season[1] == "all"){
    season <- c("Spring", "Summer", "Fall", "Winter")
  }
  if(tmap){
    tmap::qtm(x@habitatUse[[season]], ...)
  }else {
    plot(x@habitatUse[[season]], ...)
  }
  
})
