#' Rasterize line density
#'
#' Rasterize line density in meters per hectare.
#'
#' @param x an sf object containing lines and/or points
#' @param r a SpatRaster or RasterLayer object to be used as a template for the
#'   output raster
#' @param ptDensity a number giving the density to assign to points, in units of
#'   `res(r)`. A value of 1 indicates one straight line crossing of the pixel. A
#'   value of 2+2*2^0.5 is horizontal, vertical, and diagonal crossings. If
#'   NULL, points in x will be ignored.
#'
#' @return A SpatRaster object with values representing the density of lines in
#'   meters per hectare.
#' @examples
#' # create example raster
#' lc <- terra::rast(xmin = 0, xmax = 25000, ymin = 0, ymax = 25000, 
#'                      resolution = 250, crs = "EPSG:5070")
#'
#' #' # create line
#' lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10000, 10000),
#'                                                             ncol = 2, byrow = TRUE)),
#'                                    sf::st_linestring(matrix(c(1, 10001, 10001, 1),
#'                                                             ncol = 2, byrow = TRUE)),
#'                                    sf::st_linestring(matrix(c(5001, 10001, 5001, 1),
#'                                                             ncol = 2, byrow = TRUE))),
#'                                    crs = 5070))
#'
#' rastLines <- rasterizeLineDensity(lf, lc)
#'
#' plot(rastLines)
#'
#' @family habitat
#' @export
#' 
rasterizeLineDensity <- function(x, r, ptDensity = 1) {
  if(any(c("POINT", "MULTIPOINT") %in% 
         sf::st_geometry_type(x, by_geometry = TRUE))){
    lfPt <- sf::st_collection_extract(x, "POINT")
    x <- sf::st_collection_extract(x, "LINESTRING")
  } else {
    lfPt <- slice(x, 0)
  }
  
  line_len <- terra::rasterizeGeom(terra::vect(x), r, fun = "length")
  
  cell_area <- terra::cellSize(r)/10000 
  
  r <- round(line_len/cell_area, digits = 1)

  
  if(!is.null(ptDensity)){
    if(ptDensity > 2+2*2^0.5){
      warning("ptDensity is greater than the expected max of 4.828.",
              " see ?rasterizeLineDensity for details",
              call. = FALSE)
    }
      
    if(nrow(lfPt) > 0){
      lfR <- terra::rasterizeGeom(terra::vect(lfPt), r, fun = "count")    
      
      lfR <- lfR * ptDensity * terra::res(r)[1]
      
      lfR <- round(lfR / (terra::res(r)[1] * terra::res(r)[2] / 10000), digits = 1)
      r <- r + lfR
    }
  }
  
  return(r)
}
