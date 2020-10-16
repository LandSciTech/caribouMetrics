
#' Functions to align spatial objects
#'
#' Helper functions for cropping, transforming and checking the overlap of
#' spatial objects
#'
#' @param x
#' @param y
#' @param nmx
#' @param nmy
#'
#' @return
#' @export
#'
#' @examples

aggregateIf <- function(x, y, nmx, nmy){
  if(inherits(x, "sf")){
    return(x)
  }
  if(inherits(x,"RasterLayer")){
    if(all(raster::res(x) == raster::res(y))){
      return(x)
    } else {
      message(nmx, " being aggregated to have resolution matching ", nmy,
              " using the mode")
      return(raster::aggregate(x, fact = res(y)[1]/res(x)[1],
                               fun = raster::modal))
    }
  }
  
}

transformIf <- function(x, y, nmx, nmy){
  if(inherits(x, "sf")){
    if(st_crs(x) == st_crs(y)){
      return(x)
    } else {
     message(nmx, " being transformed to have crs matching ", nmy)
     return(st_transform(x, st_crs(y)))
    }
  }
  if(inherits(x,"RasterLayer")){
    if(st_crs(x) == st_crs(y)){
      return(x)
    } else {
      message(nmx, " being transformed to have crs matching ", nmy,
              " using nearest neighbour")
      return(raster::projectRaster(x, to = y, method = "ngb"))
    }
  }
}

cropIf <- function(x, y, nmx, nmy){
  if(inherits(x, "sf")){
    if(raster::extent(x) != raster::extent(y)){
      message("cropping ", nmx, " to extent of ", nmy)
      return(st_crop(x, y))
    } else {
      return(x)
    }
  }
  if(inherits(x,"RasterLayer")){
    if(!compareRaster(x, y, extent = TRUE, rowcol = FALSE, crs = FALSE, 
                      orig = FALSE, stopiffalse = FALSE)){
      message("cropping ", nmx, " to extent of ", nmy)
      return(raster::crop(x, y))
    } else {
      return(x)
    }
  }
}

extendIf <- function(x, y, nmx, nmy){
  if(inherits(x, "sf")){
    return(x)
  }
  if(inherits(x,"RasterLayer")){
    if(!compareRaster(x, y, extent = TRUE, rowcol = FALSE, crs = FALSE, 
                     orig = FALSE, stopiffalse = FALSE)){
      message("extending ", nmx, " to extent of ", nmy)
      return(raster::extend(x, y))
    } else {
      return(x)
    }
  }
}

checkOverlap <- function(x,y, nmx, nmy){
  if(!st_intersects(st_as_sfc(st_bbox(x)) %>% st_transform(st_crs(y)), 
                    st_as_sfc(st_bbox(y)),
                    sparse = FALSE)[1,1]){
    stop(nmx, " does not overlap with ", nmy, call. = FALSE)
  } else {
    return(x)
  }
}


checkAlign <- function(x, y, nmx, nmy){
  x <- transformIf(x, y, nmx, nmy)
  x <- checkOverlap(x, y, nmx, nmy)
  x <- cropIf(x, y, nmx, nmy)
  x <- extendIf(x, y, nmx, nmy)
  return(x)
}
