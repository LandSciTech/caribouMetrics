
#' Functions to align spatial objects
#'
#' Helper functions for cropping, transforming and checking the overlap of
#' spatial objects
#'
#'
#' @noRd

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
    return(x)
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
      return(raster::crop(x, y, snap = "out"))
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
  return(x)
}

checkCompRast <- function(x, y, nmx, nmy, y2 = NULL){
  chk1 <- raster::compareRaster(x, y, stopiffalse = FALSE)
  
  if(!is.null(y2)){
    chk2 <- raster::compareRaster(x, y2, stopiffalse = FALSE)
  } else {
    chk2 <- FALSE
  }

  if(all(!chk1, !chk2)){
    stop(nmx, " and ",  nmy,
         " rasters do not have the same extent,",
         " number of rows and columns, projection, resolution, or origin. ",
         "Use raster::compareRaster() to identify the problem.", call. = FALSE)
  }

  invisible(TRUE)
}

