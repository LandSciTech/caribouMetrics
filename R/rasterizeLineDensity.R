#' Rasterize line density
#'
#' Rasterize line density in meters per hectare.
#'
#' @param x an sf object containing lines and/or points
#' @param r a RasterLayer object to be used as a template for the output raster
#' @param ptDensity a number giving the density to assign to points, in units of
#'   res(r). A value of 1 indicates one straight line crossing of the pixel. A
#'   value of 2+2*2^0.5 is horizontal, vertical, and diagonal crossings. If
#'   NULL, points in linObj will be ignored.
#'
#' @return A RasterLayer object with values representing the density of lines
#'   per hectare.
#' @export
#' 
rasterizeLineDensity <- function(x, r, ptDensity = 1) {


  r[] <- 1:ncell(r)
  
  rPoly <- spex::polygonize(r) %>% set_names("ID", "geometry") %>% 
    st_set_agr("constant")
  
  rp2 <- st_intersection(rPoly, st_set_agr(x, "constant")) %>% 
    mutate(length = st_length(geometry) %>% units::set_units(NULL)) %>% 
    select(ID, length, geometry) %>% st_drop_geometry() %>% 
    group_by(ID) %>% 
    summarise(length = round(sum(length, na.rm = TRUE)/(res(r)[1]*res(r)[2]/10000), digits = 1))
  
  rp2 <- left_join(rPoly %>% st_drop_geometry(), rp2, by = "ID") %>% 
    mutate(length = replace_na(length, 0))
  
  r[] <- rp2$length
 
  lfPt <- x %>% dplyr::filter(st_is(. , "POINT"))
  
  if(!is.null(ptDensity)){
    if(ptDensity > 2+2*2^0.5){
      warning("ptDensity is greater than the expected max of 4.828.",
              " see ?rasterizeLineDensity for details",
              call. = FALSE)
    }
      
    if(nrow(lfPt)>0){
      lfR <- raster::rasterize(lfPt, r, field = "linFID")    
      
      lfR[!is.na(lfR)] <- ptDensity * res(r)[1]
      
      lfR[is.na(lfR)] <- 0
      lfR <- round(lfR / (res(r)[1] * res(r)[2] / 10000), digits = 1)
      r <- r + lfR
    }
  }
  
  return(r)
}