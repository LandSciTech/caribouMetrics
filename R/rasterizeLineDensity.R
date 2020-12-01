#' Rasterize line density
#'
#' Rasterize line density in meters per hectare.
#'
#' @param linObj an sf object containing lines
#' @param r a RasterLayer object to be used as a template for the output raster
#'
#' @return A RasterLayer object with values representing the density of lines
#'   per hectare.
#' @export
#'
rasterizeLineDensity <- function(x, r) {
  r[] <- 1:ncell(r)
  
  rPoly <- spex::polygonize(r) %>% set_names("ID", "geometry") %>% 
    st_set_agr("constant")
  
  rp2 <- st_intersection(rPoly, st_set_agr(x, "constant")) %>% 
    mutate(length = st_length(geometry) %>% units::drop_units()) %>% 
    select(ID, length, geometry) %>% st_drop_geometry() %>% 
    group_by(ID) %>% 
    summarise(length = round(sum(length, na.rm = TRUE)/(res(r)[1]*res(r)[2]/10000), digits = 1))
  
  rp2 <- left_join(rPoly %>% st_drop_geometry(), rp2, by = "ID") %>% 
    mutate(length = replace_na(length, 0))
  
  r[] <- rp2$length
  
  return(r)
}