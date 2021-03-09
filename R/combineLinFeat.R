
#' Combine linear features
#' 
#' Combine roads, rail and utilities in to one linear features object.
#'
#' @param roads 
#' @param rail 
#' @param utilities 
#'
#' @return sf object
#' @export
#'
#' @examples
combineLinFeat <- function(roads, rail, utilities){
  
  if(is.character(roads)){
    roads <- st_read(roads, quiet = TRUE, agr = "constant")
  }
  if(is.character(rail)){
    rail <- st_read(rail, quiet = TRUE, agr = "constant")
  }
  if(is.character(utilities)){
    utilities <- st_read(utilities, quiet = TRUE, agr = "constant")
  }
  
  if(is(roads, "Raster")){
      roads <- raster::rasterToPoints(roads, fun = function(x){x > 0}, 
                                        spatial = TRUE) %>% 
        sf::st_as_sf()
  }
  
  roads <- roads %>% transmute(ID = 1, Type = "road") 
  
  utilities <- utilities %>% 
    transmute(ID = 1, Type = "utility") %>% 
    st_transform(st_crs(roads))
  
  rail <- rail %>% 
    transmute(ID = 1, Type = "rail") %>% 
    st_transform(st_crs(roads))
  
  rbind(roads, utilities, rail) %>% mutate(ID = 1:n()) %>% 
    st_set_agr(st_agr(roads))
}