#' Apply disturbances to land cover class data
#'
#' This sets the proportion of forest resource types to 0 when there is more
#' than 35% disturbance in a 16ha area. Natural and anthropogenic disturbance
#' will affect all resource types except water (LGW) and natural disturbance
#' (DTN). Harvest disturbance will only affect forst resource types and in
#' addition to water and natural disturbance does not affect wetlands (ST, LGTP,
#' LGOP)
#' 
#' 
#' @param landCover 
#' @param harv 
#' @param anthroDist 
#' @param natDist 
#' 
#' @export

applyDist <- function(landCover, natDist, anthroDist, harv){
  # convert to 16 ha resolution stack of ResType proportion to match 16 ha
  # hexagons in Rempel
  tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))
  
  landCover <- raster::layerize(landCover, classes = resTypeCode$code) %>% 
    raster::resample(tmplt, method = "bilinear")
  
  allDist16ha <- raster::stack(harv, anthroDist, 
                               natDist) %>% 
    raster::resample(tmplt, method = "bilinear")
  
  # Get proportion of land in 16 ha area that has each type of disturbance
  watCode <- resTypeCode %>% 
    filter(ResourceType == "LGW") %>% 
    pull(code)
  
  # get proportion land
  land <- 1 - landCover[[watCode]]
  
  # divide prop disturbance by prop land 
  propLandDist <- Vectorize(function(dist, land){
    ifelse(land < dist, 1, 
           ifelse(land == 0, 0, dist/land))
  })
  
  allDist16ha <- raster::overlay(allDist16ha, land, fun = propLandDist)
  
  
  # make landCover have 0 forest classes when 16 ha area disturbed > 0.35 by natural
  # disturbance or anthropogenic disturbance
  toChange <- resTypeCode %>% filter(!ResourceType %in% c("DTN", "LGW")) %>% 
    pull(code)
  
  for (i in toChange) {
    landCover[[i]] <- raster::mask(landCover[[i]], max(allDist16ha[[2:3]] > 0.35),
                                   maskvalue = 1,
                                   updatevalue = 0)
  }
  
  
  # make landCover have 0 forest classes when 16 ha area disturbed > 0.35 by harvest
  # disturbance but don't change wetlands
  toChange <- resTypeCode %>% 
    filter(!ResourceType %in% c("DTN", "LGW", "LGTP", "LGOP", "ST")) %>% 
    pull(code)
  
  for (i in toChange) {
    landCover[[i]] <- raster::mask(landCover[[i]], allDist16ha[[1]] > 0.35,
                                   maskvalue = 1, 
                                   updatevalue = 0)
  }
  
  return(landCover)
}