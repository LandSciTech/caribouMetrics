#' Apply disturbances to land cover class data
#'
#' This sets the proportion of forest resource types to 0 when there is more
#' than 35% disturbance in a 16ha area. Natural and anthropogenic disturbance
#' will affect all resource types except water (LGW) and natural disturbance
#' (DTN). 
#' 
#' 
#' @param landCover 
#' @param anthroDist 
#' @param natDist 
#' 
#' @export

applyDist <- function(landCover, natDist, anthroDist, tmplt){
  # check anthroDist and natDist are real if not make dummy
  if(length(raster::unique(anthroDist)) == 0){
    anthroDist <- landCover
    anthroDist[] <- 0
  }
  if(length(raster::unique(natDist)) == 0){
    natDist <- landCover
    natDist[] <- 0
  }
  
  # Transfer DTN from natDist to landCover
  DTNcode <- resTypeCode %>%
    filter(ResourceType == "DTN") %>%
    pull(code)
  
  landCover <- mask(landCover, natDist, maskvalue = 1, 
                    updatevalue = DTNcode)

  # convert to 16 ha resolution stack of ResType proportion to match 16 ha
  # hexagons in Rempel
  #tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))
  
  landCover <- raster::layerize(landCover, classes = resTypeCode$code) %>% 
    raster::resample(tmplt, method = "bilinear")
  
  # if no distubance data provided return landCover
  if(length(raster::unique(natDist)) == 0 &&
     length(raster::unique(anthroDist)) == 0){
    return(landCover)
  }
  

  allDist16ha <- raster::stack(anthroDist, natDist) %>% 
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
    landCover[[i]] <- raster::mask(landCover[[i]], max(allDist16ha > 0.35),
                                   maskvalue = 1,
                                   updatevalue = 0)
  }
  
  return(landCover)
}