#' Apply disturbances to land cover class data
#'
#' This sets the proportion of forest resource types to 0 when there is more
#' than 35% disturbance in a 16ha area. Natural and anthropogenic disturbance
#' will affect all resource types except water (LGW) and natural disturbance
#' (DTN). 
#' 
#' 
#' @param landCover landCover raster
#' @param anthroDist binary raster
#' @param natDist binary raster
#' @param tmplt template
#' 
#' @noRd


applyDist <- function(landCover, natDist, anthroDist, tmplt){
  # check anthroDist and natDist are real if not make dummy
  anthroDummy <- terra::ncell(anthroDist) == 1
  natDummy <- terra::ncell(natDist) == 1
  
  if(anthroDummy){
    anthroDist <- terra::init(
      landCover, fun = function(x){rep(0, x)}, 
      filename = tempfile(pattern = "spat", 
                          tmpdir = terra::terraOptions(print = FALSE)$tempdir, 
                          fileext = ".grd"))
  }
  if(natDummy){
    natDist <- terra::init(
      landCover, fun = function(x){rep(0, x)}, 
      filename = tempfile(pattern = "spat", 
                          tmpdir = terra::terraOptions(print = FALSE)$tempdir, 
                          fileext = ".grd"))
  }
  
  # Transfer DTN from natDist to landCover
  DTNcode <- resTypeCode %>%
    filter(.data$ResourceType == "DTN") %>%
    pull(.data$code)
  
  landCover <- terra::mask(landCover, natDist, maskvalue = 1, 
                           updatevalue = DTNcode)
  
  # convert to 16 ha resolution stack of ResType proportion to match 16 ha
  # hexagons in Rempel
  #tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))
  
  landCover <- terra::segregate(landCover, classes = resTypeCode$code) %>% 
    terra::resample(tmplt, method = "bilinear")
  
  # if no distubance data provided return landCover
  if(anthroDummy && natDummy){
    return(landCover)
  }
  
  allDist16ha <- c(anthroDist, natDist) %>% 
    terra::resample(tmplt, method = "bilinear")
  rm(anthroDist, natDist)
  
  # Get proportion of land in 16 ha area that has each type of disturbance
  watCode <- resTypeCode %>% 
    filter(.data$ResourceType == "LGW") %>% 
    pull(.data$code)
  
  # get proportion land
  land <- 1 - landCover[[watCode]]
  
  # divide prop disturbance by prop land 
  propLandDist <- function(dist, land){
    ifelse(land < dist, 1, 
           ifelse(land == 0, 0, dist/land))
  }
  
  allDist16ha <- c(
    terra::lapp(c(allDist16ha[[1]], land), fun = propLandDist),
    terra::lapp(c(allDist16ha[[2]], land), fun = propLandDist)
  )
  
  # make landCover have 0 forest classes when 16 ha area disturbed > 0.35 by natural
  # disturbance or anthropogenic disturbance  
  anyDist35 <- max(allDist16ha > 0.35)
  rm(allDist16ha)
  toChange <- resTypeCode %>% filter(!.data$ResourceType %in% c("DTN", "LGW")) %>% 
    pull(.data$code)
  
  for (i in toChange) {
    landCover[[i]] <- raster::mask(landCover[[i]], anyDist35,
                                   maskvalue = 1,
                                   updatevalue = 0)
  }
  
  return(landCover)
}