

#' Update provincial land cover
#'
#' Update provincial land cover based on FRI, age and disturbance history
#'
#' @param plc 
#' @param fri 
#' @param age 
#' @param natDist 
#' @param anthroDist 
#' @param harv 
#'
#' @return
#' @export
#'
#' @examples
#' 
setGeneric("updatePLC", function(x, newData, ...) standardGeneric("updatePLC"))

setMethod(
  "updatePLC", 
  signature(x = "missing", newData = "missing"), 
  function(plc, fri, age, natDist, anthroDist, harv, resTypeLevels){
    # if(raster::res(plc)[1] != 250){
    #   stop("plc must have a resolution of 250m", call. = FALSE)
    # }
    
    # # Add plc class DTN to natDist
    DTNcode <- resTypeLevels %>%
      filter(resType == "DTN") %>%
      pull(code)
    # 
    # ATNcode <- resTypeLevels %>% 
    #   filter(resType == "ATN") %>% 
    #   pull(code)
    # 
    # natDist <- natDist == 1 
    
    # Find cells in PLC that meet condition for update: had natDist but have
    # age > 35 OR where FRI indicates Harvest (since PLC was created)
    toUpdate <- ((age > 35) & (natDist == 1)) | harv == 1 
    
    # only update if there is FRI data available
    toUpdate <- toUpdate - is.na(fri)
    
    # Update cells in toUpdate in PLC to value in FRI
    plc <- mask(plc, mask = toUpdate, maskvalue = 1, 
                updatevalue = NA)
    plc <- cover(plc, fri)
    
    # remove cells that were updated from natDist
    natDist <- mask(natDist, mask = toUpdate, maskvalue = 1, 
                    updatevalue = 0)
    
    # update PLC to DTN if natDist and age <= 35
    yfNatDist <- (age <= 35 | is.na(age)) & natDist == 1
    
    plc <- mask(plc, yfNatDist, maskvalue = 1, 
                updatevalue = DTNcode)
    
    # convert to 16 ha resolution stack of ResType proportion to match 16 ha
    # hexagons in Rempel
    tmplt <- raster(plc) %>% raster::`res<-`(c(400, 400))
    plc <- raster::layerize(plc, classes = resTypeLevels$code) %>% 
      raster::resample(tmplt, method = "bilinear")
    
    allDist16ha <- raster::stack(harv, anthroDist, 
                                 natDist) %>% 
      raster::resample(tmplt, method = "bilinear")
    
    # Get proportion of land in 16 ha area that has each type of disturbance
    watCode <- resTypeLevels %>% 
      filter(resType == "LGW") %>% 
      pull(code)
    
    # get proportion land
    land <- 1 - plc[[watCode]]
    
    # divide prop disturbance by prop land 
    #allDist16ha2 <- allDist16ha / land 
    
    propLandDist <- Vectorize(function(dist, land){
      ifelse(land < dist, 1, 
             ifelse(land == 0, 0, dist/land))
    })
    
    allDist16ha <- raster::overlay(allDist16ha, land, fun = propLandDist)
    
    
    # make plc have 0 forest classes when 16 ha area disturbed > 0.35 by natural
    # disturbance or anthropogenic disturbance
    toChange <- resTypeLevels %>% filter(!resType %in% c("DTN", "LGW")) %>% 
      pull(code)
    
    for (i in toChange) {
      plc[[i]] <- raster::mask(plc[[i]], max(allDist16ha[[2:3]] > 0.35),
                               maskvalue = 1,
                               updatevalue = 0)
    }
    
    
    # make plc have 0 forest classes when 16 ha area disturbed > 0.35 by harvest
    # disturbance but don't change wetlands
    toChange <- resTypeLevels %>% 
      filter(!resType %in% c("DTN", "LGW", "LGTP", "LGOP", "ST")) %>% 
      pull(code)
    
    for (i in toChange) {
      plc[[i]] <- raster::mask(plc[[i]], allDist16ha[[1]] > 0.35,
                               maskvalue = 1, 
                               updatevalue = 0)
    }

    return(lst(plc, natDist))
  }
)

# Not currently used idea was to avoid slow steps if possible when updating data
# but doesn't make sense with current version of the other updatePLC method

# setMethod(
#   "updatePLC", 
#   signature(x = "CaribouHabitat", newData = "list"), 
#   function(x, newData, resTypeLevels){
#     DTNcode <- resTypeLevels %>% 
#       filter(resType == "DTN") %>% 
#       pull(code)
#     
#     if(!is.null(newData$fri)){
#       message("updating plc based on fri")
#       toUpdate <- newData$fri %>% raster::setValues(0)
#       
#       if(!is.null(newData$natDist) && !is.null(newData$age)){
#         message("updating plc where natDist = 1 and age > 35")
#         # Add plc class DTN to natDist
#         newData$natDist <- newData$natDist == 1 | x@plc == DTNcode
#         
#         # Find cells in PLC that meet condition for update: had natDist but have
#         # age > 35 OR where FRI indicates Harvest since PLC was created
#         toUpdate <- (newData$age > 35) & (newData$natDist == 1) 
#       }
#       
#       if(!is.null(newData$harv)){
#         message("updating plc where harv = 1")
#         toUpdate <- toUpdate == 1 | newData$harv == 1
#       }
#       
#       # only update if there is FRI data available
#       toUpdate <- toUpdate - is.na(newData$fri)
#       
#       # Update cells in toUpdate in PLC to value in FRI
#       x@plc <- mask(x@plc, mask = toUpdate, maskvalue = 1, 
#                   updatevalue = NA)
#       x@plc <- cover(x@plc, newData$fri)
#       
#       if(!is.null(newData$natDist)){
#         # remove cells that were updated from natDist
#         newData$natDist <- mask(newData$natDist, mask = toUpdate, maskvalue = 1, 
#                                 updatevalue = 0)
#       }
#     }
#     
# 
#     
#     if(!is.null(newData$natDist) && !is.null(newData$age)){
#       message("updating plc to DTN where natDist = 1 and age <= 35")
#       # update PLC to DTN if natDist and age <= 35
#       yfNatDist <- newData$age <= 35 & newData$natDist == 1
#       
#       x@plc <- mask(x@plc, yfNatDist, maskvalue = 1, 
#                   updatevalue = DTNcode)
#       
#       if(is.null(newData$fri)){
#         newData$natDist <- mask(newData$natDist, mask = newData$age > 35, 
#                                 maskvalue = 1, 
#                                 updatevalue = 0)
#       }
#     }
#     
#     # Get proportion of land in 16 ha area that has each type of disturbance
#     water <- x@plc == resTypeLevels %>% 
#       filter(resType == "LGW") %>% 
#       pull(code)
#     
#     # Seems like aggregate/disaggregate is better than resample because more
#     # clearly takes mean of cells in area
#     # stack keeps names from list
#     allDist16ha <- raster::mask(
#       raster::stack(newData[which(names(newData) %in% 
#                                     c("natDist", "anthroDist", "harv"))]),
#       water, 
#       maskvalue = 1, 
#       updatevalue = NA
#     ) %>% 
#       raster::disaggregate(fact = 5) %>% 
#       raster::aggregate(fact = 8) %>% 
#       raster::disaggregate(fact = 8) %>% 
#       raster::aggregate(fact = 5) %>% 
#       raster::crop(newData[[which.max(names(newData) %in% 
#                                        c("natDist", "anthroDist", "harv"))]])
#     
#     # Another option, a bit faster but doen't seem to exclude the NA values in
#     # the same way
#     # system.time({
#     #   tmplt <- raster(harv) %>% raster::`res<-`(c(400, 400))
#     #   harv16ha2 <-raster::mask(harv, water, maskvalue = 1, 
#     #                            updatevalue = NA) %>%
#     #     raster::resample(tmplt) %>% 
#     #     raster::disaggregate(fact = 8) %>% 
#     #     raster::aggregate(fact = 5) %>% 
#     #     raster::crop(harv)
#     # })
#     
#     # set dist to 0 for cells that are water
#     allDist16ha <- raster::mask(allDist16ha, water, maskvalue = 1, updatevalue = 0)
#     
#     # make plc have DTN class when 16 ha area disturbed > 0.35 by natural
#     # disturbance
#     if(!is.null(newData$natDist)){
#       message("updating plc to ALLDist where natDist > 35% of 16ha area")
#       x@plc <- raster::mask(x@plc, allDist16ha[["natDist"]] > 0.35,
#                           maskvalue = 1, 
#                           updatevalue = DTNcode) 
#     }
#     
#     # set dist to 0 for cells that are already DTN
#     allDist16ha <- raster::mask(allDist16ha, x@plc,
#                                 maskvalue = DTNcode, updatevalue = 0)
#     
#     if(!is.null(newData$anthroDist) || !is.null(newData$harv)){
#       message("updating plc to ATN where harv or anthroDist > 35% of 16ha area")
#       # For anthropogenic disturbance give class ATN and don't change ST, LGTP, LGOP
#       ALLDISTcode <- resTypeLevels %>% 
#         filter(resType == "ATN") %>% 
#         pull(code)
#       
#       notChanged <- resTypeLevels %>% 
#         filter(resType %in% c("LGTP", "LGOP", "ST")) %>% 
#         pull(code)
#       
#       allDist16ha <- raster::mask(allDist16ha, 
#                                   x@plc == notChanged[1] | 
#                                     x@plc == notChanged[2] |
#                                     x@plc == notChanged[3],
#                                   maskvalue = 1, updatevalue = 0)
#       
#       # suppress warinings if one of the rasters is missing since we only want
#       # the max
#       x@plc <- raster::mask(
#         x@plc, 
#         suppressWarnings(max(allDist16ha[[which(names(allDist16ha) %in% 
#                                                   c("anthroDist", "harv"))]] > 0.35)),
#         maskvalue = 1, 
#         updatevalue = ALLDISTcode
#       )
#     }
# 
#     if(is.null(newData$natDist)){
#       newData$natDist <- x@natDist
#     }
#     
#     return(lst(plc = x@plc, natDist = newData$natDist))
#   }
# )