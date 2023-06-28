# pthBase <- system.file("extdata", package = "caribouMetrics")
# 
# test_that("works", {
#   oldIn <- caribouMetrics:::inputData(landCover = file.path(pthBase, "landCover.tif"),
#             esker = file.path(pthBase, "esker.shp"),
#             natDist = file.path(pthBase, "natDist.tif"),
#             anthroDist = file.path(pthBase, "anthroDist.tif"),
#             linFeat = file.path(pthBase, "roads.shp"),
#             projectPoly = file.path(pthBase, "projectPoly.shp"),
#             linFeatSave = file.path(pthBase, "linFeatTif400.tif"),
#             eskerSave = file.path(pthBase, "eskerTif400.tif"), 
#             caribouRange = data.frame(Range = "Churchill"), 
#             coefTable = coefTableHR,
#             winArea = 500)
#   newIn <- inputData(landCover = file.path(pthBase, "landCover.tif"),
#                      esker = file.path(pthBase, "esker.shp"),
#                      natDist = file.path(pthBase, "natDist.tif"),
#                      anthroDist = file.path(pthBase, "anthroDist.tif"),
#                      linFeat = file.path(pthBase, "roads.shp"),
#                      projectPoly = file.path(pthBase, "projectPoly.shp"),
#                      linFeatSave = file.path(pthBase, "linFeatTif400.tif"),
#                      eskerSave = file.path(pthBase, "eskerTif400.tif"), 
#                      caribouRange = data.frame(Range = "Churchill"), 
#                      coefTable = coefTableHR,
#                      winArea = 500)
# })
# 
# toTerra <- function(x) terra::rast(x)
# toRast <- function(x) as(x, "Raster")
# 
# terraEqual <- function(x, y){
#   if(is(x, "Raster")){
#     x <- terra::rast(x)
#   }
#   if(is(y, "Raster")){
#     y <- terra::rast(y)
#   }
#   terra::all.equal(x, y)
# }
# 
# plotDif <- function(x, y){
#   if(is(x, "SpatRaster")){
#     x <- toRast(x)
#   }
#   if(is(y, "SpatRaster")){
#     y <- toRast(y)
#   }
#   
#   raster::plot(x - y)
# }
# 
# terraEqual(newIn@landCover, oldIn@landCover)
# terraEqual(newIn@linFeat, oldIn@linFeat)
# plotDif(newIn@linFeat, oldIn@linFeat)
# terraEqual(newIn@esker, oldIn@esker)
# plotDif(newIn@esker, oldIn@esker)
# terraEqual(newIn@natDist, oldIn@natDist)
# terraEqual(newIn@anthroDist, oldIn@anthroDist)
# 
# inData <- oldIn
# oldAppD <- caribouMetrics:::applyDist(inData@landCover, inData@natDist, inData@anthroDist, 
#                      inData@attributes$tmplt)
# 
# inData <- newIn
# newAppD <- applyDist(inData@landCover, inData@natDist, inData@anthroDist, 
#                                       inData@attributes$tmplt)
# 
# plotDif(newAppD, oldAppD)
# 
# # is the difference caused by differences in bilinear interpolation?
# # seems like it is See issue:https://github.com/rspatial/terra/issues/1208
# lc_terra <- terra::segregate(newIn@landCover, classes = resTypeCode$code) %>% 
#   terra::resample(newIn@attributes$tmplt, method = "bilinear")
# 
# lc_rast <- raster::layerize(oldIn@landCover, classes = resTypeCode$code) %>% 
#   raster::resample(toRast(newIn@attributes$tmplt), method = "bilinear")
# 
# plotDif(lc_terra, lc_rast)
