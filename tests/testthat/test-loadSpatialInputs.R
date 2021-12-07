pth_base <- system.file("extdata", package = "caribouMetrics")

landCover <- raster(file.path(pth_base, "landCover.tif")) 

esker <- read_sf(file.path(pth_base, "esker.shp"))

linFeat <- list(roads = read_sf(file.path(pth_base, 
                                          "roads.shp")),
                rail = read_sf(file.path(pth_base, "rail.shp")),
                utilities = read_sf(file.path(pth_base, "utilities.shp")))

natDist <- sf::st_read(file.path(pth_base, "fireAFFES2020.shp")) 

anthroDist <- raster(file.path(pth_base, "anthroDist.tif"))

# make a polygon inside the bounding box of the rasters
singlePoly <- (raster::extent(landCover) - 30000) %>% st_bbox() %>%
  st_as_sfc() %>% st_as_sf() %>% st_set_crs(st_crs(landCover))

# paths
landCoverP <- file.path(pth_base, "landCover.tif")

eskerP <- file.path(pth_base, "esker.shp")

linFeatP <- list(roads = file.path(pth_base, 
                                          "roads.shp"),
                rail = file.path(pth_base, "rail.shp"),
                utilities = file.path(pth_base, "utilities.shp"))

natDistP <- file.path(pth_base, "natDist.tif") 

anthroDistP <- file.path(pth_base, "anthroDist.tif")

singlePolyP <- file.path(pth_base, "projectPoly.shp")

tmplt <- raster(landCover) %>% `res<-`(c(400, 400))

linFeat400 <- raster(file.path(pth_base, "linFeatTif.tif")) %>% 
  resample(y = tmplt, method = "bilinear")

# with paths
out <- loadSpatialInputs(projectPoly = singlePolyP, refRast = landCoverP, 
                         inputsList = list(esker = eskerP, 
                                           linFeat = linFeatP, 
                                           natDist = natDistP, 
                                           anthroDist = anthroDistP), 
                         convertToRast = c("esker", "linFeat"),
                         reclassOptions = list(refRast = reclassPLC, 
                                               natDist = cbind(NA, 0)))

# with spatial objects
out2 <- loadSpatialInputs(
  projectPoly = singlePoly, refRast = landCover, 
  inputsList = list(esker = esker, 
                    linFeat = linFeat, 
                    natDist = natDist, 
                    anthroDist = anthroDist), 
  convertToRast = c("esker", "linFeat"),
  reclassOptions = list(refRast = reclassPLC, 
                        natDist = list(fn = reclassDist,
                                       endYr = 2020, 
                                       numCumYrs = 30,
                                       dateField = "FIRE_YEAR")),
  bufferWidth = 500
)

# try with linFeat as a raster that is already 400 res
out3 <- loadSpatialInputs(
  projectPoly = singlePoly, refRast = landCover, 
  inputsList = list(esker = esker, 
                    linFeat = linFeat400, 
                    natDist = natDist, 
                    anthroDist = anthroDist), 
  convertToRast = c("esker"),
  reclassOptions = list(refRast = reclassPLC, 
                        natDist = list(fn = reclassDist,
                                       endYr = 2020, 
                                       numCumYrs = 30,
                                       dateField = "FIRE_YEAR")),
  altTemplate = tmplt,
  useTemplate = c("linFeat"),
  bufferWidth = 500
)
