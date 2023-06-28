pth_base <- system.file("extdata", package = "caribouMetrics")

landCover <- raster(file.path(pth_base, "landCover.tif")) 

esker <- read_sf(file.path(pth_base, "esker.shp"), agr = "constant")

linFeat <- list(roads = read_sf(file.path(pth_base, 
                                          "roads.shp"), agr = "constant"),
                rail = read_sf(file.path(pth_base, "rail.shp"), agr = "constant"),
                utilities = read_sf(file.path(pth_base, "utilities.shp"), 
                                    agr = "constant"))

natDist <- sf::read_sf(file.path(pth_base, "fireAFFES2020.shp"), agr = "constant") 

anthroDist <- raster(file.path(pth_base, "anthroDist.tif"))

# make a polygon inside the bounding box of the rasters
singlePoly <- (raster::extent(landCover) - 30000) %>% st_bbox() %>%
  st_as_sfc() %>% st_as_sf() %>% st_set_crs(st_crs(landCover)) %>% 
  st_set_agr("constant")

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

tmplt <- raster(landCover) %>% raster::`res<-`(c(400, 400))

linFeat400 <- raster(file.path(pth_base, "linFeatTif.tif")) %>% 
  raster::resample(y = tmplt, method = "bilinear")

out <- loadSpatialInputs(projectPoly = singlePolyP, refRast = landCoverP, 
                         inputsList = list(esker = eskerP, 
                                           linFeat = linFeatP, 
                                           natDist = natDistP, 
                                           anthroDist = anthroDistP), 
                         reclassOptions = list(refRast = reclassPLC, 
                                               natDist = cbind(NA, 0)))

test_that("with paths", {
  expect_type(out, "list")
})

out2 <- loadSpatialInputs(
  projectPoly = singlePoly, refRast = landCover, 
  inputsList = list(esker = esker, 
                    linFeat = linFeat, 
                    natDist = natDist, 
                    anthroDist = anthroDist), 
  convertToRastDens = c("esker", "linFeat"),
  useTemplate = c("esker", "linFeat"),
  reclassOptions = list(refRast = reclassPLC, 
                        natDist = list(fn = reclassDist,
                                       endYr = 2020, 
                                       numCumYrs = 30,
                                       dateField = "FIRE_YEAR"))
)

test_that("with spatial objects", {
  expect_type(out2, "list")
})

test_that("use pre-loaded spatial inputs in caribouHabitat or disturbanceMetrics", {
  res1 <- caribouHabitat(preppedData = out4, caribouRange = "Churchill")
  
  # should be the same as if they were used directly
  res2 <- caribouHabitat(landCover = reclassPLC(landCover),
                         projectPoly = singlePoly,
                         esker = esker, 
                         linFeat = linFeat, 
                         natDist = reclassDist(natDist, endYr = 2020, 
                                               numCumYrs = 30,
                                               dateField = "FIRE_YEAR", 
                                               template = landCover), 
                         anthroDist = anthroDist,
                         caribouRange = "Churchill")
  
  expect_true(all.equal(res1@habitatUse, res2@habitatUse))
  
  dm1 <- disturbanceMetrics(preppedData = out)
  
  dm2 <- disturbanceMetrics(projectPoly = singlePolyP, landCover = landCoverP, 
                            linFeat = linFeatP, 
                            natDist = natDistP, anthroDist = anthroDistP)
  
  expect_true(all.equal(dm1@processedData, dm2@processedData))
  
})

test_that("works when inputs have different crs", {
  
  esker2 <- esker %>% st_transform(crs = 5070)
  
  out3 <- loadSpatialInputs(
    projectPoly = singlePoly, refRast = landCover, 
    inputsList = list(esker = esker2, 
                      linFeat = linFeat, 
                      natDist = natDist, 
                      anthroDist = anthroDist), 
    convertToRastDens = c("esker", "linFeat"),
    useTemplate = c("esker", "linFeat"),
    reclassOptions = list(refRast = reclassPLC, 
                          natDist = list(fn = reclassDist,
                                         endYr = 2020, 
                                         numCumYrs = 30,
                                         dateField = "FIRE_YEAR"))
  )
  expect_type(out3, "list")
})

out4 <- loadSpatialInputs(
  projectPoly = singlePoly, refRast = landCover, 
  inputsList = list(esker = esker, 
                    linFeat = linFeat, 
                    natDist = natDist, 
                    anthroDist = anthroDist), 
  convertToRastDens = c("esker", "linFeat"),
  useTemplate = c("esker", "linFeat"),
  reclassOptions = list(refRast = reclassPLC, 
                        natDist = list(fn = reclassDist,
                                       endYr = 2020, 
                                       numCumYrs = 30,
                                       dateField = "FIRE_YEAR")),
  rastOut = "raster"
)

test_that("outputs raster when asked", {
  expect_s4_class(out4$refRast, "RasterLayer")
})
