pthBase <- system.file("extdata", package = "caribouMetrics")

# load example data
plcD = raster(file.path(pthBase, "landCover.tif")) # Defines the study area - NA values are omitted from calculation, everything else is included.
natDistD = raster(file.path(pthBase, "natDist.tif"))
anthroDistD = raster(file.path(pthBase, "anthroDist.tif"))
projectPolyD = st_read(file.path(pthBase, "projectPoly.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDshp = st_read(file.path(pthBase, "linFeat.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
roadsD = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
railD = st_read(file.path(pthBase, "rail.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
utilitiesD = st_read(file.path(pthBase, "utilities.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDras = raster(file.path(pthBase, "linFeatTif.tif"))

# newData versions to check updating. 
ext <- raster::extent(natDistD)- 15000

msk <- raster::crop(natDistD, ext) %>% raster::setValues(1) %>% 
  raster::extend(natDistD, value = 0)

natDistD2 <- raster::mask(natDistD, msk, maskvalue = 1, updatevalue = 1)

linFeatDras2 <- raster::mask(linFeatDras, 
                             raster::crop(linFeatDras, ext) %>% 
                               raster::setValues(1) %>% 
                               raster::extend(linFeatDras, value = 0), 
                             maskvalue = 1, updatevalue = 0)

anthroDistD2 <- raster::mask(anthroDistD, msk, maskvalue = 1, updatevalue = 0)

dm <- disturbanceMetrics(
  landCover = plcD,
  natDist = natDistD, 
  anthroDist = anthroDistD, 
  linFeat = linFeatDshp, 
  projectPoly = projectPolyD,
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)

test_that("updateDisturbance works", {
  expect_message(
    dm_upd_lf <- updateDisturbance(dm, newData = list(linFeat = linFeatDras2)),
    "buffering"
  )
  
  expect_gt(results(dm)$Anthro, results(dm_upd_lf)$Anthro)
  
  expect_message(
    dm_upd_an <- updateDisturbance(dm, newData = list(anthroDist = anthroDistD2)),
    "buffering"
  )
  
  expect_gt(results(dm)$Anthro, results(dm_upd_an)$Anthro)
  

  dm_upd_nd <- updateDisturbance(dm, newData = list(natDist = natDistD2))
  
  expect_lt(results(dm)$Fire, results(dm_upd_nd)$Fire)
  
})




