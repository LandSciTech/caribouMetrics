
pthBase <- system.file("extdata", package = "caribouMetrics")

# this makes sure raster package is attached (not loaded) so that as(x, "Raster") works
raster::raster()

# load example data
plcD = terra::rast(file.path(pthBase, "landCover.tif")) # Defines the study area - NA values are omitted from calculation, everything else is included.
natDistD = terra::rast(file.path(pthBase, "natDist.tif"))
anthroDistD = terra::rast(file.path(pthBase, "anthroDist.tif"))
projectPolyD = st_read(file.path(pthBase, "projectPoly.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDshp = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
roadsD = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
railD = st_read(file.path(pthBase, "rail.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
utilitiesD = st_read(file.path(pthBase, "utilities.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDras = terra::rast(file.path(pthBase, "linFeatTif.tif"))

dm <- disturbanceMetrics(
  landCover = plcD,
  natDist = natDistD, 
  anthroDist = anthroDistD, 
  linFeat = linFeatDshp, 
  projectPoly = projectPolyD,
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)

test_that("results not na",{
  dm@disturbanceMetrics %>% unlist() %>% is.na() %>% any() %>% `!`() %>% 
    expect_true()
})

# dm@disturbanceMetrics
# plot(dm@processedData)

dm_path <- disturbanceMetrics(
  landCover = file.path(pthBase, "landCover.tif"),
  natDist = file.path(pthBase, "natDist.tif"), 
  anthroDist = file.path(pthBase, "anthroDist.tif"), 
  linFeat = file.path(pthBase, "roads.shp"), 
  projectPoly = file.path(pthBase, "projectPoly.shp"),
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)

dm_lflist <- disturbanceMetrics(
  landCover = plcD,
  natDist = natDistD, 
  anthroDist = anthroDistD, 
  linFeat = list(roads = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
                   st_set_agr("constant"),
                 rail = st_read(file.path(pthBase, "rail.shp"), quiet = TRUE) %>% 
                   st_set_agr("constant"),
                 utilities = st_read(file.path(pthBase, "utilities.shp"), quiet = TRUE)%>% 
                   st_set_agr("constant")),
  projectPoly = projectPolyD,
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)


test_that("results match when input is paths or data",{
  expect_equal(dm@disturbanceMetrics, dm_path@disturbanceMetrics)
})

test_that("error if rasters don't align",{
  expect_error(disturbanceMetrics(
    landCover = plcD,
    natDist = natDistD, 
    anthroDist = anthroDistD, 
    linFeat = terra::`res<-`(linFeatDras, 400), 
    projectPoly = projectPolyD
  ), "rasters do not have the same")
})

dm_lfrast <- disturbanceMetrics(
  landCover = plcD,
  natDist = natDistD, 
  anthroDist = anthroDistD, 
  linFeat = linFeatDras, 
  projectPoly = projectPolyD,
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)

test_that("results are similar with rast vs sf linFeat",{
  expect_equal(dm@disturbanceMetrics, dm_lfrast@disturbanceMetrics,
               tolerance = 5)
})

dm_noAnthro <- disturbanceMetrics(
  landCover = plcD,
  natDist = natDistD, 
  #anthroDist = anthroDistD, 
  linFeat = linFeatDras, 
  projectPoly = projectPolyD,
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)

dm_noNat <- disturbanceMetrics(
  landCover = plcD,
  #natDist = natDistD, 
  anthroDist = anthroDistD, 
  linFeat = linFeatDras, 
  projectPoly = projectPolyD,
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)

dm_noDist <- disturbanceMetrics(
  landCover = plcD,
  #natDist = natDistD, 
  #anthroDist = anthroDistD, 
  linFeat = linFeatDras, 
  projectPoly = projectPolyD,
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)

test_that("results are different without disturbances",{
  expect_false(identical(dm_noAnthro,  dm))
  expect_false(identical(dm_noNat,  dm))
  expect_false(identical(dm_noDist,  dm))
})

test_that("fire_excl_anthro lt fire",
          expect_lt(dm@disturbanceMetrics$fire_excl_anthro, 
                    dm@disturbanceMetrics$Fire)
)

test_that("line sf with points works", {
  # make points
  pts <- sf::st_sf(linFID = 1:10, geometry = sf::st_sample(projectPolyD, 10))
  
  linFeatDpts <- linFeatDshp %>% bind_rows(pts) %>% st_set_agr("constant")
  
  dm_pts <- disturbanceMetrics(
    landCover = plcD,
    natDist = natDistD, 
    anthroDist = anthroDistD, 
    linFeat = linFeatDpts, 
    projectPoly = projectPolyD,
    padFocal = FALSE, # assume data outside area is 0 for all variables
    bufferWidth = 500
  )
  
  expect_gt(dm_pts@disturbanceMetrics$Anthro, 
            dm@disturbanceMetrics$Anthro)
  
})

test_that("linBuffMethod sf works",{
  dm_sf <- disturbanceMetrics(
    landCover = plcD,
    natDist = natDistD, 
    anthroDist = anthroDistD, 
    linFeat = linFeatDshp, 
    projectPoly = projectPolyD,
    bufferWidth = 500,
    linBuffMethod = "sf"
  )
  
  expect_s4_class(dm_sf, "DisturbanceMetrics")
})

test_that("RasterLayer input works",{
  dm_rast <- disturbanceMetrics(
    landCover = as(plcD, "Raster"),
    natDist = as(natDistD, "Raster"), 
    anthroDist = as(anthroDistD, "Raster"), 
    linFeat = linFeatDshp, 
    projectPoly = projectPolyD,
    bufferWidth = 500,
    linBuffMethod = "sf"
  )
  
  expect_s4_class(dm_rast, "DisturbanceMetrics")
})

test_that("NAs handled correctly", {
  natDistD[natDistD == 0] <- NA
  anthroDistD[anthroDistD == 0] <- NA
  
  dm_na <- disturbanceMetrics(
    landCover = plcD,
    natDist = natDistD, 
    anthroDist = anthroDistD, 
    linFeat = linFeatDshp, 
    projectPoly = projectPolyD,
    bufferWidth = 500
  )
  
  # Nas in disturbance should have no impact
  expect_equal(dm@disturbanceMetrics, dm_na@disturbanceMetrics)
  
  plcD[plcD == 3] <- NA
  
  dm_na_lc <- disturbanceMetrics(
    landCover = plcD,
    natDist = natDistD, 
    anthroDist = anthroDistD, 
    linFeat = linFeatDshp, 
    projectPoly = projectPolyD,
    bufferWidth = 500
  )
  ext <- terra::ext(749233.8, 757072.9, 12564948, 12572160)
  na_out <- terra::crop(is.na(dm_na_lc@processedData[[1]]), ext)
  na_lc <- terra::crop(is.na(dm_na_lc@landCover), ext)
  expect_true(terra::global(terra::compare(na_out, na_lc, oper = "=="), min)[1,1] == 1)
})

# Compare output to previous run. This will raise a flag if the result has
# changed. Update the stored result if the change was expected.
resultCompare <- readRDS(file.path(test_path("data"), "dm_resultCompare.rds"))

# To update
# saveRDS(dm@disturbanceMetrics, file.path("tests/testthat/data", "dm_resultCompare.rds"),
# version = 2)

testthat::test_that("results match previous results",{
  testthat::expect_equal(dm@disturbanceMetrics, resultCompare)
})
