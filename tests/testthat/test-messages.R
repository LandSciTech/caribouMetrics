context("caribouHabitat messages")

# Load data

# Note working directory is tests/testthat when tests are run
# to run interactively use
# pthBase <- "tests/testthat/data/"
pthBase <- system.file("extdata", package = "caribouMetrics")

# load data for tests
landCoverD = raster(file.path(pthBase, "landCover.tif")) %>% 
  reclassPLC()
eskerDras = raster(file.path(pthBase, "eskerTif.tif"))
eskerDshp = st_read(file.path(pthBase, "esker.shp"), quiet = TRUE)
natDistD = raster(file.path(pthBase, "natDist.tif"))
anthroDistD = raster(file.path(pthBase, "anthroDist.tif"))
linFeatDras = raster(file.path(pthBase, "linFeatTif.tif"))
projectPolyD = st_read(file.path(pthBase, "projectPoly.shp"), quiet = TRUE)
linFeatDshp = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE)

test_that("error if a required arg is missing", {
  expect_error(caribouHabitat(
    #landCover = file.path(pthBase, "landCover.tif"),
    esker = file.path(pthBase, "esker.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    linFeat = file.path(pthBase, "roads.shp"),
    projectPoly = file.path(pthBase, "projectPoly.shp"),
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_error(caribouHabitat(
    landCover = file.path(pthBase, "landCover.tif"),
    #esker = file.path(pthBase, "esker.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    linFeat = file.path(pthBase, "roads.shp"),
    projectPoly = file.path(pthBase, "projectPoly.shp"),
    caribouRange = "Churchill"
  ), "is missing, with no default")

  expect_error(caribouHabitat(
    landCover = file.path(pthBase, "landCover.tif"),
    esker = file.path(pthBase, "esker.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    #linFeat = file.path(pthBase, "roads.shp"),
    projectPoly = file.path(pthBase, "projectPoly.shp"),
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_error(caribouHabitat(
    landCover = file.path(pthBase, "landCover.tif"),
    esker = file.path(pthBase, "esker.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    linFeat = file.path(pthBase, "roads.shp"),
    #projectPoly = file.path(pthBase, "projectPoly.shp"),
    caribouRange = "Churchill"
  ), "is missing, with no default")
 
  expect_error(caribouHabitat(
    landCover = file.path(pthBase, "landCover.tif"),
    esker = file.path(pthBase, "esker.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    linFeat = file.path(pthBase, "roads.shp"),
    projectPoly = file.path(pthBase, "projectPoly.shp"),  
    #caribouRange = "Churchill"
  ), "is missing, with no default")
})

test_that("error if a path file is not found", {
  expect_error(caribouHabitat(
    landCover = file.path(pthBase, "landCovertest.tif"),
    esker = file.path(pthBase, "esker.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    linFeat = file.path(pthBase, "roads.shp"),
    projectPoly = file.path(pthBase, "projectPoly.shp"), 
    caribouRange = "Churchill"
  ), "do.es. not exist")
  expect_error(caribouHabitat(
    landCover = file.path(pthBase, "landCovertest.tif"),
    esker = file.path(pthBase, "eskertest.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    linFeat = file.path(pthBase, "roads.shp"),
    projectPoly = file.path(pthBase, "projectPoly.shp"),
    caribouRange = "Churchill"
  ), "do.es. not exist")
})

test_that("error if file is in wrong format",{
  expect_error(caribouHabitat(
    landCover = file.path(pthBase, "roads.shp"),
    esker = file.path(pthBase, "esker.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    linFeat = file.path(pthBase, "roads.shp"),
    projectPoly = file.path(pthBase, "projectPoly.shp"),
    caribouRange = "Churchill"
  ), "raster file must be provided for")
  expect_error(caribouHabitat(
    landCover = file.path(pthBase, "landCover.tif"),
    esker = file.path(pthBase, "esker.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    linFeat = file.path(pthBase, "roads.shp"),
    projectPoly = file.path(pthBase, "landCover.tif"),
    caribouRange = "Churchill"
  ), "shapefile must be provided for")
})

test_that("error if caribouRange missing or doesn't match", {
  expect_error(caribouHabitat(
    landCover = landCoverD, esker = eskerDras, natDist = natDistD, 
    anthroDist = anthroDistD,
    linFeat = linFeatDras, projectPoly = projectPolyD,
    winArea = 500),
    "is missing, with no default")
  expect_error(caribouHabitat(
    landCover = landCoverD, 
    esker = eskerDras, 
    natDist = natDistD, 
    anthroDist = anthroDistD, 
    linFeat = linFeatDras, 
    projectPoly = projectPolyD,
    caribouRange = "test",
    winArea = 500),
    "Caribou Range must match one of")
}) 

test_that("error if spatial data doesn't align", {
  landCoverD2 <- raster::shift(landCoverD, 100)

  expect_error(caribouHabitat(
    landCover = landCoverD2, 
    esker = eskerDras, 
    natDist = natDistD, 
    anthroDist = anthroDistD,
    linFeat = linFeatDras, 
    projectPoly = projectPolyD, 
    caribouRange = "Churchill"),
    "rasters do not have the same")
  
    projectPolyD2 <- mutate(projectPolyD, geometry = geometry + 200000) %>% 
    st_set_crs(st_crs(projectPolyD))
  
  expect_error(caribouHabitat(
    landCover = landCoverD, 
    esker = eskerDras, 
    natDist = natDistD,  
    anthroDist = anthroDistD, 
    linFeat = linFeatDras, 
    projectPoly = projectPolyD2, 
    caribouRange = "Churchill"),
    "landCover does not overlap with")
}) 

test_that("error if landCover is in lonlat", {
  landCoverD3 <- landCoverD %>% projectRaster(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                                  method = "ngb")
  
  expect_error(caribouHabitat(
    landCover = landCoverD3, 
    esker = eskerDras, 
    natDist = natDistD,  
    anthroDist = anthroDistD, 
    linFeat = linFeatDras, 
    projectPoly = projectPolyD,
    caribouRange = "Churchill"),
    "all raster data sets must have matching CRS")
})

test_that("error if landCover is not in resource types", {
  landCoverD4 <- raster(file.path(pthBase, "landCover.tif"))
  
  expect_error(caribouHabitat(
    landCover = landCoverD4,
    esker = eskerDras, 
    natDist = natDistD,  
    anthroDist = anthroDistD, 
    linFeat = linFeatDras, 
    projectPoly = projectPolyD,
    caribouRange = "Churchill"), 
    "resource types")

})
