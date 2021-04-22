context("caribouHabitat messages")

# Load data

# Note working directory is tests/testthat when tests are run
# to run interactively use
# pthBase <- "tests/testthat/data/"
pthBase <- "data/"

# load data for tests
landCoverD = raster(paste0(pthBase, "plc", ".tif")) %>% 
  reclassPLC()
eskerDras = raster(paste0(pthBase, "eskerTif", ".tif"))
eskerDshp = st_read(paste0(pthBase, "esker", ".shp"), quiet = TRUE)
natDistD = raster(paste0(pthBase, "natDist", ".tif"))
anthroDistD = raster(paste0(pthBase, "anthroDist", ".tif"))
linFeatDras = raster(paste0(pthBase, "linFeatTif", ".tif"))
projectPolyD = st_read(paste0(pthBase, "projectPoly", ".shp"), quiet = TRUE)
linFeatDshp = st_read(paste0(pthBase, "linFeat", ".shp"), quiet = TRUE)

test_that("error if a required arg is missing", {
  expect_error(caribouHabitat(
    #landCover = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_error(caribouHabitat(
    landCover = paste0(pthBase, "plc", ".tif"),
    #esker = paste0(pthBase, "esker", ".shp"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    caribouRange = "Churchill"
  ), "is missing, with no default")

  expect_error(caribouHabitat(
    landCover = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    #linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_error(caribouHabitat(
    landCover = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    #projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"),  
    caribouRange = "Churchill"
  ), "is missing, with no default")
 
  expect_error(caribouHabitat(
    landCover = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"),  
    #caribouRange = "Churchill"
  ), "is missing, with no default")
})

test_that("error if a path file is not found", {
  expect_error(caribouHabitat(
    landCover = paste0(pthBase, "plctest", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"),  
    caribouRange = "Churchill"
  ), "do.es. not exist")
  expect_error(caribouHabitat(
    landCover = paste0(pthBase, "plctest", ".tif"),
    esker = paste0(pthBase, "eskertest", ".shp"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    caribouRange = "Churchill"
  ), "do.es. not exist")
})

test_that("error if some args are paths and others spatial", {
  expect_error(caribouHabitat(
    landCover = paste0(pthBase, "plc", ".tif"),
    esker = st_read(paste0(pthBase, "esker", ".shp"), quiet = TRUE),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"),  
    caribouRange = "Churchill"
  ), ".*supplied as sf or raster objects or character.*")
  expect_error(caribouHabitat(
    landCover = landCoverD,
    esker = paste0(pthBase, "eskertest", ".shp"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    caribouRange = "Churchill"
  ), ".*supplied as sf or raster objects or character.*")
})

test_that("error if file is in wrong format",{
  expect_error(caribouHabitat(
    landCover = paste0(pthBase, "linFeat", ".shp"),
    esker = paste0(pthBase, "esker", ".shp"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    caribouRange = "Churchill"
  ), "raster file must be provided for")
  expect_error(caribouHabitat(
    landCover = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "plc", ".tif"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"),  
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
    "different extent")
  
    projectPolyD2 <- mutate(projectPolyD, geometry = geometry + 20000) %>% 
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
    "landCover must have a projected CRS")
})

test_that("error if landCover is not in resource types", {
  landCoverD4 <- raster(paste0(pthBase, "plc", ".tif"))
  
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
