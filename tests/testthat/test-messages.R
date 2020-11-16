context("caribouHabitat messages")

# Load data

# Note working directory is tests/testthat when tests are run
# to run interactively use
# pthBase <- "tests/testthat/data/"
pthBase <- "data/"

# load data for tests
plcD = raster(paste0(pthBase, "plc", ".tif"))
eskerDras = raster(paste0(pthBase, "eskerTif", ".tif"))
eskerDshp = st_read(paste0(pthBase, "esker", ".shp"), quiet = TRUE)
friD = raster(paste0(pthBase, "fri", ".tif"))
ageD = raster(paste0(pthBase, "age", ".tif"))
natDistD = raster(paste0(pthBase, "natDist", ".tif"))
anthroDistD = raster(paste0(pthBase, "anthroDist", ".tif"))
harvD = raster(paste0(pthBase, "harv", ".tif"))
linFeatDras = raster(paste0(pthBase, "linFeatTif", ".tif"))
projectPolyD = st_read(paste0(pthBase, "projectPoly", ".shp"), quiet = TRUE)
hexgridD = st_read(paste0(pthBase, "hexgrid", ".shp"), quiet = TRUE)
linFeatDshp = st_read(paste0(pthBase, "linFeat", ".shp"), quiet = TRUE)

test_that("error if a required arg is missing", {
  expect_error(caribouHabitat(
    #plc = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "plc", ".tif"),
    #esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    #fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    #age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    #natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    #linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    #projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "is missing, with no default")
  expect_warning(caribouHabitat(
    plc = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    #friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "friLU is required to process the data")
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    #caribouRange = "Churchill"
  ), "is missing, with no default")
})

test_that("error if a path file is not found", {
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "plctest", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "do.es. not exist")
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "plctest", ".tif"),
    esker = paste0(pthBase, "eskertest", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "do.es. not exist")
})

test_that("error if file is in wrong format",{
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "linFeat", ".shp"),
    esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "projectPoly", ".shp"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "raster file must be provided for")
  expect_error(caribouHabitat(
    plc = paste0(pthBase, "plc", ".tif"),
    esker = paste0(pthBase, "esker", ".shp"),
    fri = paste0(pthBase, "fri", ".tif"),
    age = paste0(pthBase, "age", ".tif"),
    natDist = paste0(pthBase, "natDist", ".tif"),
    anthroDist = paste0(pthBase, "anthroDist", ".tif"),
    harv = paste0(pthBase, "harv", ".tif"),
    linFeat = paste0(pthBase, "linFeat", ".shp"),
    projectPoly = paste0(pthBase, "age", ".tif"),
    linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
    eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"
  ), "shapefile must be provided for")
})

test_that("error if friLU is missing values in fri raster", {
  friLUD <- read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE)
  friLUD <- filter(friLUD, ID != 33)
  
  expect_error(caribouHabitat(
    plc = plcD, esker = eskerDras, fri = friD, age = ageD, natDist = natDistD, 
    anthroDist = anthroDistD, harv = harvD,
    linFeat = linFeatDras, projectPoly = projectPolyD,
    hexgridSave = paste0(pthBase, "hexgrid.shp"),
    friLU = friLUD, winArea = 500,
    caribouRange = "Churchill"),
    "All unique values in fri must be present in friLU")
}) 

test_that("error if caribouRange missing or doesn't match", {
  expect_error(caribouHabitat(
    plc = plcD, esker = eskerDras, fri = friD, age = ageD, natDist = natDistD, 
    anthroDist = anthroDistD, harv = harvD,
    linFeat = linFeatDras, projectPoly = projectPolyD,
    hexgridSave = paste0(pthBase, "hexgrid.shp"),
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE),
    winArea = 500),
    "is missing, with no default")
  expect_error(caribouHabitat(
    plc = plcD, esker = eskerDras, fri = friD, age = ageD, natDist = natDistD, 
    anthroDist = anthroDistD, harv = harvD,
    linFeat = linFeatDras, projectPoly = projectPolyD,
    caribouRange = "test",
    hexgridSave = paste0(pthBase, "hexgrid.shp"),
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE),
    winArea = 500),
    "caribouRange must match one of")
}) 

test_that("error if friLU doesn't match rfuToResType", {
  friLUD <- read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE)
  friLUD <- mutate(friLUD, RFU = ifelse(RFU == "BFDOM", "test",
                                        ifelse(RFU == "BWDOM", "test2", RFU)))
  
  expect_error(caribouHabitat(
    plc = plcD, esker = eskerDras, fri = friD, age = ageD, natDist = natDistD, 
    anthroDist = anthroDistD, harv = harvD,
    linFeat = linFeatDras, projectPoly = projectPolyD,
    hexgridSave = paste0(pthBase, "hexgrid.shp"),
    winArea = 500,
    friLU = friLUD, 
    caribouRange = "Churchill"),
    "The second column of friLU must match a regional forest unit")
}) 

test_that("error if spatial data doesn't align", {
  plcD2 <- raster::shift(plcD, 100)

  expect_error(caribouHabitat(
    plc = plcD2, esker = eskerDras, fri = friD, age = ageD, natDist = natDistD, 
    anthroDist = anthroDistD, harv = harvD,
    linFeat = linFeatDras, projectPoly = projectPolyD,
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"),
    "different extent")
  
    projectPolyD2 <- mutate(projectPolyD, geometry = geometry + 20000) %>% 
    st_set_crs(st_crs(projectPolyD))
  
  expect_error(caribouHabitat(
    plc = plcD, esker = eskerDras, fri = friD, age = ageD, natDist = natDistD,  
    anthroDist = anthroDistD, harv = harvD,
    linFeat = linFeatDras, projectPoly = projectPolyD2,
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"),
    "plc does not overlap with")
}) 

test_that("error if plc is in lonlat", {
  plcD3 <- plcD %>% projectRaster(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                                  method = "ngb")
  
  expect_error(caribouHabitat(
    plc = plcD3, esker = eskerDras, fri = friD, age = ageD, natDist = natDistD,  
    anthroDist = anthroDistD, harv = harvD,
    linFeat = linFeatDras, projectPoly = projectPolyD,
    hexgrid = hexgridD,
    friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
    caribouRange = "Churchill"),
    "plc must have a projected CRS")
})
