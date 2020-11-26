context("Overall process to run caribouHabitat")

# create a small dataset for testing #==========================================
# Only needs to happen once
# 
# allData <- lst(
#   plcP = "./inputNV/intermediaryData/plc_aligned.tif",
#   eskerP = "inputNV/HornsethRempelInfo/Eskers_Ontario/Eskers_Ontario.shp",
#   friP = "inputNV/intermediaryData/fri.tif",
#   ageP = "inputNV/intermediaryData/age.tif",
#   natDistP = "./inputNV/intermediaryData/natDistOFRI.tif",
#   anthroDistP = "./inputNV/intermediaryData/dummyAnthroDist.tif",
#   harvP = "inputNV/intermediaryData/harvPost2000OFRI.tif",
#   linFeatP = "./inputNV/intermediaryData/linFeatVect.shp",
#   projectPolyP = "./inputNV/intermediaryData/projectPoly.shp",
#   linFeatTifP = "inputNV/intermediaryData/roadRastDen.tif",
#   eskerTifP = "inputNV/intermediaryData/eskerRastDen.tif",
#   roadsP = "inputNV/intermediaryData/roads.shp",
#   railP = "inputNV/intermediaryData/rail.shp",
#   utilitiesP = "inputNV/intermediaryData/utilities.shp"
# )
# plcD <- raster(allData$plcP)
# 
# # Interactively get box extent
# # plot(plcD)
# # ext <- raster::drawExtent()
# 
# ext <- raster::extent(c(xmin = 371624.6, xmax = 391476.2,
#                         ymin = 12719130, ymax = 12735164))
# box <- st_bbox(ext) %>% st_as_sfc() %>% st_as_sf() %>% st_set_crs(st_crs(plcD))
# 
# # crop each data layer with box and save in tests
# 
# cropAllData <- function(pth, savePth, cropBox) {
#   savePth <- paste0("tests/testthat/data/", savePth)
#   if(grepl(".shp$", pth)){
#     dat <- st_read(pth)
#     cropBox <- st_transform(cropBox, st_crs(dat))
#     dat <- dat %>% 
#       st_crop(cropBox)
#     st_write(dat, paste0(savePth, ".shp"), append = FALSE)
#   }
#   if(grepl(".tif$", pth)){
#     dat <- raster(pth) 
#     cropBox <- st_transform(cropBox, st_crs(dat))
#     dat <- dat %>% raster::crop(cropBox)
#     raster::writeRaster(dat, paste0(savePth, ".tif"), overwrite = TRUE)
#   }
# }
# 
# purrr::walk2(allData, stringr::str_remove(names(allData),"P$"),
#              cropAllData, cropBox = box)
# 
# sSimData <- getSyncSimData(ssimLib = "./inputNV/Libraries/ChurchillBC",
#                            ssimProject = "2008/09 Plans",
#                            rfuLUPath = "inputTables/OLTBorealCaribouFURegionalForestUnitMapping.csv",
#                            scnNum = 6)
# 
# friLU <- sSimData$friLU %>% filter(!is.na(RFU))
# 
# write.csv(friLU, "tests/testthat/data/friLU.csv", row.names = FALSE)

# Test all different ways to run from paths #===================================

# Note working directory is tests/testthat when tests are run
# to run interactively use
# pthBase <- "tests/testthat/data/"
pthBase <- "data/"

paths_eskshp_linFshp <- caribouHabitat(
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
  caribouRange = "Churchill", 
  winArea = 500
)

paths_list_linF <- caribouHabitat(
  plc = paste0(pthBase, "plc", ".tif"),
  esker = paste0(pthBase, "esker", ".shp"),
  fri = paste0(pthBase, "fri", ".tif"),
  age = paste0(pthBase, "age", ".tif"),
  natDist = paste0(pthBase, "natDist", ".tif"),
  anthroDist = paste0(pthBase, "anthroDist", ".tif"),
  harv = paste0(pthBase, "harv", ".tif"),
  linFeat = list(roads = paste0(pthBase, "roads.shp"),
                 rail = paste0(pthBase, "rail.shp"),
                 utilities = paste0(pthBase, "utilities.shp")), 
  projectPoly = paste0(pthBase, "projectPoly", ".shp"),
  linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
  eskerSave = paste0(pthBase, "eskerTif", ".tif"), 
  friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
  caribouRange = "Churchill", 
  winArea = 500
)


test_that("results match when input is paths",{
  expect_equal(paths_eskshp_linFshp@habitatUse, paths_list_linF@habitatUse)
})


# Test all different ways to run with data #====================================
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
linFeatDshp = st_read(paste0(pthBase, "linFeat", ".shp"), quiet = TRUE)

data_esktif_linFtif <- caribouHabitat(
  plc = plcD, esker = eskerDras, fri = friD, age = ageD, natDist = natDistD,
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = projectPolyD,
  friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
  caribouRange = "Churchill", 
  winArea = 500
)

data_eskshp_linFshp <- caribouHabitat(
  plc = plcD, esker = eskerDshp, fri = friD, age = ageD, natDist = natDistD, 
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDshp, projectPoly = projectPolyD,
  eskerSave = paste0(pthBase, "eskerTif", ".tif"),
  linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
  friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
  caribouRange = "Churchill", 
  winArea = 500
)

data_list_linF <- caribouHabitat(
  plc = plcD, esker = eskerDshp, fri = friD, age = ageD, natDist = natDistD, 
  anthroDist = anthroDistD, harv = harvD,
  linFeat = list(roads = st_read(paste0(pthBase, "roads.shp"), quiet = TRUE),
                 rail = st_read(paste0(pthBase, "rail.shp"), quiet = TRUE),
                 utilities = st_read(paste0(pthBase, "utilities.shp"), quiet = TRUE)), 
  projectPoly = projectPolyD,
  eskerSave = paste0(pthBase, "eskerTif", ".tif"),
  linFeatSave = paste0(pthBase, "linFeatTif", ".tif"),
  friLU = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE), 
  caribouRange = "Churchill", 
  winArea = 500
)

test_that("results match when input is paths or data, or list of either for linFeat",{
  expect_equal(data_esktif_linFtif@habitatUse, paths_eskshp_linFshp@habitatUse)
  expect_equal(data_eskshp_linFshp@habitatUse, data_esktif_linFtif@habitatUse)
  expect_equal(data_list_linF@habitatUse, data_esktif_linFtif@habitatUse)
})


# Compare output table to previous run. Use table and not whole object in case
# we want to change the object in future
# test_that("results match previous results", {
#   expect_known_value(data_esktif_linFtif@habitatUse, 
#                       file = paste0(pthBase, "resultCompare.rds"), 
#                       update = FALSE)
# })




