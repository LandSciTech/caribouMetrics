context("test disturbanceMetrics")

# Note working directory is tests/testthat when tests are run
# to run interactively use
# pthBase <- "tests/testthat/data/"
pthBase <- "data/"

# load example data
plcD = raster(paste0(pthBase, "landCover", ".tif")) # Defines the study area - NA values are omitted from calculation, everything else is included.
natDistD = raster(paste0(pthBase, "natDist", ".tif"))
anthroDistD = raster(paste0(pthBase, "anthroDist", ".tif"))
projectPolyD = st_read(paste0(pthBase, "projectPoly", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDshp = st_read(paste0(pthBase, "linFeat", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
roadsD = st_read(paste0(pthBase, "roads.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
railD = st_read(paste0(pthBase, "rail.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
utilitiesD = st_read(paste0(pthBase, "utilities.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDras = raster(paste0(pthBase, "linFeatTif", ".tif"))


dm <- disturbanceMetrics(
  landCover = plcD,
  natDist = natDistD, 
  anthroDist = anthroDistD, 
  linFeat = linFeatDshp, 
  projectPoly = projectPolyD,
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)

dm@disturbanceMetrics
plot(dm@processedData)

dm_path <- disturbanceMetrics(
  landCover = paste0(pthBase, "landCover", ".tif"),
  natDist = paste0(pthBase, "natDist", ".tif"), 
  anthroDist = paste0(pthBase, "anthroDist", ".tif"), 
  linFeat = paste0(pthBase, "linFeat", ".shp"), 
  projectPoly = paste0(pthBase, "projectPoly", ".shp"),
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)
dm_lflist <- disturbanceMetrics(
  landCover = plcD,
  natDist = natDistD, 
  anthroDist = anthroDistD, 
  linFeat = list(roads = st_read(paste0(pthBase, "roads.shp"), quiet = TRUE) %>% 
                   st_set_agr("constant"),
                 rail = st_read(paste0(pthBase, "rail.shp"), quiet = TRUE) %>% 
                   st_set_agr("constant"),
                 utilities = st_read(paste0(pthBase, "utilities.shp"), quiet = TRUE)%>% 
                   st_set_agr("constant")),
  projectPoly = projectPolyD,
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)


test_that("results match when input is paths or data or list for linFeat",{
  expect_equal(dm@disturbanceMetrics, dm_path@disturbanceMetrics)
  expect_equal(dm@disturbanceMetrics, dm_lflist@disturbanceMetrics)
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
               tolerance = 0.01)
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
# things to test
# linFeat sf with some points and lines?
# anthro or nat Dist missing

# Compare output to previous run. This will raise a flag if the result has
# changed. Update the stored result if the change was expected.
resultCompare <- readRDS(paste0(pthBase, "dm_resultCompare.rds"))

# To update
# saveRDS(dm@disturbanceMetrics, paste0(pthBase, "dm_resultCompare.rds"))

testthat::test_that("results match previous results",{
  testthat::expect_equal(dm@disturbanceMetrics, resultCompare)
})
