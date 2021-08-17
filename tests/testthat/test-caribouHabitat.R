context("Overall process to run caribouHabitat")

# create a small dataset for testing #==========================================
# Only needs to happen once
# pth_base <- "vignettes/Example_data/"
# allData <- lst(
#   landCoverP = paste0(pth_base, "plc250.tif"),
#   eskerP = paste0(pth_base, "esker.shp"),
#   natDistP = paste0(pth_base, "natDist250.tif"),
#   anthroDistP = paste0(pth_base, "anthroDist250.tif"),
#   roadsP = paste0(pth_base, "road_ORNMNRFROF2010.shp"),
#   railP = paste0(pth_base, "rail.shp"),
#   utilitiesP = paste0(pth_base, "util2010.shp")
# )
# anthroDistD <- raster(allData$anthroDistP)
# 
# # Interactively get box extent
# # plot(anthroDistD)
# # ext <- raster::drawExtent()
# 
# ext <- raster::extent(c(xmin = 682454.7, xmax = 770621,
#                         ymin = 12563591, ymax = 12640381))
# box <- st_bbox(ext) %>% st_as_sfc() %>% st_as_sf() %>% st_set_crs(st_crs(landCoverD))
# 
# # convert all to modern crs
# 
# use_crs <- st_read(allData$eskerP) %>% st_crs()
# 
# # crop each data layer with box and save in tests
# 
# cropAllData <- function(pth, savePth, cropBox, use_crs) {
#   savePth <- paste0("tests/testthat/data/", savePth)
#   if(grepl(".shp$", pth)){
#     dat <- st_read(pth)
#     cropBox <- st_transform(cropBox, st_crs(dat))
#     dat <- dat %>%
#       st_crop(cropBox) %>%
#       st_transform(use_crs)
#     st_write(dat, paste0(savePth, ".shp"), append = FALSE)
#   }
#   if(grepl(".tif$", pth)){
#     dat <- raster(pth)
#     cropBox <- st_transform(cropBox, st_crs(dat))
#     dat <- dat %>% raster::crop(cropBox) %>%
#       raster::projectRaster(crs = use_crs$wkt)
# 
#     raster::writeRaster(dat, paste0(savePth, ".tif"), overwrite = TRUE)
#   }
# }
# 
# purrr::walk2(allData, stringr::str_remove(names(allData),"P$"),
#              cropAllData, cropBox = box, use_crs = use_crs)
# 
# st_write(box %>% st_transform(use_crs),
#          "tests/testthat/data/projectPoly.shp", append = FALSE)
# 
# linFeatD <- combineLinFeat(list(roads = "tests/testthat/data/roads.shp",
#                                 rail = "tests/testthat/data/rail.shp",
#                                 utilities = "tests/testthat/data/utilities.shp"))
# st_write(linFeatD, "tests/testthat/data/linFeat.shp", append = FALSE)
# 
# linFeatDras <- rasterizeLineDensity(linFeatD, anthroDistD)
# writeRaster(linFeatDras, "tests/testthat/data/linFeatTif.tif", overwrite = TRUE)
# 
# eskerDras <- rasterizeLineDensity(st_read("tests/testthat/data/esker.shp"),
#                                   anthroDistD)
# 
# writeRaster(eskerDras, "tests/testthat/data/eskerTif.tif", overwrite = TRUE)

# Test all different ways to run from paths #===================================

# Note working directory is tests/testthat when tests are run
# to run interactively use
# pthBase <- "tests/testthat/data/"
pthBase <- "data/"

paths_eskshp_linFshp <- caribouHabitat(
  landCover = paste0(pthBase, "landCover", ".tif"),
  esker = paste0(pthBase, "esker", ".shp"),
  natDist = paste0(pthBase, "natDist", ".tif"),
  anthroDist = paste0(pthBase, "anthroDist", ".tif"),
  linFeat = paste0(pthBase, "linFeat", ".shp"),
  projectPoly = paste0(pthBase, "projectPoly", ".shp"),
  linFeatSave = paste0(pthBase, "linFeatTif400", ".tif"),
  eskerSave = paste0(pthBase, "eskerTif400", ".tif"), 
  caribouRange = "Churchill", 
  winArea = 500
)

paths_list_linF <- caribouHabitat(
  landCover = paste0(pthBase, "landCover", ".tif"),
  esker = paste0(pthBase, "esker", ".shp"),
  natDist = paste0(pthBase, "natDist", ".tif"),
  anthroDist = paste0(pthBase, "anthroDist", ".tif"),
  linFeat = list(roads = paste0(pthBase, "roads.shp"),
                 rail = paste0(pthBase, "rail.shp"),
                 utilities = paste0(pthBase, "utilities.shp")), 
  projectPoly = paste0(pthBase, "projectPoly", ".shp"), 
  caribouRange = "Churchill", 
  winArea = 500
)


test_that("results match when input is paths",{
  expect_equal(paths_eskshp_linFshp@habitatUse, paths_list_linF@habitatUse)
})


# Test all different ways to run with data #====================================
landCoverD = raster(paste0(pthBase, "landCover", ".tif")) %>% 
  reclassPLC()
eskerDras = raster(paste0(pthBase, "eskerTif400", ".tif"))
eskerDshp = st_read(paste0(pthBase, "esker", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
natDistD = raster(paste0(pthBase, "natDist", ".tif"))
anthroDistD = raster(paste0(pthBase, "anthroDist", ".tif"))
linFeatDras = raster(paste0(pthBase, "linFeatTif400", ".tif"))
projectPolyD = st_read(paste0(pthBase, "projectPoly", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDshp = st_read(paste0(pthBase, "linFeat", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")

data_esktif_linFtif <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, natDist = natDistD,
  anthroDist = anthroDistD, linFeat = linFeatDras, projectPoly = projectPolyD,
  caribouRange = "Churchill", 
  winArea = 500
)

data_eskshp_linFshp <- caribouHabitat(
  landCover = landCoverD, esker = eskerDshp, natDist = natDistD, 
  anthroDist = anthroDistD, linFeat = linFeatDshp, projectPoly = projectPolyD,
  caribouRange = "Churchill", 
  winArea = 500
)

data_list_linF <- caribouHabitat(
  landCover = landCoverD, esker = eskerDshp, natDist = natDistD, 
  anthroDist = anthroDistD, 
  linFeat = list(roads = st_read(paste0(pthBase, "roads.shp"), quiet = TRUE) %>% 
                   st_set_agr("constant"),
                 rail = st_read(paste0(pthBase, "rail.shp"), quiet = TRUE) %>% 
                   st_set_agr("constant"),
                 utilities = st_read(paste0(pthBase, "utilities.shp"), quiet = TRUE)%>% 
                   st_set_agr("constant")), 
  projectPoly = projectPolyD,
  caribouRange = "Churchill", 
  winArea = 500
)


test_that("raster road input works as expected", {
  data_list_linFrdRast1 <- caribouHabitat(
    landCover = landCoverD, esker = eskerDshp, natDist = natDistD, 
    anthroDist = anthroDistD, 
    linFeat = list(roads = linFeatDras > 10,
                   rail = st_read(paste0(pthBase, "rail.shp"), quiet = TRUE) %>% 
                     st_set_agr("constant"),
                   utilities = st_read(paste0(pthBase, "utilities.shp"), quiet = TRUE)%>% 
                     st_set_agr("constant")), 
    projectPoly = projectPolyD,
    caribouRange = "Churchill", 
    winArea = 500, 
    ptDensity = 1
  )
  data_list_linFrdRast2 <- caribouHabitat(
    landCover = landCoverD, esker = eskerDshp, natDist = natDistD, 
    anthroDist = anthroDistD, 
    linFeat = list(roads = linFeatDras > 10,
                   rail = st_read(paste0(pthBase, "rail.shp"), quiet = TRUE) %>% 
                     st_set_agr("constant"),
                   utilities = st_read(paste0(pthBase, "utilities.shp"), quiet = TRUE)%>% 
                     st_set_agr("constant")), 
    projectPoly = projectPolyD,
    caribouRange = "Churchill", 
    winArea = 500, 
    ptDensity = 2
  )
  expect_gt(data_list_linFrdRast2@linFeat %>% raster::cellStats(max),
            data_list_linFrdRast1@linFeat %>% raster::cellStats(max))
})

test_that("results match when input is paths or data, or list of either for linFeat",{
  expect_equal(data_esktif_linFtif@habitatUse, paths_eskshp_linFshp@habitatUse)
  expect_equal(data_eskshp_linFshp@habitatUse, data_esktif_linFtif@habitatUse)
  expect_equal(data_list_linF@habitatUse, data_esktif_linFtif@habitatUse)
})

test_that("results are different when disturbance is missing", {

  data_noAnthro <- caribouHabitat(
    landCover = landCoverD, 
    esker = eskerDras, 
    #anthroDist = anthroDistD,
    natDist = natDistD, 
    linFeat = linFeatDras, 
    projectPoly = projectPolyD,
    caribouRange = "Churchill", 
    winArea = 500
  )
  
  data_noNatDist <- caribouHabitat(
    landCover = landCoverD, 
    esker = eskerDras, 
    anthroDist = anthroDistD,
    #natDist = natDistD, 
    linFeat = linFeatDras, 
    projectPoly = projectPolyD,
    caribouRange = "Churchill", 
    winArea = 500
  )
  
  data_noDist <- caribouHabitat(
    landCover = landCoverD, 
    esker = eskerDras, 
    #anthroDist = anthroDistD,
    #natDist = natDistD, 
    linFeat = linFeatDras, 
    projectPoly = projectPolyD,
    caribouRange = "Churchill", 
    winArea = 500
  )
  
  expect_false(identical(data_noAnthro,  data_esktif_linFtif))
  expect_false(identical(data_noNatDist,  data_esktif_linFtif))
  expect_false(identical(data_noDist,  data_esktif_linFtif))
  
  if(interactive()){
    plot(data_esktif_linFtif, main = "Orig")
    plot(data_noAnthro, main = "No Anthro")
    plot(data_noNatDist, raster.title = "No NatDist")
    plot(data_noDist, raster.title = "No Dist")
  }

})

# Compare output raster to previous run. Use raster and not whole object in case
# we want to change the object in future
# Do to changes in CRS order we simply check that the two rasters are
# equivalent to each other
resultCompare <- readRDS(paste0(pthBase, "resultCompare.rds"))

# To update
# saveRDS(data_esktif_linFtif@habitatUse, paste0(pthBase, "resultCompare.rds"))

testthat::test_that("results match previous results",{
  testthat::expect_true(
    raster::all.equal(resultCompare, 
                      data_esktif_linFtif@habitatUse)
  )
})

