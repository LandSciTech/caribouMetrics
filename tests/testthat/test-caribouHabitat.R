
# create a small dataset for testing #==========================================
# Only needs to happen once
# pth_base <- "vignettes/Example_data/"
# allData <- lst(
#   landCoverP = file.path(pth_base, "plc250.tif"),
#   eskerP = file.path(pth_base, "esker.shp"),
#   natDistP = file.path(pth_base, "natDist250.tif"),
#   anthroDistP = file.path(pth_base, "anthroDist250.tif"),
#   roadsP = file.path(pth_base, "road_ORNMNRFROF2010.shp"),
#   railP = file.path(pth_base, "rail.shp"),
#   utilitiesP = file.path(pth_base, "util2010.shp")
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
# #st_write(linFeatD, "tests/testthat/data/linFeat.shp", append = FALSE)
#
# make roads smaller file size by simplifying
# linFeatDshp = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
#   st_set_agr("constant")
# st_write(linFeatDshp %>% st_simplify(dTolerance = 100), 
# file.path(pthBase, "roads.shp"))
# linFeatDras <- rasterizeLineDensity(linFeatD, anthroDistD)
# writeRaster(linFeatDras, "tests/testthat/data/linFeatTif.tif", overwrite = TRUE)
# 
# eskerDras <- rasterizeLineDensity(st_read("tests/testthat/data/esker.shp"),
#                                   anthroDistD)
# 
# writeRaster(eskerDras, "tests/testthat/data/eskerTif.tif", overwrite = TRUE)

# Test all different ways to run from paths #===================================
pthBase <- system.file("extdata", package = "caribouMetrics")

paths_eskshp_linFshp <- caribouHabitat(
  landCover = file.path(pthBase, "landCover.tif"),
  esker = file.path(pthBase, "esker.shp"),
  natDist = file.path(pthBase, "natDist.tif"),
  anthroDist = file.path(pthBase, "anthroDist.tif"),
  linFeat = file.path(pthBase, "roads.shp"),
  projectPoly = file.path(pthBase, "projectPoly.shp"),
  linFeatSave = file.path(pthBase, "linFeatTif400.tif"),
  eskerSave = file.path(pthBase, "eskerTif400.tif"), 
  caribouRange = "Churchill", 
  winArea = 500
)

paths_list_linF <- caribouHabitat(
  landCover = file.path(pthBase, "landCover.tif"),
  esker = file.path(pthBase, "esker.shp"),
  natDist = file.path(pthBase, "natDist.tif"),
  anthroDist = file.path(pthBase, "anthroDist.tif"),
  linFeat = list(roads = file.path(pthBase, "roads.shp"),
                 rail = file.path(pthBase, "rail.shp"),
                 utilities = file.path(pthBase, "utilities.shp")), 
  projectPoly = file.path(pthBase, "projectPoly.shp"), 
  caribouRange = "Churchill", 
  winArea = 500
)



# Test all different ways to run with data #====================================
landCoverD = terra::rast(file.path(pthBase, "landCover.tif")) %>% 
  reclassPLC()
eskerDras = terra::rast(file.path(pthBase, "eskerTif400.tif"))
eskerDshp = st_read(file.path(pthBase, "esker.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
natDistD = terra::rast(file.path(pthBase, "natDist.tif"))
anthroDistD = terra::rast(file.path(pthBase, "anthroDist.tif"))
linFeatDras = terra::rast(file.path(pthBase, "linFeatTif400.tif"))
projectPolyD = st_read(file.path(pthBase, "projectPoly.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDshp = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
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
  linFeat = list(roads = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
                   st_set_agr("constant"),
                 rail = st_read(file.path(pthBase, "rail.shp"), quiet = TRUE) %>% 
                   st_set_agr("constant"),
                 utilities = st_read(file.path(pthBase, "utilities.shp"), quiet = TRUE)%>% 
                   st_set_agr("constant")), 
  projectPoly = projectPolyD,
  caribouRange = "Churchill", 
  winArea = 500
)

test_that("saveOutput works",{
  f <- file.path(tempdir(), "test.tif")
  caribouHabitat(landCover = landCoverD, esker = eskerDshp, natDist = natDistD, 
                 anthroDist = anthroDistD, 
                 linFeat = linFeatDras, 
                 projectPoly = projectPolyD,
                 caribouRange = "Churchill", 
                 winArea = 500, 
                 saveOutput = f) %>% 
    expect_s4_class("CaribouHabitat")
  unlink(f)
  
  f <- file.path(tempdir(), "test.asc")
  caribouHabitat(landCover = landCoverD, esker = eskerDshp, natDist = natDistD, 
                 anthroDist = anthroDistD, 
                 linFeat = linFeatDras, 
                 projectPoly = projectPolyD,
                 caribouRange = "Churchill", 
                 winArea = 500, 
                 saveOutput = f) %>% 
    expect_s4_class("CaribouHabitat")
  unlink(f)
  
})


test_that("raster road input works as expected", {
  data_list_linFrdRast1 <- caribouHabitat(
    landCover = landCoverD, esker = eskerDshp, natDist = natDistD, 
    anthroDist = anthroDistD, 
    linFeat = list(roads = linFeatDras > 10,
                   rail = st_read(file.path(pthBase, "rail.shp"), quiet = TRUE) %>% 
                     st_set_agr("constant"),
                   utilities = st_read(file.path(pthBase, "utilities.shp"), quiet = TRUE)%>% 
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
                   rail = st_read(file.path(pthBase, "rail.shp"), quiet = TRUE) %>% 
                     st_set_agr("constant"),
                   utilities = st_read(file.path(pthBase, "utilities.shp"), quiet = TRUE)%>% 
                     st_set_agr("constant")), 
    projectPoly = projectPolyD,
    caribouRange = "Churchill", 
    winArea = 500, 
    ptDensity = 2
  )
  expect_gt(data_list_linFrdRast2@linFeat %>% terra::global(max),
            data_list_linFrdRast1@linFeat %>% terra::global(max))
})

test_that("results match when input is paths or data",{
  expect_equal(data_esktif_linFtif@habitatUse, paths_eskshp_linFshp@habitatUse)
  expect_equal(data_eskshp_linFshp@habitatUse, data_esktif_linFtif@habitatUse)
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
    plot(data_esktif_linFtif, raster.title = "Orig")
    plot(data_noAnthro, raster.title = "No Anthro")
    plot(data_noNatDist, raster.title = "No NatDist")
    plot(data_noDist, raster.title = "No Dist")
  }

})

test_that("different coefTable works", {
  w_coefHR <- caribouHabitat(
    landCover = file.path(pthBase, "landCover.tif"),
    esker = file.path(pthBase, "esker.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    linFeat = file.path(pthBase, "roads.shp"),
    projectPoly = file.path(pthBase, "projectPoly.shp"),
    linFeatSave = file.path(pthBase, "linFeatTif400.tif"),
    eskerSave = file.path(pthBase, "eskerTif400.tif"), 
    caribouRange = "Churchill", 
    winArea = 500
  )
  
  coefMod <- coefTableHR %>%
    mutate(Coefficient = ifelse(Variable %in% c("ST", "CON"), -5, Coefficient))
  
  w_coefMod <- caribouHabitat(
    landCover = file.path(pthBase, "landCover.tif"),
    esker = file.path(pthBase, "esker.shp"),
    natDist = file.path(pthBase, "natDist.tif"),
    anthroDist = file.path(pthBase, "anthroDist.tif"),
    linFeat = file.path(pthBase, "roads.shp"),
    projectPoly = file.path(pthBase, "projectPoly.shp"),
    linFeatSave = file.path(pthBase, "linFeatTif400.tif"),
    eskerSave = file.path(pthBase, "eskerTif400.tif"), 
    caribouRange = "Churchill", 
    coefTable = coefMod,
    winArea = 500
  )
  
  expect_false(isTRUE(all.equal(w_coefMod@habitatUse, w_coefHR@habitatUse)))
})

# Compare output raster to previous run. Use raster and not whole object in case
# we want to change the object in future
# Do to changes in CRS order we simply check that the two rasters are
# equivalent to each other
resultCompare <- readRDS(file.path(test_path("data"), "resultCompare.rds"))

# To update
# saveRDS(terra::wrap(data_esktif_linFtif@habitatUse),
#         file.path("tests/testthat/data", "resultCompare.rds"),
#         version = 2)

testthat::test_that("results match previous results",{
  testthat::expect_true(
    terra::all.equal(terra::rast(resultCompare), 
                      data_esktif_linFtif@habitatUse)
  )
})

#tidy created files
file.remove(file.path(pthBase, "linFeatTif400.tif"))
file.remove(file.path(pthBase, "eskerTif400.tif"))



