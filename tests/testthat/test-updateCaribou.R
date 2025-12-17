# test updateLC
context("Test updating process")

pthBase <- system.file("extdata", package = "caribouMetrics")


landCoverD = terra::rast(file.path(pthBase, "landCover.tif")) %>% 
  reclassPLC()
eskerDras = terra::rast(file.path(pthBase, "eskerTif.tif"))
eskerDshp = st_read(file.path(pthBase, "esker.shp"), quiet = TRUE)
natDistD = terra::rast(file.path(pthBase, "natDist.tif"))
anthroDistD = terra::rast(file.path(pthBase, "anthroDist.tif"))
linFeatDras = terra::rast(file.path(pthBase, "linFeatTif.tif"))
projectPolyD = st_read(file.path(pthBase, "projectPoly.shp"), quiet = TRUE)
linFeatDshp = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE)


# newData versions to check updating. 
ext <- terra::ext(natDistD)- 15000

msk <- terra::crop(natDistD, ext) %>% terra::setValues(1) %>% 
  terra::extend(natDistD, fill = 0)

natDistD2 <- terra::mask(natDistD, msk, maskvalue = 1, updatevalue = 1)

linFeatDras2 <- terra::mask(linFeatDras, 
                             terra::crop(linFeatDras, ext) %>% 
                               terra::setValues(1) %>% 
                               terra::extend(linFeatDras, fill = 0), 
                             maskvalue = 1, updatevalue = 0)


# need to change these so can see the difference
anthroDistD <- terra::mask(anthroDistD, msk, maskvalue = 1, updatevalue = 0)
# make a point in the middle to compare
ext2 <- terra::ext(natDistD) %>% as.vector()

pointCompare <- st_sf(ID = 1, 
                      geometry = st_sfc(st_point(c((ext2["xmax"] - ext2["xmin"])/2 +ext2["xmin"], 
                                            (ext2["ymax"] - ext2["ymin"])/2 + ext2["ymin"]))),
                      crs = st_crs(natDistD2))

# process data to pass to updateLC
procedData <- caribouHabitat(landCover = landCoverD,
                             esker = eskerDras,
                             natDist = natDistD,
                             anthroDist = anthroDistD,
                             linFeat = linFeatDras, 
                             projectPoly = projectPolyD, 
                             caribouRange = "Churchill", 
                             winArea = 500)     

test_that("processData works for updated data", {
  expect_warning(updted <- processData(procedData, 
                                       newData = list(natDist = natDistD2)))
  
  expect_true(terra::global(procedData@processedData$DTN != 
                                  updted@processedData$DTN, max, na.rm = TRUE)[1,1] == 1)
  expect_true(terra::global(procedData@processedData$TDENLF == 
                                  updted@processedData$TDENLF, min, na.rm = TRUE)[1,1] == 1)
  
  updted2 <- processData(procedData, newData = list(linFeat = linFeatDras2))
  expect_true(terra::global(procedData@processedData$DTN == 
                                  updted2@processedData$DTN, min, na.rm = TRUE)[1,1] == 1)
  expect_true(terra::global(procedData@processedData$TDENLF != 
                                  updted2@processedData$TDENLF, max, na.rm = TRUE)[1,1] ==1)
    
  expect_true(terra::extract(procedData@linFeat, pointCompare)$linFeatTif != 
                terra::extract(updted2@linFeat, pointCompare)$linFeatTif)
  
  expect_error(processData(procedData, newData = list(natDist2 = natDistD2)),
               "newData must be a named list")
  
})

test_that("process data works for updated data that is not aligned", {
  expect_warning({
    updted <- processData(procedData, 
                        newData = list(natDist = natDistD %>% 
                                         terra::extend(terra::ext(natDistD)+251)))
  })
  expect_equal(procedData@habitatUse, updted@habitatUse)

  expect_error(processData(procedData, 
                           newData = list(natDist = natDistD %>% 
                                            terra::aggregate(4))), 
               "all raster data sets must have matching resolution")

})

test_that("update process result is same with same data", {
  expect_error(updateCaribou(procedData, newData = list(landCoverD)), 
               "newData must be a named list")
  
  update_data <- updateCaribou(procedData, 
                               newData = list(landCover = landCoverD,
                                              natDist = natDistD,
                                              anthroDist = anthroDistD, 
                                              linFeat = linFeatDras))

  #TODO: this is different by a small amount (< 0.01) but only along the edges.
  #Not sure why but also not very important. expVars are the same before the
  #movingWindow but different afterwards. The call is the same though. Have not
  #figured out but not serious
  expect_true(terra::all.equal(procedData@habitatUse$Spring, update_data@habitatUse$Spring))
  
})
