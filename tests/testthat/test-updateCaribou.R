# test updateLC
context("Test updating process")

pthBase <- system.file("extdata", package = "caribouMetrics")


landCoverD = raster(file.path(pthBase, "landCover.tif")) %>% 
  reclassPLC()
eskerDras = raster(file.path(pthBase, "eskerTif.tif"))
eskerDshp = st_read(file.path(pthBase, "esker.shp"), quiet = TRUE)
natDistD = raster(file.path(pthBase, "natDist.tif"))
anthroDistD = raster(file.path(pthBase, "anthroDist.tif"))
linFeatDras = raster(file.path(pthBase, "linFeatTif.tif"))
projectPolyD = st_read(file.path(pthBase, "projectPoly.shp"), quiet = TRUE)
linFeatDshp = st_read(file.path(pthBase, "linFeat.shp"), quiet = TRUE)


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


# need to change these so can see the difference
anthroDistD <- raster::mask(anthroDistD, msk, maskvalue = 1, updatevalue = 0)
# make a point in the middle to compare
ext2 <- raster::extent(natDistD)

pointCompare <- st_sf(ID = 1, 
                      geometry = st_sfc(st_point(c((ext2@xmax - ext2@xmin)/2 + ext2@xmin, 
                                            (ext2@ymax - ext2@ymin)/2 + ext2@ymin))),
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
  
  expect_true(raster::cellStats(procedData@processedData$DTN != 
                                  updted@processedData$DTN, max) == 1)
  expect_true(raster::cellStats(procedData@processedData$TDENLF == 
                                  updted@processedData$TDENLF, min) == 1)
  
  updted2 <- processData(procedData, newData = list(linFeat = linFeatDras2))
  expect_true(raster::cellStats(procedData@processedData$DTN == 
                                  updted2@processedData$DTN, min) == 1)
  expect_true(raster::cellStats(procedData@processedData$TDENLF != 
                                  updted2@processedData$TDENLF, max) ==1)
    
  expect_true(raster::extract(procedData@linFeat, pointCompare) != 
                raster::extract(updted2@linFeat, pointCompare))
  
  expect_error(processData(procedData, newData = list(natDist2 = natDistD2)),
               "newData must be a named list")
  
})

test_that("process data works for updated data that is not aligned", {
  expect_warning({
    updted <- processData(procedData, 
                        newData = list(natDist = natDistD %>% 
                                         raster::extend(raster::extent(natDistD)+251)))
  })
  expect_equal(procedData@habitatUse, updted@habitatUse)

  expect_error(processData(procedData, 
                           newData = list(natDist = natDistD %>% 
                                            raster::aggregate(4))), 
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
  #TODO: this is different by a small amount but only along the edges. Not sure
  #why but also not very important
  expect_equal(procedData@habitatUse,
               update_data@habitatUse,
               label = "",
               expected.label = "different by a small amount but only along the edges. Not sure why but also not very important")
  
})
