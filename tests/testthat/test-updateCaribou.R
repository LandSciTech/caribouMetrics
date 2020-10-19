# test updatePLC
context("Test updating process")
#pthBase <- "tests/testthat/data/"
pthBase <- "data/"


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
friLUD = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE)

# newData versions to check updating. 
ext <- raster::extent(natDistD)- 15000

msk <- raster::crop(natDistD, ext) %>% raster::setValues(1) %>% 
  raster::extend(natDistD, value = 0)

natDistD2 <- raster::mask(natDistD, msk, maskvalue = 1, updatevalue = 1)
ageD2 <- raster::mask(ageD, msk, maskvalue = 1, updatevalue = 1)
ageD3 <- raster::mask(ageD, msk, maskvalue = 1, updatevalue = 50)
harvD2 <- raster::mask(harvD, msk, maskvalue = 1, updatevalue = 1)
linFeatDras2 <- raster::mask(linFeatDras, msk, maskvalue = 1, updatevalue = 0)


# need to change these so can see the difference with harvD2
harvD <- raster::mask(harvD, msk, maskvalue = 1, updatevalue = 0)
anthroDistD <- raster::mask(anthroDistD, msk, maskvalue = 1, updatevalue = 0)
# make a point in the middle to compare
ext2 <- raster::extent(natDistD)

pointCompare <- st_sf(ID = 1, 
                      geometry = st_sfc(st_point(c((ext2@xmax - ext2@xmin)/2 + ext2@xmin, 
                                            (ext2@ymax - ext2@ymin)/2 + ext2@ymin))),
                      crs = st_crs(natDistD2))

# process data to pass to updatePLC
procedData <- caribouHabitat(plcD, eskerDras, friD, ageD, natDistD, 
                             anthroDistD, harvD, linFeatDras, projectPolyD, 
                             friLU = friLUD, 
                             caribouRange = "Churchill")     

# Not using for now
# test_that("updatePLC works with different newData", {
# 
#   # get resTypeLevels from procedData@plc RAT
#   resTypeLevels <- procedData@plc %>% raster::levels() %>% .[[1]] %>% 
#     set_names(c("code", "resType")) %>% 
#     mutate(resType = as.character(resType))
#   
#   updted <- updatePLC(procedData, newData = list(natDist = natDistD2,
#                                                  age = ageD2),
#                       resTypeLevels = resTypeLevels)
#   
#   expect_equal(resTypeLevels %>% filter(resType == "DTN") %>% pull(code),
#                raster::extract(updted$plc, pointCompare))
#   
#   updted2 <- updatePLC(procedData, newData = list(harv = harvD2), 
#                        resTypeLevels = resTypeLevels)
#   # TODO: update based on changes
#   # expect_equal(resTypeLevels %>% filter(resType == "ATN") %>% pull(code),
#   #              raster::extract(updted2$plc, pointCompare))
#   
#   updted3 <- updatePLC(procedData, newData = list(natDist = natDistD2,
#                                                   age = ageD3,
#                                                   fri = procedData@fri), 
#                        resTypeLevels = resTypeLevels)
#   
#   expect_equal(resTypeLevels %>% filter(resType == "CON") %>% pull(code),
#                raster::extract(updted3$plc, pointCompare))
# 
# })

test_that("processData works for updated data", {
  updted <- processData(procedData, newData = list(natDist = natDistD2,
                                                   age = ageD2), 
                        caribouRange = "Churchill", 
                        friLU = friLUD)
  expect_true(raster::cellStats(procedData@processedData$DTN != 
                                  updted@processedData$DTN, max) == 1)
  expect_true(raster::cellStats(procedData@processedData$TDENLF == 
                                  updted@processedData$TDENLF, min) == 1)
  
  updted2 <- processData(procedData, newData = list(linFeat = linFeatDras2), 
                        caribouRange = "Churchill", 
                        friLU = friLUD)
  expect_true(raster::cellStats(procedData@processedData$DTN == 
                                  updted2@processedData$DTN, min) == 1)
  expect_true(raster::cellStats(procedData@processedData$TDENLF != 
                                  updted2@processedData$TDENLF, max) ==1)
    
  expect_true(raster::extract(procedData@linFeat, pointCompare) != 
                raster::extract(updted2@linFeat, pointCompare))
  
  expect_error(processData(procedData, newData = list(natDist = natDistD2,
                                                      age2 = ageD2), 
                           caribouRange = "Churchill", 
                           friLU = friLUD),
               "newData must be a named list")
  
  expect_error(processData(procedData, newData = list(natDist = natDistD2), 
                           caribouRange = "Churchill", 
                           friLU = friLUD),
               "both natDist and age must be provided")
  
  expect_error(processData(procedData, newData = list(fri = natDistD2), 
                           caribouRange = "Churchill", 
                           friLU = friLUD),
               "to use fri data either harv")
})

test_that("process data works for updated data that is not aligned", {
  updted <- processData(procedData, 
                        newData = list(harv = harvD %>% 
                                         raster::extend(raster::extent(harvD)+251)),
                        caribouRange = "Churchill", 
                        friLU = friLUD)
  expect_equal(procedData@habitatUse, updted@habitatUse)

  expect_error(processData(procedData, 
                           newData = list(harv = harvD %>% 
                                            raster::aggregate(4)),
                           caribouRange = "Churchill", 
                           friLU = friLUD), 
               "different extent")
})

test_that("update process result is same with same data", {
  expect_error(updateCaribou(procedData, newData = list(friD), 
                             friLU = read.csv(paste0(pthBase, "friLU", ".csv"), 
                                              stringsAsFactors = FALSE),
                             caribouRange = "Churchill"), 
               "newData must be a named list")
  
  update_data <- updateCaribou(procedData, 
                               newData = list(fri = friD, harv = harvD, 
                                              age = ageD, natDist = natDistD, 
                                              anthroDist = anthroDistD, 
                                              linFeat = linFeatDras), 
                               friLU = read.csv(paste0(pthBase, "friLU", ".csv"), 
                                                stringsAsFactors = FALSE),
                               caribouRange = "Churchill")
  
  expect_equal(procedData@habitatUse, 
               update_data@habitatUse)
  
})