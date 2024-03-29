
pth_base <- system.file("extdata", package = "caribouMetrics")

plcD <- terra::rast(file.path(pth_base, "landCover.tif")) 

fireYr <- sf::read_sf(file.path(pth_base, "fireAFFES2020.shp")) 

test_that("all options work", {
  natDist <- reclassDist(fireYr,
                         endYr = 2010,
                         numCumYrs = 30,
                         template = plcD,
                         dateField = "FIRE_YEAR")
  
  expect_s4_class(natDist, "SpatRaster")
  
  fireYrRast <- terra::rasterize(terra::vect(fireYr), plcD, field = "FIRE_YEAR", 
                                     background = 0)
  # with raster year input
  natDist2 <- reclassDist(fireYrRast,
                          endYr = 2010,
                          numCumYrs = 30,
                          template = plcD,
                          dateField = "FIRE_YEAR")
  
  expect_true(terra::global(natDist - natDist2, "max") == 0)
  
  # time since dist polygons
  natDist3 <- reclassDist(fireYr %>% mutate(FIRE_YEAR = 2010 - FIRE_YEAR ),
                          endYr = 0,
                          numCumYrs = 30,
                          template = plcD,
                          dateField = "FIRE_YEAR")
  
  expect_true(terra::global(natDist - natDist3, "max") == 0)
  
  # time since dist raster
  natDist4 <- reclassDist(2010 - fireYrRast,
                          endYr = 0,
                          numCumYrs = 30,
                          template = plcD,
                          dateField = "FIRE_YEAR")
  
  expect_true(terra::global(natDist - natDist4, "max") == 0)
  
  # if it was na in input raster should come out as 0
  fireYrRast2 <- terra::rasterize(terra::vect(fireYr), plcD, field = "FIRE_YEAR", 
                       background = NA)
  
  natDist5 <- reclassDist(2010 - fireYrRast2,
                          endYr = 0,
                          numCumYrs = 30,
                          template = plcD,
                          dateField = "FIRE_YEAR")
  
  expect_true(terra::global(is.na(natDist5), "max") == 0)
})


