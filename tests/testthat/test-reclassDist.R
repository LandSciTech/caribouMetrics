
pth_base <- system.file("extdata", package = "caribouMetrics")

plcD <- raster::raster(file.path(pth_base, "landCover.tif")) 

fireYr <- sf::st_read(file.path(pth_base, "fireAFFES2020.shp")) 



test_that("simple example works", {
  natDist <- reclassDist(fireYr,
                         endYr = 2010,
                         numCumYrs = 30,
                         template = plcD,
                         dateField = "FIRE_YEAR")
  expect_s4_class(natDist, "Raster")
})

test_that("raster works", {
  fireYrRast <- fasterize::fasterize(fireYr, raster = plcD, field = "FIRE_YEAR")
  #Not working
  natDist2 <- reclassDist(fireYrRast,
                         endYr = 2010,
                         numCumYrs = 30,
                         template = plcD,
                         dateField = "FIRE_YEAR")
  
})
