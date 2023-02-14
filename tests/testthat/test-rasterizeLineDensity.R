test_that("Line density is correct", {
  # create example raster
  lc <- raster::raster(xmn = 0, xmx = 25000, ymn = 0, ymx = 25000,
                       resolution = 250, crs = 5070)
  lc[] <- 1:ncell(lc)

  # create line
  lf <- sf::st_as_sf(sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10000, 10000),
                                                              ncol = 2, byrow = TRUE)),
                                     sf::st_linestring(matrix(c(1, 10001, 10001, 1),
                                                              ncol = 2, byrow = TRUE)),
                                     sf::st_linestring(matrix(c(5001, 10001, 5001, 1),
                                                              ncol = 2, byrow = TRUE))),
                                     crs = 5070))

  rastLines <- rasterizeLineDensity(lf, lc)

  plot(rastLines)
  
  spst_im <- spatstat.geom::pixellate(x = spatstat.geom::as.psp(sf::st_geometry(lf)), 
                                      W = maptools::as.im.RasterLayer(lc),
                                      DivideByPixelArea = F)
  spst_rast <- raster::raster(spst_im)/(res(lc)[1]*res(lc)[2]/10000) 
  spst_rast <- round(spst_rast, digits = 1)
  
  bm <- bench::mark(spex = rasterizeLineDensity(lf, lc),
                    spatstat = {  spst_im <- spatstat.geom::pixellate(x = spatstat.geom::as.psp(st_geometry(lf)), 
                                                                      W = maptools::as.im.RasterLayer(lc),
                                                                      DivideByPixelArea = F)
                    spst_rast <- raster::raster(spst_im)/(res(lc)[1]*res(lc)[2]/10000) 
                    spst_rast <- round(spst_rast, digits = 1)
                    spst_rast <- raster::`crs<-`(spst_rast, value = raster::crs(lc))}, check = F)
})
