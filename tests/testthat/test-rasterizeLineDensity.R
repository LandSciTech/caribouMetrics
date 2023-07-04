# create example raster
lc <- terra::rast(xmin = 0, xmax = 25000, ymin = 0, ymax = 25000,
                  resolution = 250, crs = "EPSG:5070")
lc[] <- 1:terra::ncell(lc)

# create line
lf <- sf::st_sf(geometry = sf::st_sfc(list(sf::st_linestring(matrix(c(0, 0, 10000, 10000),
                                                            ncol = 2, byrow = TRUE)),
                                   sf::st_linestring(matrix(c(1, 10001, 10001, 1),
                                                            ncol = 2, byrow = TRUE)),
                                   sf::st_linestring(matrix(c(5001, 10001, 5001, 1),
                                                            ncol = 2, byrow = TRUE))),
                              crs = 5070))


test_that("Line density is correct", {

  linDens <- rasterizeLineDensity(lf, lc)
  maxDens <- terra::global(linDens, max) %>% unlist() %>% unname()
  
  expect_equal(maxDens, 97)
})

test_that("can handle points", {
  lf2 <- lf %>% bind_rows(st_as_sf(data.frame(x = 7000, y = 4000), 
                                   coords = c("x", "y"), crs = 5070))
  
  linDens2 <- rasterizeLineDensity(lf2, lc)
  ptDens <- linDens2[terra::cellFromXY(linDens2, matrix(c(7000, 4000), ncol = 2))][1,1]
  
  expect_equal(ptDens, 40)
})



# comparing methods #=========================
# 
# lc_ras <- as(lc, "Raster")
# 
# Doesn't work any more
# rastLines <- rasterizeLineDensity(lf, lc_ras)
# 
# plot(rastLines)
# 
# # from https://rdrr.io/github/jeffreyevans/spatialEco/src/R/pseudo.absence.R
# raster.as.im <- function(im) {
#   r <- terra::res(im)[1]
#   orig <- as.numeric(sf::st_bbox(im)) + 0.5 * r
#   dm <- dim(im)[2:1]
#   xx <- unname(orig[1] + cumsum(c(0, rep(r[1], dm[1] - 1))))
#   yy <- unname(orig[2] + cumsum(c(0, rep(r[1], dm[2] - 1))))
#   return(spatstat.geom::im(matrix(terra::values(im), ncol = dm[1], 
#                                   nrow = dm[2], byrow = TRUE)[dm[2]:1, ], 
#                            xcol = xx, yrow = yy))
# }
# 
# rastLD_spst <- function(line, ras){
#   spst_im <- spatstat.geom::pixellate(x = spatstat.geom::as.psp(sf::st_geometry(line)),
#                                       W = raster.as.im(ras),
#                                       DivideByPixelArea = F)
#   spst_rast <- terra::rast(spst_im/(terra::res(ras)[1]*terra::res(ras)[2]/10000))
#   spst_rast <- round(spst_rast, digits = 1)
#   terra::crs(spst_rast) <- terra::crs(ras)
#   spst_rast
# }
# 
# plot(rastLD_spst(lf, lc))
# 
# rastLD_terra <- function(line, ras){
#   line_len <- terra::rasterizeGeom(terra::vect(lf), lc, fun = "length")
#   
#   cell_area <- terra::cellSize(lc)/10000 
#   
#   line_len/cell_area
# }
# 
# plot(rastLD_terra(lf, lc))
# 
# bm <- bench::mark(spex = rasterizeLineDensity(lf, lc_ras),
#                   spatstat = rastLD_spst(lf, lc), 
#                   terra = rastLD_terra(lf, lc),
#                   check = F)
# 
# # Conclusion: spatstat.geom method is much faster but terra uses the least
# # memory and is twice as fast as old version. So just keep terra version but
# # keep this code around in case we need a faster version and decide it is
# # worth the extra dependencies
# # # A tibble: 3 × 13
# # expression      min   median itr/se…¹ mem_a…² gc/se…³ n_itr  n_gc total…⁴ result memory     time
# # <bch:expr> <bch:tm> <bch:tm>    <dbl> <bch:b>   <dbl> <int> <dbl> <bch:t> <list> <list>     <list>
# #   1 spex        625.8ms  625.8ms     1.60 16.96MB    1.60     1     1   626ms <NULL> <Rprofmem> <bench_tm>
# #   2 spatstat     46.6ms   48.9ms    19.8   1.38MB    1.80    11     1   555ms <NULL> <Rprofmem> <bench_tm>
# #   3 terra       297.5ms    300ms     3.33 26.08KB    0        2     0   600ms <NULL> <Rprofmem> <bench_tm>
