context("Test movingWindowAvg")

test_that("movingWindowAvg gives expected result on dummy data", {
  mat <- matrix(c(rep(1,5000), rep(2, 5000)), nrow = 100, ncol = 100)
  exps <- raster(mat, xmn = 0, xmx = 100, ymn = 0, ymx = 100)
  lyrs <- layerize(exps)
  pt <- st_sfc(st_point(c(49.5, 50))) %>% st_sf(PID = 1) %>%
    set_names(c("PID", "geometry")) %>% st_set_geometry("geometry")
  
  # Internal of function helpful for interactive testing
  # cf2 <- focalWeight(lyrs, 2.5, "circle")
  # #cf2[which(cf2 > 0)] <- 1
  # cf2ras <- raster(cf2, xmn = 0.2, xmx = 0.7, ymn = 0.2, ymx = 0.7)
  # 
  # # Need to use sum focal with weights that will give mean instead of mean
  # # because mean counts cells with weight 0 in the denominator and so doesn't
  # # work
  # rast <- velox::velox(lyrs)
  # rast$sumFocal(cf2, bands = 1:2)
  # 
  # out <- rast$extract_points(sp = pt)
  # outRast <- rast$as.RasterLayer()

  res <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = FALSE)
  expect_equal(raster::extract(res$X1, pt), 13/21)
  
  res2 <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE)
  expect_gt(raster::extract(res$X1, pt), raster::extract(res2$X1, pt))
  
  res3 <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = 9, na.rm = TRUE)
  expect_equal(raster::extract(res2$X1, pt), raster::extract(res3$X1, pt))
  
  res4 <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = 16, na.rm = TRUE)
  expect_gt(raster::extract(res2$X1, pt), raster::extract(res4$X1, pt))
  
  res5 <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = FALSE, na.rm = TRUE, pad = TRUE)
  expect_equal(raster::extract(res3$X1, pt), 13/21)
  
  res6 <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE, na.rm = TRUE, pad = TRUE)
  expect_gt(raster::extract(res$X1, pt), raster::extract(res4$X1, pt))
 
  res7 <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = FALSE, na.rm = TRUE, 
                          pad = TRUE, padValue = 0)
  expect_equal(raster::extract(res3$X1, pt), 13/21)
  
  res8 <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE, na.rm = TRUE, 
                          pad = TRUE, padValue = 0)
  expect_gt(raster::extract(res$X1, pt), raster::extract(res4$X1, pt))
  })


