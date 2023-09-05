
mat <- matrix(c(rep(1,5000), rep(2, 5000)), nrow = 100, ncol = 100)
exps <- terra::rast(mat, crs = "EPSG:5070")
lyrs <- terra::segregate(exps)
lyrs[55:60, 20:80] <- NA
pt <- st_sfc(st_point(c(49.5, 50))) %>% st_sf(PID = 1) %>%
  set_names(c("PID", "geometry")) %>% st_set_geometry("geometry") %>% 
  terra::vect()

# res <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = FALSE)
# 
# res2 <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE)

res_rast <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = FALSE)

res2_rast <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE)

test_that("defaults give expected result on dummy data", {
  # result matches expected
  expect_equal(terra::extract(res_rast$X1, pt)[1,2], 13/21)
  
  # Offset results in larger window
  expect_gt(terra::extract(res_rast$X1, pt)[1,2], 
            terra::extract(res2_rast$X1, pt)[1,2])
  
  # edge should be 1
  expect_equal(res_rast[[1]][1,5][1,1], 1)
  #outside NA is not affected
  expect_equal(res_rast[[1]][54, 21][1,1], 1)
  # inside NA is still NA
  expect_true(is.na(res_rast[[1]][55,20]))
})

test_that("other options work",{
  
  res2_rast <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE, 
                               naExternal = "NA", naInternal = "interpolate")
  
  # NA edges leads to edge cells becoming NA and internal NAs are not
  # catching but stay NA
  expect_true(is.na(res2_rast[[1]][1,5]))
  expect_true(!is.na(res2_rast[[1]][54, 21]))
  expect_true(is.na(res2_rast[[1]][55,20]))
  
  res3_rast <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE, 
                               naExternal = "expand")
  
  # expand edges leads to edge cells not being NA and internal NAs are not
  # catching but stay NA
  expect_equal(res3_rast[[1]][1,5][1,1], 1)
  expect_equal(res3_rast[[1]][54,21][1,1], 1)
  expect_true(is.na(res3_rast[[1]][55,20]))
  
  
  res4_rast <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE, 
                               naInternal = "interpolate")
  # there is a little "bleeding" from the interpolation
  expect_gt(res4_rast[[1]][54,55][1,1], 0)
  expect_true(is.na(res4_rast[[1]][55,20]))

  
  res5_rast <- movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE, 
                               naInternal = "zero")
  # NA assumed to be 0
  expect_lt(res5_rast[[1]][54,45][1,1], 1)
  expect_true(is.na(res5_rast[[1]][55,20]))
  
})

# use bench::mark to compare speed
# bm <- bench::mark(pfocal = movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE),
#                   rasterfocal = movingWindowAvg(lyrs, 2.5, c("X1", "X2"), offset = TRUE, 
#                                                 usePfocal = FALSE),
#                   iterations = 10)
# plot(bm)

# # A tibble: 2 x 13
#   expression     min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result
#   <bch:expr>  <bch:> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>
# 1 pfocal      32.8ms 33.2ms      29.9    2.33MB     3.32     9     1      301ms <Rstr~
# 2 rasterfocal 37.7ms 38.7ms      24.1    1.41MB     2.67     9     1      374ms <Rstr~
# # ... with 3 more variables: memory <list>, time <list>, gc <list>

# difference in speed is larger for a bigger raster but the same if there is a
# level2 gc
