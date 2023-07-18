pthBase <- system.file("extdata", package = "caribouMetrics")

# load example data
plcD = terra::rast(file.path(pthBase, "landCover.tif")) # Defines the study area - NA values are omitted from calculation, everything else is included.
natDistD = terra::rast(file.path(pthBase, "natDist.tif"))
anthroDistD = terra::rast(file.path(pthBase, "anthroDist.tif"))
projectPolyD = st_read(file.path(pthBase, "projectPoly.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDshp = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
roadsD = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
railD = st_read(file.path(pthBase, "rail.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
utilitiesD = st_read(file.path(pthBase, "utilities.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDras = terra::rast(file.path(pthBase, "linFeatTif.tif"))

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

anthroDistD2 <- terra::mask(anthroDistD, msk, maskvalue = 1, updatevalue = 0)

dm <- disturbanceMetrics(
  landCover = plcD,
  natDist = natDistD, 
  anthroDist = anthroDistD, 
  linFeat = linFeatDras, 
  projectPoly = projectPolyD,
  padFocal = FALSE, # assume data outside area is 0 for all variables
  bufferWidth = 500
)

test_that("updateDisturbance works", {
  expect_message(
    dm_upd_lf <- updateDisturbance(dm, newData = list(linFeat = linFeatDras2)),
    "buffering"
  )
  
  expect_gt(results(dm)$Anthro, results(dm_upd_lf)$Anthro)
  
  expect_message(
    dm_upd_an <- updateDisturbance(dm, newData = list(anthroDist = anthroDistD2)),
    "buffering"
  )
  
  expect_gt(results(dm)$Anthro, results(dm_upd_an)$Anthro)
  

  dm_upd_nd <- updateDisturbance(dm, newData = list(natDist = natDistD2))
  
  expect_lt(results(dm)$Fire, results(dm_upd_nd)$Fire)
  expect_equal(results(dm)$Anthro, results(dm_upd_nd)$Anthro)
  
  dm_upd_all <- updateDisturbance(dm, newData = list(linFeat = linFeatDras2,
                                                   anthroDist = anthroDistD2,
                                                   natDist = natDistD2))
  
})

test_that("updateDisturbance works with nondefault inputs", {
  expect_message(
    dm_upd_lfshp <- updateDisturbance(dm, newData = list(linFeat = linFeatDshp)),
    "buffering"
  )
  expect_error(
    updateDisturbance(dm, newData = list(linFeat2 = linFeatDras2, 
                                         natDist = natDistD)),
    "named list"
  )
  
  expect_message(
    dm_upd_lfnotPerc <- updateDisturbance(dm,
                                          newData = list(linFeat = linFeatDras2),
                                          isPercent = FALSE),
    "buffering anthropogenic disturbance"
  )
  
  expect_true(all(dm_upd_lfnotPerc@disturbanceMetrics[1,2:5] < 1))
  
  expect_message(
    dm_upd_lf_bufsf <- updateDisturbance(dm, 
                                         newData = list(linFeat = linFeatDshp), 
                                         linBuffMethod = "sf"),
    "buffering linear features"
  )
  
})

# run benchmark to see if it is actually speeding things up

# bm3 <- bench::mark(
#   distMet = dm2 <- disturbanceMetrics(
#     landCover = plcD,
#     natDist = natDistD,
#     anthroDist = anthroDistD,
#     linFeat = linFeatDras,
#     projectPoly = projectPolyD,
#     padFocal = FALSE, # assume data outside area is 0 for all variables
#     bufferWidth = 500
#   ),
#   updateAnthro = updateDisturbance(dm2, newData = list(anthroDist = anthroDistD)),
#   updateLF = updateDisturbance(dm2, newData = list(linFeat = linFeatDras)),
#   updateNat = updateDisturbance(dm2, newData = list(natDist = natDistD)),
#   updateAll = updateDisturbance(dm2, newData = list(linFeat = linFeatDras,
#                                                    anthroDist = anthroDistD,
#                                                    natDist = natDistD)),
#   min_iterations = 5,
#   filter_gc = FALSE,
#   check = FALSE
# )

# bm3
# # A tibble: 5 x 13
# expression          min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
# <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
# 1 distMet         2.12s    2.16s     0.463    69.8MB        0     5     0     10.81s
# 2 updateAnthro 773.72ms 799.51ms     1.23     59.3MB        0     5     0      4.05s
# 3 updateLF     774.43ms 780.72ms     1.28     59.2MB        0     5     0      3.91s
# 4 updateNat    704.12ms 709.01ms     1.41     51.9MB        0     5     0      3.55s
# 5 updateAll       1.71s    1.73s     0.573    66.3MB        0     5     0      8.73s