context("test combineLinFeat function")

pthBase <- system.file("extdata", package = "caribouMetrics")

roads <- st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")

# these are actually empty
rail <- st_read(file.path(pthBase, "rail.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
utilities <- st_read(file.path(pthBase, "utilities.shp"), quiet = TRUE)%>% 
  st_set_agr("constant")

rail <- st_sf(geometry = st_sfc(st_bbox(roads) %>% matrix(ncol = 2, byrow = T) %>%
                     st_linestring())) %>% 
  st_set_crs(st_crs(roads))
utilities <- st_sf(geometry = st_sfc(st_bbox(roads) %>% 
                                       {c(.["xmin"], .["ymax"],
                                         .["xmax"], .["ymin"])} %>% 
                                       matrix(ncol = 2, byrow = T) %>%
                                       st_linestring()))%>% 
  st_set_crs(st_crs(roads))

linFeatDras <- terra::rast(file.path(pthBase, "linFeatTif.tif"))

roadsSp <- sf::as_Spatial(roads)


test_that("results are same with different input formats",{
  out1 <- combineLinFeat(lst(roads, rail, utilities))
  
  out2 <- combineLinFeat(lst(roadsSp, rail, utilities))
  
  out3 <- combineLinFeat(lst(file.path(pthBase, "roads.shp"),
                             rail,
                             utilities))
  expect_equal(out1, out2)
  expect_equal(out1, out3)
  
  denOut1 <- rasterizeLineDensity(out1 %>% st_transform(st_crs(linFeatDras)),
                                                       linFeatDras)
  
  if(interactive()){
    plot(out1)
    plot(denOut1)
  }
})

test_that("raster input works as expected", {
  out4 <- combineLinFeat(lst(linFeatDras >10, rail, utilities))
  
  expect_s3_class(out4, "sf")
  
  denOut4 <- rasterizeLineDensity(out4, linFeatDras)
  
  if(interactive()){
    plot(out4)
    
    plot(denOut4)
  }
})
