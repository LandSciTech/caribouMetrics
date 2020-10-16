#' Create nested hexagonal grid
#'
#' Takes an existing hexagonal grid and builds larger grids by joining the
#' neighbouring hexs to a center hex to make a larger hexagon.
#'
#' @param x sf object of starting grid
#' @param scale the number of neighbouring hexagons to join, ie 1 for the six
#'   neighbouring hexs and 2 for their neighbours
#'   

createNestedHexGrid <- function(x, scale){
  
  # Find out how many hexs are in a row ie ncol
  # a = side length of hexs
  A <- st_area(x %>% slice(1))
  a <- sqrt((2*A)/(3*sqrt(3)))
  
  # width of grid = w
  bbox <- st_bbox(x)
  w <- bbox["xmax"] - bbox["xmin"]
  h <- bbox["ymax"] - bbox["ymin"]
  
  # One grid cell is 3*a
  ncol1 <- trunc(w/(3*a)) %>% units::drop_units()
  ncol2 <- ifelse(round(w/(3*a)) < w/(3*a), ncol1, ncol1+1)
  ncell2Rows <- ncol1 + ncol2 
  nrowx <- trunc(h/(sqrt(3)/2 * a)) %>% units::drop_units()
  
  # Figure out number of grid cells between each hex that will be a centre
  nCellSpace <- 1+ 6*sum(1:scale)
  
  if(ncol1 < nCellSpace){
    stop("the hex grid supplied is not wide enough for the selected scale")
  }
  
  # pattern is 1 + 6 + 2*6 + 3*6 + 4*6
  # c(7, 19, 37, 61)
  
  # width of 1 big hex
  expan <- scale * 2 + 1
  
  polyIDsA <- vector(length = ncell(x)/nCellSpace/scale)
  initial <- seq(1, by = nCellSpace, length.out = ncol1/nCellSpace)
  polyIDsA[1:length(initial)] <- initial
  for (i in 2:nrowx) {
    even <- !as.logical(i%%2)
    endCellpre <- round(ncell2Rows*((i-1)/1.999999999))
    endCell <- ncell2Rows*((i)/2)
    if(even){
      startCell <- nCellSpace + endCellpre + 1 - scale - expan*(i/2 - 1) - 
        abs(ncol1 - ncol2)
      while(startCell <= endCellpre){
        startCell <- startCell+nCellSpace
      }
      polyIDsi <- seq(startCell,
                      by = nCellSpace,
                      length.out = (endCell - startCell)/nCellSpace )
    } else {
      startCell <- nCellSpace + endCellpre -expan*(i/2-1) -(scale-0.5)
      while(startCell <= endCellpre){
        startCell <- startCell+nCellSpace
      }
      polyIDsi <- seq(startCell,
                      by = nCellSpace,
                      length.out =(endCell - startCell)/nCellSpace )
    }
    
    if(tail(polyIDsi, 1)+nCellSpace <= endCell){
      polyIDsi <- c(polyIDsi, tail(polyIDsi, 1)+nCellSpace)
    }
    lastID <- max(which(polyIDsA>0))
    nextID <- length(polyIDsi)+lastID
    polyIDsA[(lastID+1):nextID] <- polyIDsi
  }
  
  startHexs <- x %>% slice(polyIDsA) %>% select(PID)
  
  #qtm(x)+qtm(startHexs, borders = "red")+qtm(out2, fill = NULL, borders = "blue")
  
  if(scale == 1){
    bigHexs <- {st_join(x, startHexs, st_intersects)} %>% 
      filter(!is.na(PID.y)) %>% 
      group_by(PID.y) %>% 
      summarise(PID1 = sum(PID.y, na.rm = TRUE)) %>% 
      transmute(PID1 = 1:n())
  }
  if(scale == 2){
    bigHexs <- st_join(x, startHexs, st_intersects) %>% 
      filter(!is.na(PID.y)) %>% 
      {st_join(x, ., st_intersects)} %>%
      filter(!is.na(PID.y)) %>% 
      group_by(PID.y) %>% 
      summarise(PID2 = sum(PID.y, na.rm = TRUE)) %>% 
      transmute(PID2 = 1:n())
  }
  if(scale == 3){
    bigHexs <- st_join(x, startHexs, st_intersects) %>%
      filter(!is.na(PID.y)) %>%
      {st_join(x, ., st_intersects)} %>%
      filter(!is.na(PID.y)) %>%
      {st_join(x, ., st_intersects)} %>%
      filter(!is.na(PID.y)) %>%
      group_by(PID.y.1) %>%
      summarise(PID3 = sum(PID.y, na.rm = TRUE)) %>%
      transmute(PID3 = 1:n())
  }
  
  return(bigHexs)
}

# # testing createNestedHexGrid
# {
#   library(raster); library(tidyverse); library(sf); library(tmap)
#   tmap_mode("view")
#   baseDir =  "./inputNV"
#   
#   HRDataPath = paste0(baseDir,"/HornsethRempelInfo")
#   plc = raster(paste0(HRDataPath,"/Provincial-Landcover-2000/studyAreaMerged250.tif"))
#   plot(plc)
#   #plot(testGrid)
#   ext <- drawExtent()
#   box <- st_bbox(ext) %>% st_as_sfc() %>% st_as_sf()
#   
#   hexArea <- 160000 # in map units. i.e. metres^2
#   hexShortDia <-  sqrt(2*hexArea/sqrt(3))
#   
#   testGrid <- st_make_grid(box,
#                            cellsize = hexShortDia,
#                            offset = box %>% st_bbox() %>%
#                              .[c("xmin", "ymin")] + c(126,148),
#                            square = FALSE,
#                            flat_topped = TRUE) %>%
#     {data.frame(PID = 1:length(.), geometry = .)} %>%
#     st_as_sf(sf_column_name = "geometry") %>%
#     st_set_crs(st_crs(plc))
#   
#   # grid can some times end up with different ncols in odd vs even rows
#   # make second grid slightly bigger to test
#   
#   ext2 <- ext
#   ext2@xmin <- ext2@xmin-(hexShortDia*1.25)
#   box2 <- st_bbox(ext2) %>% st_as_sfc() %>% st_as_sf()
#   testGrid2 <- st_make_grid(box2,
#                             cellsize = hexShortDia,
#                             offset = box %>% st_bbox() %>%
#                               .[c("xmin", "ymin")] + c(126,148),
#                             square = FALSE,
#                             flat_topped = TRUE)%>%
#     {data.frame(PID = 1:length(.), geometry = .)} %>%
#     st_as_sf(sf_column_name = "geometry") %>%
#     st_set_crs(st_crs(plc))
#   
#   x <- testGrid2
#   scale <- 2 # ie combine 3 cells out from centre ie 592 ha
#   
#   out1 <- createNestedHexGrid(x, 1)
#   out2 <- createNestedHexGrid(x, 2)
#   out3 <- createNestedHexGrid(x, 3)
#   beepr::beep()
#   qtm(out1, fill = NULL)
#   qtm(out2, fill = NULL)
#   qtm(out3, fill = NULL)
#   
#   # make grid then summarize results through spatial join. This makes more sense
#   # than summarising the data when making the larger hex grid because the
#   # spatial join is faster than building the hexgrid
#   outReal1 <- createNestedHexGrid(caribouResultsSF, 1)
#   
#   caribiouResults112 <- st_join(outReal1, caribouResultsSF, join = st_contains) %>% 
#     st_drop_geometry() %>% 
#     group_by(PID1) %>% 
#     summarise_at(.vars = vars(Fall, Winter, Spring, Summer),
#                  mean, na.rm =TRUE)
# }
