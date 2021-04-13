
# path to use interactively
# pthBase <- "tests/testthat/data/"
pthBase <- "data/"


landCoverD = raster(paste0(pthBase, "plc", ".tif")) %>% 
  reclassPLC()
eskerDras = raster(paste0(pthBase, "eskerTif", ".tif"))
eskerDshp = st_read(paste0(pthBase, "esker", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
friLUD = read.csv(paste0(pthBase, "friLU", ".csv"), stringsAsFactors = FALSE)
updatedLCD = raster(paste0(pthBase, "fri", ".tif")) %>% 
  reclassFRI(friLUD)
ageD = raster(paste0(pthBase, "age", ".tif"))
natDistD = raster(paste0(pthBase, "natDist", ".tif"))
anthroDistD = raster(paste0(pthBase, "anthroDist", ".tif"))
harvD = raster(paste0(pthBase, "harv", ".tif"))
linFeatDras = raster(paste0(pthBase, "linFeatTif", ".tif"))
projectPolyD = st_read(paste0(pthBase, "projectPoly", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDshp = st_read(paste0(pthBase, "linFeat", ".shp"), quiet = TRUE) %>% 
  st_set_agr("constant")

data_esktif_linFtif <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, updatedLC = updatedLCD, 
  age = ageD, natDist = natDistD,
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = projectPolyD,
  caribouRange = "Churchill", 
  winArea = 500
)

projPolyPts <- projectPolyD$geometry %>% st_cast("POINT") %>% st_coordinates()
twoRange <- st_sfc(st_polygon(list(projPolyPts[c(1, 2, 3 ,1),])), 
                   st_polygon(list(projPolyPts[c(1, 4, 3 ,1),]))) %>%
  st_as_sf() %>% 
  mutate(ID = 1:n(), 
         Range = c("Missisa", "Nipigon")) %>% 
  st_set_crs(st_crs(projectPolyD))

landCoverD2 <- raster::merge(landCoverD, 
                             raster::shift(landCoverD, dx = 19851.6, 
                                           dy = raster::nrow(landCoverD)*
                                             raster::xres(landCoverD)))
plot(landCoverD2)
# side by side polygons that might end up with different extent/origin
# poly1 <- raster::drawPoly()
# poly1 <- st_as_sf(poly1)

poly1 <-
  structure(
    list(
      geometry = structure(
        list(structure(
          list(structure(
            c(
              375883.834871147,
              382980.752147705,
              388658.285968951,
              376416.103666888,
              375883.834871147,
              12729145.8750217,
              12730387.8355451,
              12726661.9539749,
              12723290.9182686,
              12729145.8750217
            ),
            .Dim = c(5L, 2L)
          )), class = c("XY", "POLYGON",
                        "sfg")
        )),
        class = c("sfc_POLYGON", "sfc"),
        precision = 0,
        bbox = structure(
          c(
            xmin = 375883.834871147,
            ymin = 12723290.9182686,
            xmax = 388658.285968951,
            ymax = 12730387.8355451
          ),
          class = "bbox"
        ),
        crs = structure(list(input = NA_character_,
                             wkt = NA_character_), class = "crs"),
        n_empty = 0L
      )
    ),
    row.names = 1L,
    class = c("sf",
              "data.frame"),
    sf_column = "geometry",
    agr = structure(
      integer(0),
      class = "factor",
      .Label = c("constant",
                 "aggregate", "identity"),
      .Names = character(0)
    )
  )
# poly2 <- raster::drawPoly()
# poly2 <- st_as_sf(poly2)
poly2 <-
  structure(
    list(
      geometry = structure(
        list(structure(
          list(structure(
            c(
              396819.740836993,
              402497.274658239,
              406578.00209226,
              402142.428794411,
              395045.511517853,
              396819.740836993,
              12748839.8204642,
              12749726.9351238,
              12748130.1287365,
              12739968.6738685,
              12741920.3261195,
              12748839.8204642
            ),
            .Dim = c(6L,
                     2L)
          )), class = c("XY", "POLYGON", "sfg")
        )),
        class = c("sfc_POLYGON",
                  "sfc"),
        precision = 0,
        bbox = structure(
          c(
            xmin = 395045.511517853,
            ymin = 12739968.6738685,
            xmax = 406578.00209226,
            ymax = 12749726.9351238
          ),
          class = "bbox"
        ),
        crs = structure(list(input = NA_character_,
                             wkt = NA_character_), class = "crs"),
        n_empty = 0L
      )
    ),
    row.names = 1L,
    class = c("sf",
              "data.frame"),
    sf_column = "geometry",
    agr = structure(
      integer(0),
      class = "factor",
      .Label = c("constant",
                 "aggregate", "identity"),
      .Names = character(0)
    )
  )
twoRange2 <- rbind(poly1, poly2) %>% mutate(Range = c("Missisa", "Nipigon")) %>% 
  st_set_crs(st_crs(landCoverD2))
# supply polygon with multiple ranges

eskerD2 <- eskerDras %>% raster::merge(raster::shift(eskerDras, dx = 400*49, 
                                                     dy = 400*40)) %>% 
  raster::projectRaster(to = landCoverD2)
linFeatD2 <- linFeatDras %>% 
  raster::merge(raster::shift(linFeatDras, dx = 400*49, dy = 400*40)) %>% 
  raster::projectRaster(to = landCoverD2)


# same coefficients as range
resTwoRange <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, updatedLC = updatedLCD, 
  age = ageD, natDist = natDistD,
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = twoRange,
  caribouRange = data.frame(Range = c("Missisa", "Nipigon"), 
                            coefRange = c("Missisa", "Nipigon"), 
                            stringsAsFactors = FALSE), 
  winArea = 500
)

# different coefficients as range
resTwoRangeDif <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, updatedLC = updatedLCD, 
  age = ageD, natDist = natDistD,
  anthroDist = anthroDistD, harv = harvD,
  linFeat = linFeatDras, projectPoly = twoRange,
  caribouRange = data.frame(Range = c("Missisa", "Nipigon"), 
                            coefRange = c("Nipigon", "Missisa"), 
                            stringsAsFactors = FALSE), 
  winArea = 500
)

# different winArea need to supply dif coef table to work on small example
coefTableSmall <- coefTableHR %>% mutate(WinArea = WinArea/10)

resTwoRangeDifWin <- caribouHabitat(
  landCover = landCoverD2, esker = eskerD2, 
  linFeat = linFeatD2, projectPoly = twoRange2,
  caribouRange = data.frame(Range = c("Missisa", "Nipigon"), 
                            coefRange = c("Nipigon", "Missisa"), 
                            stringsAsFactors = FALSE), 
  coefTable = coefTableSmall
)

resRange1 <- caribouHabitat(
  landCover = landCoverD2, esker = eskerD2, 
  linFeat = linFeatD2, projectPoly = twoRange2 %>% slice(1),
  caribouRange = data.frame(Range = c("Missisa"), 
                            coefRange = c("Nipigon"), 
                            stringsAsFactors = FALSE), 
  coefTable = coefTableSmall,
  padProjPoly = TRUE
)

resRange2 <- caribouHabitat(
  landCover = landCoverD2, esker = eskerD2, 
  linFeat = linFeatD2, projectPoly = twoRange2 %>% slice(2),
  caribouRange = data.frame(Range = c("Nipigon"), 
                            coefRange = c("Missisa"), 
                            stringsAsFactors = FALSE), 
  coefTable = coefTableSmall,
  padProjPoly = TRUE
)

test_that("results for two ranges done separately same as done together",{
  ext1 <- raster::extent(resRange1@habitatUse)
  pointCompare <- st_sf(ID = 1, 
                        geometry = st_sfc(st_point(c((ext1@xmax - ext1@xmin)/2 + ext1@xmin, 
                                                     (ext1@ymax - ext1@ymin)/2 + ext1@ymin))),
                        crs = st_crs(resRange1@habitatUse))
  expect_lt(
  raster::extract(resRange1@habitatUse$Fall, pointCompare)-
    raster::extract(resTwoRangeDifWin@habitatUse$Fall, pointCompare),
  0.01
  )
  
  ext2 <- raster::extent(resRange2@habitatUse)
  pointCompare <- st_sf(ID = 1, 
                        geometry = st_sfc(st_point(c((ext2@xmax - ext2@xmin)/2 + ext2@xmin, 
                                                     (ext2@ymax - ext2@ymin)/2 + ext2@ymin))),
                        crs = st_crs(resRange2@habitatUse))
  expect_lt(
    raster::extract(resRange2@habitatUse$Fall, pointCompare)-
      raster::extract(resTwoRangeDifWin@habitatUse$Fall, pointCompare),
    0.01
  )
  

})

# compare plots
if(interactive()){
  plot(resTwoRangeDifWin@processedData$CON)
  plot(resRange1@processedData$CON, add = T)
  plot(resRange2@processedData$CON, add = T)
}

