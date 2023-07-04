
# path to use interactively
# pthBase <- "tests/testthat/data/"
pthBase <- system.file("extdata", package = "caribouMetrics")


landCoverD = terra::rast(file.path(pthBase, "landCover.tif")) %>% 
  reclassPLC()
eskerDras = terra::rast(file.path(pthBase, "eskerTif.tif"))
eskerDshp = st_read(file.path(pthBase, "esker.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")

natDistD = terra::rast(file.path(pthBase, "natDist.tif"))
anthroDistD = terra::rast(file.path(pthBase, "anthroDist.tif"))

linFeatDras = terra::rast(file.path(pthBase, "linFeatTif.tif"))
projectPolyD = st_read(file.path(pthBase, "projectPoly.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")
linFeatDshp = st_read(file.path(pthBase, "roads.shp"), quiet = TRUE) %>% 
  st_set_agr("constant")

data_esktif_linFtif <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras,  natDist = natDistD,
  anthroDist = anthroDistD, linFeat = linFeatDras, projectPoly = projectPolyD,
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

landCoverD2 <- terra::merge(landCoverD, 
                             terra::shift(landCoverD, 
                                           dx = terra::nrow(landCoverD)*
                                             terra::xres(landCoverD), 
                                           dy = terra::nrow(landCoverD)*
                                             terra::xres(landCoverD)))
# plot(landCoverD2)
# side by side polygons that might end up with different extent/origin
# poly3 <- terra::draw()
# poly3 <- st_as_sf(poly3)
# dput(poly3)

poly1 <- {structure(
    list(
      geometry = structure(
        list(structure(
          list(structure(
            c(
              700170.862790608,
              751689.250606358,
              749481.319699969,
              701642.816728201,
              700170.862790608,
              12619084.4908842,
              12617612.5369466,
              12581549.6654756,
              12584493.5733508,
              12619084.4908842
            ),
            .Dim = c(5L, 2L)
          )), class = c("XY", "POLYGON",
                        "sfg")
        )),
        class = c("sfc_POLYGON", "sfc"),
        precision = 0,
        bbox = structure(
          c(
            xmin = 700170.862790608,
            ymin = 12581549.6654756,
            xmax = 751689.250606358,
            ymax = 12619084.4908842
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
  )}

poly2 <- {structure(
    list(
      geometry = structure(
        list(structure(
          list(structure(
            c(
              778184.42148303,
              804679.592359701,
              830438.786267577,
              823814.993548409,
              778184.42148303,
              778184.42148303,
              12707401.7271398,
              12712553.5659214,
              12691946.2107951,
              12656619.3162928,
              12670602.8787,
              12707401.7271398
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
            xmin = 778184.42148303,
            ymin = 12656619.3162928,
            xmax = 830438.786267577,
            ymax = 12712553.5659214
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
  )}

poly3 <- {
  structure(
    list(
      geometry = structure(
        list(structure(
          list(structure(
            c(
              770242.699723159,
              825231.135383378,
              763577.434794648,
              770242.699723159,
              12671265.8150613,
              12652936.3365079,
              12645715.6328353,
              12671265.8150613
            ),
            .Dim = c(4L,
                     2L)
          )), class = c("XY", "POLYGON", "sfg")
        )),
        class = c("sfc_POLYGON",
                  "sfc"),
        precision = 0,
        bbox = structure(
          c(
            xmin = 763577.434794648,
            ymin = 12645715.6328353,
            xmax = 825231.135383378,
            ymax = 12671265.8150613
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
}

twoRange2 <- rbind(poly1, poly2) %>% mutate(Range = c("Missisa", "Nipigon")) %>% 
  st_set_crs(st_crs(landCoverD2))

threeRange <- rbind(poly1, poly2, poly3) %>% 
  mutate(Range = c("Missisa", "Nipigon", "Pagwachuan")) %>% 
  st_set_crs(st_crs(landCoverD2))
# supply polygon with multiple ranges

eskerD2 <- eskerDras %>%
  terra::merge(terra::shift(eskerDras, 
                              dx = terra::nrow(eskerDras)*
                                terra::xres(eskerDras), 
                              dy = terra::nrow(eskerDras)*
                                terra::xres(eskerDras))) 
linFeatD2 <- linFeatDras %>% 
  terra::merge(terra::shift(linFeatDras, 
                              dx = terra::nrow(linFeatDras)*
                                terra::xres(linFeatDras), 
                              dy = terra::nrow(linFeatDras)*
                                terra::xres(linFeatDras))) 


# same coefficients as range
resTwoRange <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras,  natDist = natDistD,
  anthroDist = anthroDistD, linFeat = linFeatDras, projectPoly = twoRange,
  caribouRange = data.frame(Range = c("Missisa", "Nipigon"), 
                            coefRange = c("Missisa", "Nipigon"), 
                            stringsAsFactors = FALSE), 
  winArea = 500
)

# different coefficients as range
resTwoRangeDif <- caribouHabitat(
  landCover = landCoverD, esker = eskerDras, natDist = natDistD,
  anthroDist = anthroDistD, linFeat = linFeatDras, projectPoly = twoRange,
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
  ext1 <- terra::ext(resRange1@habitatUse) %>% as.vector()
  pointCompare <- st_sf(ID = 1, 
                        geometry = st_sfc(st_point(c((ext1["xmax"] - ext1["xmin"])/2 + ext1["xmin"], 
                                                     (ext1["ymax"] - ext1["ymin"])/2 + ext1["ymin"]))),
                        crs = st_crs(resRange1@habitatUse))
  expect_lt(
  terra::extract(resRange1@habitatUse$Fall, pointCompare)[1,1]-
    terra::extract(resTwoRangeDifWin@habitatUse$Fall, pointCompare)[1,1],
  0.01
  )
  
  ext2 <- terra::ext(resRange2@habitatUse) %>% as.vector()
  pointCompare <- st_sf(ID = 1, 
                        geometry = st_sfc(st_point(c((ext2["xmax"] - ext2["xmin"])/2 + ext2["xmin"], 
                                                     (ext2["ymax"] - ext2["ymin"])/2 + ext2["ymin"]))),
                        crs = st_crs(resRange2@habitatUse))
  expect_lt(
    terra::extract(resRange2@habitatUse$Fall, pointCompare)[1,1]-
      terra::extract(resTwoRangeDifWin@habitatUse$Fall, pointCompare)[1,1],
    0.01
  )
  

})

# compare plots
if(interactive()){
  plot(resTwoRangeDifWin@processedData$CON)
  plot(resRange1@processedData$CON, add = T)
  plot(resRange2@processedData$CON, add = T)
}

# Try with three ranges so both apply with same window and with different
resThreeRange <- caribouHabitat(
  landCover = landCoverD2, esker = eskerD2, 
  linFeat = linFeatD2, projectPoly = twoRange2,
  caribouRange = data.frame(Range = c("Missisa", "Nipigon", "Pagwachuan"), 
                            coefRange = c("Missisa", "Nipigon", "Pagwachuan"),
                            stringsAsFactors = FALSE), 
  coefTable = coefTableSmall
)

