# Explore multirange
pth_base <- "../ChurchillAnalysis/inputNV/"
devtools::load_all(".")

caribouRanges <- c("Pagwachuan", "Missisa", "Ozhiski")
caribouRangeCoefs <- rev(caribouRanges)

projectPoly <- st_read(paste0(pth_base,
                              "caribouRanges/Caribou_Range_Boundary.shp")) %>% 
  filter(RANGE_NAME %in% caribouRanges) %>% 
  rename(Range = RANGE_NAME)

landCover <- raster(paste0(pth_base, "ROFData/plc50.tif")) %>% 
  reclassPLC() %>% 
  raster::aggregate(5, fun = raster::modal)

esker <- st_read(paste0(pth_base, "ROFData/esker.shp"))

linFeat <- list(roads = st_read(paste0(pth_base, 
                                       "ROFData/road_ORNMNRFROF2010.shp")),
                rail = st_read(paste0(pth_base, "ROFData/rail.shp")),
                utilities = st_read(paste0(pth_base, "ROFData/util2010.shp")))

caribouRange <- data.frame(Range = caribouRanges, 
                           coefRange = caribouRanges, stringsAsFactors = FALSE)

sameCoefMultRange <- caribouHabitat(landCover = landCover,
                                    esker = esker, 
                                    linFeat = linFeat,  
                                    projectPoly = projectPoly, 
                                    caribouRange = caribouRange, 
                                    padProjPoly = TRUE)
plot(sameCoefMultRange)

# switch around the coefficients 
caribouRange <- data.frame(Range = caribouRanges, 
                           coefRange = sort(caribouRanges), 
                           stringsAsFactors = FALSE)

# can just change the caribouRange attribute rather than recomputing all the
# inputs. If we think this will be done often could do it inside updateCaribou
sameCoefMultRange@attributes$caribouRange <- caribouRange

swapCoefMultRange <- updateCaribou(sameCoefMultRange)

plot(swapCoefMultRange)

# Use just the Missisa coefs for all
sameCoefMultRange@attributes$caribouRange$coefRange <- "Missisa"

misCoefMultRange <- updateCaribou(sameCoefMultRange)

plot(misCoefMultRange)

# Try Nipigon (should error because different winArea)
caribouRange <- data.frame(Range = caribouRanges, 
                           coefRange = c("Missisa", "Nipigon", "Missisa"), 
                           stringsAsFactors = FALSE)

nipCoefMultRange <- caribouHabitat(landCover = landCover,
                                    esker = esker, 
                                    linFeat = linFeat,  
                                    projectPoly = projectPoly, 
                                    caribouRange = caribouRange, 
                                    padProjPoly = TRUE)
 
