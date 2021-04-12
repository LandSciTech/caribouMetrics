# Explore multirange
pth_base <- "../ChurchillAnalysis/inputNV/"
devtools::load_all(".")

caribouRanges <- c("Pagwachuan", "Missisa", "Ozhiski", "Nipigon", "James Bay")
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
beepr::beep()
plot(sameCoefMultRange)

# switch around the coefficients 
caribouRange <- data.frame(Range = caribouRanges, 
                           coefRange = sort(caribouRanges), 
                           stringsAsFactors = FALSE)

# can just change the caribouRange attribute rather than recomputing all the
# inputs. This should be done by the updateCaribou function so that it considers
# changes in winArea that mean the whole area should be recalculated


swapCoefMultRange <- updateCaribou(sameCoefMultRange, caribouRange)

plot(swapCoefMultRange)

# example of just changing the attribute, **NOT RECOMMENDED**
# Use just the Missisa coefs for all
sameCoefMultRange@attributes$caribouRange$coefRange <- "Missisa"

misCoefMultRange <- updateCaribou(sameCoefMultRange)

plot(misCoefMultRange)

