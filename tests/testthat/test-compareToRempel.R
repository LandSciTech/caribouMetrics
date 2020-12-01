# context("Compare results to Rempel's results")  
# 
# pthBase <- "../../"
# #pthBase <- "" 
# 
# # Our results
# sSimData <- getSyncSimData(ssimLib = paste0(pthBase, "inputNV/Libraries/ChurchillBC"),
#                            ssimProject = "2008/09 Plans",
#                            rfuLUPath = paste0(pthBase, "inputTables/OLTBorealCaribouFURegionalForestUnitMapping.csv"),
#                            scnNum = 6)
# 
# # data as rasters 
# plcD <- raster(paste0(pthBase,"inputNV/intermediaryData/plc_aligned.tif"))
# #eskerD <- st_read(paste0(pthBase, "inputNV/HornsethRempelInfo/Eskers_Ontario/Eskers_Ontario.shp"))
# eskerD <- raster(paste0(pthBase, "inputNV/intermediaryData/eskerRastDen.tif"))
# natDistD <- raster(paste0(pthBase, "inputNV/intermediaryData/natDistOFRI.tif"))
# # dummy anthroDist based on linFeats > 0 
# anthroDistD <- raster(paste0(pthBase, "inputNV/intermediaryData/dummyAnthroDist.tif"))
# harvD <- raster(paste0(pthBase, "inputNV/intermediaryData/harvPost2000OFRI.tif"))
# linFeatD <- raster(paste0(pthBase, "inputNV/intermediaryData/linFeatRastDen.tif"))
# #linFeatD <- st_read(paste0(pthBase, "inputNV/Roads/TroutLakeRoads_ornsegad.shp"))
# projectPolyD <- st_read(paste0(pthBase, "inputNV/caribouRanges/Caribou_Range_Boundary.shp"), quiet = TRUE) %>% 
#   filter(RANGE_NAME == "Churchill")
# 
# caribouResults <- caribouHabitat(plc = plcD, esker = eskerD, fri = sSimData$mySC, 
#                                  age = sSimData$myAge, natDist = natDistD, 
#                                  anthroDist = anthroDistD,
#                                  harv = harvD,
#                                  linFeat = linFeatD, projectPoly = projectPolyD, 
#                                  friLU = sSimData$friLU, 
#                                  winArea = 5000,
#                                  caribouRange = "Churchill")
# 
# # Rempel's results
# rempelResults592 <- st_read(paste0(pthBase, "inputNV/HornsethRempelInfo", 
#                                    "/Churchill_LSL_Files/Churchill_4_592_1.shp"), quiet = TRUE) %>% 
#   st_transform( crs = st_crs(caribouResults@plc))
# 
# rempelResults592 <- rempelResults592 %>% mutate(LGMD_S5 = DEC_S5+MIX_S5)
# 
# # Summarise to hexs using the rempel Hexs
# caribouResult592 <- raster::extract(raster::stack(caribouResults@habitatUse, 
#                                                   caribouResults@processedData),
#                                     rempelResults592 %>% select(HEXID), 
#                                     fun = mean, df = TRUE, sp = TRUE) %>% 
#   st_as_sf() 
# 
# 
# compar <- full_join(rempelResults592 %>% 
#                       st_drop_geometry() %>% 
#                       select(HEXID, contains("_S5"), matches("P.._USE")), 
#                     caribouResult592 %>% st_drop_geometry(), by = "HEXID", 
#                     suffix = c("_remp", "endi")) %>% 
#   filter(!is.na(CON)) %>% 
#   select(-CONST, -other)
# 
# comparLong <- compar %>% 
#   select(HEXID, matches("P.._USE"), Spring, Summer, Fall, Winter) %>% 
#   gather(season, Pred, -HEXID) %>% 
#   mutate(Source = ifelse(stringr::str_detect(season, "P.._USE"), 
#                          "Rempel", "Endicott"),
#          season = case_when(season == "PFA_USE" ~ "Fall",
#                             season == "PSP_USE" ~ "Spring",
#                             season == "PSU_USE" ~ "Summer",
#                             season == "PWI_USE" ~ "Winter",
#                             TRUE ~ season)) %>% 
#   spread(Source, Pred) 
# 
# comparLM <- comparLong %>% group_by(season) %>% nest() %>% 
#   mutate(linMod  = map(data, ~lm(Rempel ~ Endicott, data = .x)), 
#          r.sq = map_dbl(linMod, ~broom::glance(.x)$r.squared), 
#          slope = map_dbl(linMod, ~broom::tidy(.x)$estimate[2]))
# 
# test_that("slope and rsq within 0.25 of 1 and r ",{
#   expect_true(max(abs(comparLM$slope-1)) < 0.30)
#   expect_true(max(abs(comparLM$r.sq-1)) < 0.30)
# })


