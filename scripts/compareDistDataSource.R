# Compare caribou habitat results for different disturbance data sources #======
devtools::load_all()
library(tmap)
library(purrr)

# if it has not already been run 
#source("scripts/RempelHabitatDefs/prepareDistDataSources.R")
# otherwise just load the saved outputs

# load data
plcD <- raster("./inputNV/intermediaryData/plc_aligned.tif")
eskerD <- raster("inputNV/intermediaryData/eskerRastDen.tif")
linFeatD <- raster("./inputNV/intermediaryData/linFeatRastDen.tif")
projectPolyD <- st_read("./inputNV/caribouRanges/Caribou_Range_Boundary.shp", 
                        quiet = TRUE) %>% 
  filter(RANGE_NAME == "Churchill") %>% 
  st_transform(st_crs(plcD))
sSimData <- getSyncSimData(ssimLib = "./inputNV/Libraries/ChurchillBC",
                           ssimProject = "2008/09 Plans",
                           rfuLUPath = "inputTables/OLTBorealCaribouFURegionalForestUnitMapping.csv",
                           scnNum = 6)

plcD <- raster::mask(plcD, projectPolyD, updatevalue = NA)




pthBase <- "inputNV/intermediaryData/"
distDataSetNms <- c("OFRID", "MNRFAFFES", "MNRFSRB", "NRCAN", "bead") 
# Run comparison for one year #=================================================
# Set year to run comparison for
yr <- 2020

distDataSets <- map(distDataSetNms %>% 
                      paste0(".*", yr), 
                    ~list.files(pthBase, pattern = .x, full.names = TRUE)) %>% 
  unlist() %>% 
  map(raster) %>% 
  map(~raster::mask(.x, projectPolyD, updatevalue = NA))

distDataSets <- map(distDataSetNms, 
                     ~distDataSets[stringr::str_which(map(distDataSets, names),
                                                      .x)]) %>% 
  set_names(distDataSetNms)

distDataSets <- distDataSets[which(map_lgl(distDataSets, ~length(.x)>0))]

# Need to change age so that new disturbances aren't updated to FRI because of
# being 35 years or older
ageNew <- map(distDataSets,
               ~raster::mask(sSimData$myAge, (.x[[1]]+.x[[2]]) > 0, maskvalue = 1,
                             updatevalue = 0))

# run model
distResults <- map2(distDataSets, ageNew,
                   ~caribouHabitat(plc = plcD, esker = eskerD, fri = sSimData$mySC, 
                                   age = .y, natDist = .x[[1]], 
                                   anthroDist = .x[[2]],
                                   harv = .x[[2]],
                                   linFeat = linFeatD, 
                                   projectPoly = projectPolyD, 
                                   friLU = sSimData$friLU,
                                   caribouRange = "Churchill"))

beepr::beep()

# Map comparison
distHabUse <- map(distResults, ~.x@habitatUse) %>% 
  map(~raster::mask(.x, projectPolyD, updatevalue = NA))

compMap <- function(x, y, nmx, nmy, rstStyl = "pretty", mdpnt = NULL){
  qtm(x - y, title = paste0(nmx, " - ", nmy), raster.style = rstStyl, 
                            midpoint = mdpnt) 
} 
nDS <- length(distHabUse) 
compBEADMaps <- map2(distHabUse[1:(nDS-1)], 
                     names(distResults)[1:(nDS-1)],
                    ~compMap(.x, distHabUse[[nDS]], .y, 
                             names(distHabUse)[nDS])) 

compExpBEADMaps <- map2(distDataSets[1:(nDS-1)], names(distDataSets)[1:(nDS-1)],
                        ~compMap(raster::stack(.x), 
                                 raster::stack(distDataSets[[nDS]]),
                                 .y, names(distHabUse)[nDS]))


compHabMaps <- map2(distHabUse, names(distResults) %>% toupper(),
                    ~qtm(.x, title = .y)) %>% 
  tmap_arrange(ncol = 1)


# explanatory variables harv and natDist (fire)
expCompMaps <- map(distDataSets,
                   ~qtm(raster::stack(.x))+ tm_layout(legend.show = FALSE)) %>% 
  tmap_arrange(ncol = 1)

# Graphs
library(ggplot2)


compGraph <- function(x, y, xnm, ynm){
  data.frame(x = x %>% raster::values(),
             y = y %>% raster::values(),
             ID = 1:raster::ncell(x)) %>%
    filter(!is.na(x.Spring)) %>% 
    gather(lab, prUSE, -ID) %>% 
    separate(lab, into = c("dataSource", "Season"), sep = "\\.") %>% 
    spread(dataSource, prUSE) %>% 
    ggplot(aes(x, y))+
    geom_hex(aes(fill = log(..count..)))+
    geom_abline(slope = 1, intercept = 0, col = "black")+
    geom_smooth(method = "lm", col = "red")+
    facet_wrap(~Season)+
    theme_bw()+
    labs(x = xnm, y = ynm)
}
  
  
compBEADGrphs <- map2(distHabUse[1:(nDS-1)], names(distResults)[1:(nDS-1)], 
                      ~compGraph(distHabUse[[nDS]], .x, "BEAD", .y))

# Coefficients graph
coefGrph <- coefTableHR %>% filter(Range == "Churchill") %>% 
  ggplot(aes(Variable, Coefficient, fill = Season))+
  geom_col(position = "dodge")+
  coord_flip()

# Simplified comparison

# category 2 habitat
distBinUse <- map(distResults, ~calcBinaryUse(.x, caribouRange = "Churchill")) %>% 
  map(~raster::mask(.x, projectPolyD, updatevalue = NA))

# area cat2
distAreaCat2 <- map(distBinUse, ~raster::cellStats(.x, sum)*raster::res(.x)[1]^2)

distAreaCat2 <- tibble(dataSource = names(distAreaCat2) %>% toupper(), 
                           areaCat2 = (unlist(distAreaCat2)/10000)/1000000)

sumCat2Plot <- ggplot(distAreaCat2, aes(dataSource, areaCat2))+
  geom_col()+
  theme_bw()+
  labs(x = "Data Source", y = "Category 2 Habitat Area (million ha)")
  
distBinUseAll <- raster::stack(distBinUse) %>% sum()

sumCat2Map <- tmapCatRast(distBinUseAll, lbls = 0:5 %>% as.character(), 
                          ttle = "", pal = "YlOrBr")

maxDif <- distAreaCat2 %>% mutate(dif = max(areaCat2) - areaCat2) %>% 
  filter(dif == max(dif)) %>% pull(dif) %>% round(3)
  
rmarkdown::render("scripts/RempelHabitatDefs/compareDistDataSource.Rmd",
                  output_file = paste0("../../outputs/compareDistDataSource",
                                       yr, ".html"))
beepr::beep()

# Compare 2010 and 2020 #======================================================

distDataSets2010 <- map(distDataSetNms %>% 
                      paste0(".*", 2010), 
                    ~list.files(pthBase, pattern = .x, full.names = TRUE)) %>% 
  unlist() %>% 
  map(raster) %>% 
  map(~raster::mask(.x, projectPolyD, updatevalue = NA))

distDataSets2010 <- map(distDataSetNms, 
                    ~distDataSets2010[stringr::str_which(map(distDataSets2010, names),
                                                     .x)]) %>% 
  set_names(distDataSetNms)

distDataSets2010 <- distDataSets2010[which(map_lgl(distDataSets2010, ~length(.x)>0))]

age2010 <- map(distDataSets2010,
               ~raster::mask(sSimData$myAge, (.x[[1]]+.x[[2]]) > 0, maskvalue = 1,
                             updatevalue = 0))

# run model
distResults2010 <- map2(distDataSets2010, age2010,
                   ~caribouHabitat(plc = plcD, esker = eskerD, fri = sSimData$mySC, 
                                   age = .y, natDist = .x[[1]], 
                                   anthroDist = .x[[2]],
                                   harv = .x[[2]],
                                   linFeat = linFeatD, 
                                   projectPoly = projectPolyD, 
                                   friLU = sSimData$friLU,
                                   caribouRange = "Churchill"))

beepr::beep()

distDataSets2020 <- map(distDataSetNms %>% 
                      paste0(".*", 2020), 
                    ~list.files(pthBase, pattern = .x, full.names = TRUE)) %>% 
  unlist() %>% 
  map(raster) %>% 
  map(~raster::mask(.x, projectPolyD, updatevalue = NA))

distDataSets2020 <- map(distDataSetNms, 
                    ~distDataSets2020[stringr::str_which(map(distDataSets2020, names),
                                                     .x)]) %>% 
  set_names(distDataSetNms)

distDataSets2020 <- distDataSets2020[which(map_lgl(distDataSets2020, ~length(.x)>0))]

age2020 <- map(distDataSets2020,
               ~raster::mask(sSimData$myAge, (.x[[1]]+.x[[2]]) > 0, maskvalue = 1,
                             updatevalue = 0))

# run model
distResults2020 <- map2(distDataSets2020, age2020,
                   ~caribouHabitat(plc = plcD, esker = eskerD, fri = sSimData$mySC, 
                                   age = .y, natDist = .x[[1]], 
                                   anthroDist = .x[[2]],
                                   harv = .x[[2]],
                                   linFeat = linFeatD, 
                                   projectPoly = projectPolyD, 
                                   friLU = sSimData$friLU,
                                   caribouRange = "Churchill"))

beepr::beep()

distBinUse2010 <- map(distResults2010, 
                      ~calcBinaryUse(.x, caribouRange = "Churchill", 
                                     bySeason = TRUE)) %>% 
  map(~raster::mask(.x, projectPolyD, updatevalue = NA))

distBinUse2020 <- map(distResults2020, 
                      ~calcBinaryUse(.x, caribouRange = "Churchill",
                                     bySeason = TRUE)) %>% 
  map(~raster::mask(.x, projectPolyD, updatevalue = NA))

# area cat2
distAreaCat22010 <- map(distBinUse2010, 
                    ~raster::cellStats(.x, sum)*raster::res(.x)[1]^2)

distAreaCat22020 <- map(distBinUse2020, 
                        ~raster::cellStats(.x, sum)*raster::res(.x)[1]^2)

cat2Both <- bind_rows(
  tibble(dataSource = names(unlist(distAreaCat22010)) %>% toupper(), 
         Year = 2010,
         areaCat2 = (unlist(distAreaCat22010)/10000)/1000000),
  tibble(dataSource = names(unlist(distAreaCat22020)) %>% toupper(), 
         Year = 2020,
         areaCat2 = (unlist(distAreaCat22020)/10000)/1000000)
  ) %>% 
  separate(dataSource, c("dataSource", "Season"))

pltCat2Both <- cat2Both %>% 
  ggplot(aes(dataSource, areaCat2, fill = as.factor(Year)))+
  geom_col(position = position_dodge2(preserve = "single"))+
  theme_bw()+
  labs(x = "Data Source", y = "Category 2 Habitat Area (million ha)",
       fill = "Year")+
  facet_wrap(~Season)

# differences between different data sources in the same year
cat2Both %>% filter(dataSource != "MNRFSRB") %>% group_by(Year, Season) %>% 
  summarise(maxdif = (max(areaCat2) - min(areaCat2))/min(areaCat2)) %>% 
  summarise(avgMaxDif = mean(maxdif))

cat2Both %>% group_by(dataSource, Season) %>% 
  summarise(maxdif = (max(areaCat2) - min(areaCat2))/min(areaCat2)) %>% 
  summarise(avgMaxDif = mean(maxdif))

(cat2Both %>% filter(dataSource == "BEAD", Season == "SPRING", Year == 2010) %>% 
  pull(areaCat2) - 
  cat2Both %>% filter(dataSource == "OFRID", Season == "SPRING", Year == 2010) %>% 
  pull(areaCat2))/cat2Both %>% filter(dataSource == "BEAD", Season == "SPRING", Year == 2010) %>% 
  pull(areaCat2)

# change to a binary map of Category 2 in any season which is how Rempel does it
distBinUse2010Any <- map(distBinUse2010, ~any(.x == 1))

distBinUse2020Any <- map(distBinUse2020, ~any(.x == 1))

# Compare different years same data source which year
mapCat2Each <- map2(distBinUse2010Any[1:4], distBinUse2020Any,
                    ~raster::stack(.x, .y) %>% sum()) %>% 
  map2(distBinUse2010Any[1:4], ~raster::mask(.x, (.x + .y) == 2, maskvalue = 1, 
                                             updatevalue = 2010)) %>% 
  map(~raster::mask(.x, .x, maskvalue = 1, updatevalue = 2020)) %>% 
   map2(names(distBinUse2020),
       ~tm_shape(.x)+
         tm_raster(style = "cat", 
                   palette = tmaptools::get_brewer_pal("Paired", n = 8, plot = FALSE)[c(3,4,7,8)], 
                   labels = c("Neither Cat2", "Both Cat2", "2010 Cat2", 
                              "2020 Cat2"), title = .y)+
         tm_legend(bg.color = "white", legend.outside = FALSE, 
                   legend.position = c("right", "bottom"))) %>% 
  tmap_arrange(ncol = 2)

# Compare same year different data source 
sumCat2Map2010 <- tmapCatRast(raster::stack(distBinUse2010Any) %>% sum(), 
                              lbls = 0:5 %>% as.character(), 
                              ttle = "2010", pal = "YlOrBr")+
  tm_layout(legend.outside = FALSE, 
            legend.position = c("right", "bottom"))

sumCat2Map2020 <- tmapCatRast(raster::stack(distBinUse2020Any) %>% sum(), 
                              lbls = 0:5 %>% as.character(), 
                              ttle = "2020", pal = "YlOrBr")+
  tm_layout(legend.outside = FALSE, 
            legend.position = c("right", "bottom"))

# total area the same or different
cellsAgree2010 <- raster::freq(raster::stack(distBinUse2010Any) %>% sum(),
                               useNA = "no") %>%
  as.data.frame() 

cellsAgree2020 <- raster::freq(raster::stack(distBinUse2020Any) %>% sum(),
                               useNA = "no") %>%
  as.data.frame() 

areaTotal <- sum(cellsAgree2020$count) * 400^2 / 10000
areaAgree2020 <- (cellsAgree2020[1,2] + cellsAgree2020[5,2]) * 400^2 / 10000

areaAgree2020/areaTotal

areaAgree2010 <- (cellsAgree2010[1,2] + cellsAgree2010[6,2]) * 400^2 / 10000

areaAgree2010/areaTotal

tmap_save(tmap_arrange(sumCat2Map2010, sumCat2Map2020, nrow = 1),
          "outputs/distCat2NumMap2010and2020.png", width = 7, height = 4)


ggsave("outputs/distCat2Area2010vs2010.png", pltCat2Both, width = 9.11,
       height = 4.38)
tmap_save(mapCat2Each, "outputs/distCat2Maps2010vs2020.png", width = 7, height = 4)

layout1 <- tm_layout(legend.show = FALSE,
                     title.position = c("LEFT", "TOP"), 
                     title.size = 1)

map2(distDataSets2010[1:4], distDataSets2020, 
     ~tmap_arrange(qtm(.x[[1]], title = names(.x[[1]]))+layout1, 
                   qtm(.x[[2]], title = names(.x[[2]]))+layout1, 
                   qtm(.y[[1]], title = names(.y[[1]]))+layout1, 
                   qtm(.y[[2]], title = names(.y[[2]]))+layout1,
                   ncol = 2)) %>% 
  map2(names(distDataSets2020), 
       ~tmap_save(.x, paste0("outputs/fireHarv2010vs2020", .y, ".png"),
                             width = 7, height = 4))

tmap_arrange(qtm(distDataSets2010[[5]][[1]], title = "BEAD_fire2010")+ layout1, 
             qtm(distDataSets2010[[5]][[2]], title = "BEAD_harv2010") + layout1, 
             ncol = 1) %>% 
  tmap_save("outputs/fireHarv2010BEAD.png", width = 3.5, height = 4)
