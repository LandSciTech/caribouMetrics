## code to prepare `threshTable` dataset goes here
library(pdftools)
library(tidyverse)

# load pdf text
raw_text <- pdf_text("./inputNV/RempelHornseth-CaribouHabitatUseModels (002) SRB-IFR-01-.pdf")

# Get the thershold table from pdf
colNms <- c("Range", "Season", "Sensitivity",
            "Specificity", "Threshold", "AUC")

table <- str_split(raw_text[61], "\n", simplify = FALSE) %>% unlist()
table_start <- stringr::str_which(table, "Range")
table_end <- stringr::str_which(table, "1. Sensitivity")-1
table2 <- table[(table_start+1):(table_end)]

table3 <- str_replace_all(table2, "\\s{2,}", "|")

threshTable <-  read.csv(textConnection(table3), sep = "|", header = FALSE, 
                         stringsAsFactors = FALSE) %>%
  setNames(colNms)

# Some adjustments need to be made by hand
threshTable[5,1] <- paste(threshTable[5,1], threshTable[6,1], collapse = " ")
threshTable[6,1] <- NA_character_
threshTable[17,1] <- paste(threshTable[17,1], threshTable[18,1], collapse = "")
threshTable[18,1] <- NA_character_

threshTable <- threshTable %>% 
  mutate(Range = ifelse(Range == "", NA_character_, Range)) %>% 
  fill(Range, .direction = "down") %>% 
  separate_rows(Range, sep = " and|/") %>% 
  mutate(Range = Range %>% str_trim()) %>% 
  arrange(Range)

# Kinloch is messed up because it had different model for summer, fix
threshTable[22, 5] <- threshTable[25, 5]

threshTable <- threshTable %>% filter(Range != "Kinloch Summer")

usethis::use_data(threshTable)
