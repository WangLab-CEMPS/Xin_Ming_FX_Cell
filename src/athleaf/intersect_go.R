library(org.At.tair.db)
library(clusterProfiler)
library(pigRutils)
library(ggplot2)
library(dplyr)
library(readxl)

# read data
up_sheet <- excel_sheets("./dataLib/athLeaf/LC_LPZ_wound_distinct_up.xlsx")
unique_up <- purrr::map(
  up_sheet,
  ~ read_excel("./dataLib/athLeaf/LC_LPZ_wound_distinct_up.xlsx", sheet = ., col_names = FALSE) %>% pull(1)
)
names(unique_up) <- up_sheet

down_sheet <- excel_sheets("./dataLib/athLeaf/LC_LPZ_wound_distinct_down.xlsx")
unique_down <- purrr::map(
  down_sheet,
  ~ read_excel("./dataLib/athLeaf/LC_LPZ_wound_distinct_down.xlsx", sheet = ., col_names = FALSE) %>% pull(1)
)
names(unique_down) <- down_sheet

# go
orgdb <- org.At.tair.db

res_up <- batch_GO(unique_up, orgdb = orgdb, keyType = "TAIR")
export_GO(res_up, file = "./dataLib/athLeaf/go/LC_LPZ_wound_distinct_up_go.xlsx")

res_down <- batch_GO(unique_down, orgdb = orgdb, keyType = "TAIR")
export_GO(res_down, file = "./dataLib/athLeaf/go/LC_LPZ_wound_distinct_down_go.xlsx")
