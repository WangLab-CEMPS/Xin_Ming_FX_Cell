library(dplyr)
library(tidyr)
library(readxl)


fp <- "./dataLib/athRoot/FX_TQ_protolasted_only_pos_TRUE.xlsx"
sheets <- excel_sheets(fp)
sheets <- sheets[1:12]

de_data <- lapply(sheets, function(x) {
  read_excel(fp, sheet = x) %>%
    filter(abs(avg_log2FC) >= 2)
})
names(de_data) <- sheets
sapply(de_data, length)
head(de_data[[1]])


new_de_data <- lapply(sheets, function(x) {
  dd <- de_data[[x]]
  others <- de_data[setdiff(sheets, x)]
  others <- others %>%
    lapply(function(x) x$gene) %>%
    unlist() %>%
    unique()

  dd %>%
    mutate(Cell_type_specific = ifelse(gene %in% others, "No", "Yes")) %>%
    mutate(Cell_type = x)
})

names(new_de_data) <- sheets

sapply(new_de_data, function(x) {
  table(x$Cell_type_specific)
})


new_de_data <- do.call(rbind, new_de_data)

openxlsx::write.xlsx(new_de_data, "./dataLib/athRoot/FX_TQ_protolasted_label.xlsx")
