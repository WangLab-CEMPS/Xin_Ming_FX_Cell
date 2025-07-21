library(Seurat)
library(tidyr)
library(dplyr)

obj <- readRDS("./dataLib/maize/maize.rds")

markers <- FindAllMarkers(
  object = obj,
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 1
)

table(markers$cluster)

anno <- read.csv("./dataLib/maize/genes_all.csv")
anno <- anno[!duplicated(anno[, 1:2]), ]

markers <- left_join(markers, anno, by = c("gene" = "v5GeneModelID"))
head(markers)

markers <- markers %>%
  nest(data = everything(), .by = cluster)

names(markers$data) <- paste0("C", markers$cluster)

openxlsx::write.xlsx(markers$data, "./dataLib/maize/maize_markers.xlsx")
