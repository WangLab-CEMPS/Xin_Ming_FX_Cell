library(Seurat)
library(tidyr)
library(dplyr)

obj <- readRDS("./dataLib/tillering_node_rhizome/tillering.rds")

markers <- FindAllMarkers(
  object = obj, only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 1.32
)
table(markers$cluster)

anno <- read.csv("~/ref/Oryza_sativa/annotation/IRGSP_OsAt_Orthologous.tsv", sep = "\t")
anno <- distinct(anno, RAP_ID, .keep_all = TRUE)
markers <- left_join(markers, anno, by = c("gene" = "RAP_ID"))
head(markers)

go <- markers %>%
  nest(data = everything(), .by = cluster) %>%
  rowwise() %>%
  mutate(tair = list(data$TAIR))

names(go$tair) <- paste0("Cluster", go$cluster)
names(go$data) <- paste0("Cluster", go$cluster)

openxlsx::write.xlsx(go$data, "./dataLib/tillering_node_rhizome/tillering_markers.xlsx")

# GO ---------
library(pigRutils)
library(org.At.tair.db)
library(clusterProfiler)

orgdb <- org.At.tair.db

go_res <- batch_GO(go$tair, orgdb = orgdb, keyType = "TAIR")

export_GO(go_res, "./dataLib/tillering_node_rhizome/tillering_markers_GO.xlsx")
