library(Seurat)
library(dplyr)
library(tidyr)

obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public.Rds")
Layers(obj)
obj <- JoinLayers(obj)

markers <- FindAllMarkers(
  object = obj,
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 1.32
)

table(markers$cluster)

markers_celltype <- FindAllMarkers(
  object = obj,
  group.by = "celltype",
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 2
)

table(markers_celltype$cluster)

go <- markers %>%
  nest(data = everything(), .by = cluster) %>%
  rowwise() %>%
  mutate(tair = list(data$gene))

names(go$tair) <- paste0("Cluster", go$cluster)
names(go$data) <- paste0("Cluster", go$cluster)

go_celltype <- markers_celltype %>%
  nest(data = everything(), .by = cluster) %>%
  rowwise() %>%
  mutate(tair = list(data$gene))

names(go_celltype$tair) <- go_celltype$cluster
names(go_celltype$data) <- go_celltype$cluster

openxlsx::write.xlsx(go$data, "./dataLib/athLeaf/leaf_LC_LPZ_public_markers.xlsx")
openxlsx::write.xlsx(go_celltype$data, "./dataLib/athLeaf/leaf_LC_LPZ_public_markers_celltype.xlsx")

# GO ---------
library(pigRutils)
library(org.At.tair.db)
library(clusterProfiler)

orgdb <- org.At.tair.db

go_res <- batch_GO(go$tair, orgdb = orgdb, keyType = "TAIR")
export_GO(go_res, "./dataLib/athLeaf/leaf_LC_LPZ_public_markers_GO.xlsx")

go_res <- batch_GO(go_celltype$tair, orgdb = orgdb, keyType = "TAIR")
names(go_res) <- gsub("\\*", "_", names(go_res))
export_GO(go_res, "./dataLib/athLeaf/leaf_LC_LPZ_public_markers_celltype_GO.xlsx")
