library(dplyr)
library(Seurat)
library(ggplot2)

obj <- readRDS("./dataLib/osRoot/ggm_dgm_gdm.Rds")

obj2 <- readRDS("./dataLib/osRoot/ggm/ggm2_TQ2_strata.Rds")
obj2 <- obj2[, obj2$orig.ident == "ggm2"]

colnames(obj) <- gsub("ggm_ggm[1-2]_|gdm_gdm[1-2]_|dgm_dgm[12]_", "", colnames(obj))

mm <- obj2@meta.data %>%
  mutate(cell = rownames(.)) %>%
  select(cell, celltype)
meta <- obj@meta.data %>%
  mutate(cell = rownames(.)) %>%
  left_join(., mm, by = "cell")

rownames(meta) <- meta$cell
table(meta$celltype)
obj@meta.data <- meta

colnames(obj@meta.data)[8] <- "cc"

rm(obj2)

obj0 <- subset(obj, subset = seurat_clusters == 0)
obj0 <- FindClusters(obj0, resolution = 0.2)

levels(obj$seurat_clusters) <- c(levels(obj$seurat_clusters), "20")
obj$seurat_clusters[Cells(obj) %in% Cells(obj0)[obj0$seurat_clusters %in% c(3, 4)]] <- 20
Idents(obj) <- obj$seurat_clusters
DimPlot(obj)

saveRDS(obj, "./dataLib/osRoot/ggm_dgm_gdm.Rds")




