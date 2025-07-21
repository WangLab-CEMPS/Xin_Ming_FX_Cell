library(Seurat)
library(dplyr)


obj <- readRDS("./dataLib/maize/maize_pub.Rds")
DimPlot(obj, label = TRUE)
obj$cc <- obj$celltype
obj$celltype[obj$seurat_clusters == 4 | obj$celltype == "Root Cap"] <- "Root Cap & Epidermis"

saveRDS(obj, "./dataLib/maize/maize_pub.Rds")
