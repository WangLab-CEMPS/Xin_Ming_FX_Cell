library(Seurat)

obj <- readRDS("./dataLib/athRoot/bigRoots_FX.Rds")
obj <- UpdateSeuratObject(obj)
saveRDS(obj, "./dataLib/athRoot/bigRoots_FX_v5.rds")

DimPlot(obj)
