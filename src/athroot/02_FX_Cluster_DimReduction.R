library(Seurat)
library(magrittr)
library(ggplot2)


fl <- c("./dataLib/athRoot/atroot1.rds", "./dataLib/athRoot/atroot2.rds")
n <- c("atroot1", "atroot2")

obj_list <- lapply(fl, function(x) {
  readRDS(x) %>%
    GetAssayData(layer = "counts") %>%
    CreateSeuratObject()
})

(obj <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(fl)],
  add.cell.ids = n
))

obj <- JoinLayers(obj)
obj

rm(obj_list)
gc()

obj %<>%
  NormalizeData() %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = round(0.1 * length(rownames(.)))
  ) %>%
  ScaleData(features = rownames(.), var.to.regress = "nCount_RNA") %>%
  RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)

ElbowPlot(obj, reduction = "pca", ndims = 100)


dims <- 1:40
obj <- FindNeighbors(obj, dims = dims)

obj %<>%
  FindClusters(resolution = 0.6) %>%
  RunTSNE(dims = dims) %>%
  RunUMAP(dims = dims)

DimPlot(
  obj,
  reduction = "umap",
  label = TRUE,
  pt.size = 1
)


saveRDS(obj, "./dataLib/athRoot/atroot12_fx.rds")
