library(Seurat)
library(magrittr)
library(ggplot2)

fp <- c("./dataLib/athLeaf/LC.rds", "./dataLib/athLeaf/LPZ.rds")

for (p in fp) {
  obj <- readRDS(p)
  obj %<>%
    NormalizeData() %>%
    FindVariableFeatures(
      selection.method = "vst",
      nfeatures = round(0.1 * length(rownames(.)))
    ) %>%
    ScaleData(features = rownames(.), var.to.regress = "nCount_RNA") %>%
    RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)

  dims <- 1:50
  obj <- FindNeighbors(obj, dims = dims)

  obj %<>%
    FindClusters(resolution = 0.7) %>%
    RunTSNE(dims = dims) %>%
    RunUMAP(dims = dims, metric = "correlation")

  saveRDS(obj, p)
}
