library(magrittr)
library(Seurat)
library(dplyr)


fl <- list.files("./dataLib/athLeaf", full.names = TRUE, pattern = "rds$")
n <- gsub(".rds", "", basename(fl))
n[2] <- "JC"
n[3] <- "WB"


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

rm(obj_list)
gc()

obj$orig.ident <- stringr::str_split(Cells(obj), "_", simplify = TRUE)[, 1]
table(obj$orig.ident)


obj %<>%
  NormalizeData() %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = round(0.1 * nrow(.))
  ) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)

obj <- IntegrateLayers(
  object = obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE,
  theta = 4
)

dims <- 1:30
mm <- "harmony"

obj <- FindNeighbors(obj, dims = dims, reduction = mm)

obj %<>%
  FindClusters(resolution = 0.6, cluster.name = "harmony_clusters") %>%
  RunUMAP(dims = dims, n.neighbors = 50, reduction = mm)

DimPlot(obj, group.by = "orig.ident")

saveRDS(obj, "./dataLib/athLeaf/LCLPZ_public_JCWB_harmony.Rds")
