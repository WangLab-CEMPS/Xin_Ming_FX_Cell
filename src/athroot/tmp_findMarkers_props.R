library(Seurat)
library(readxl)
library(tidyr)
library(dplyr)

obj <- readRDS("./dataLib/athRoot/FX_TQ_v5.rds")
table(obj$tech)

FX <- obj[, obj$tech == "atroot"]

FX %<>%
  NormalizeData() %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = round(0.1 * length(rownames(.)))
  ) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)


dims <- 1:40
FX <- FindNeighbors(FX, dims = dims)

FX %<>%
  FindClusters(resolution = 0.6) %>%
  RunUMAP(dims = dims)


TQ <- obj[, obj$tech == "tq"]
TQ

TQ %<>%
  NormalizeData() %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = round(0.1 * length(rownames(.)))
  ) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)


dims <- 1:40
TQ <- FindNeighbors(TQ, dims = dims)

TQ %<>%
  FindClusters(resolution = 0.6) %>%
  RunUMAP(dims = dims)



markers_FX <- FindAllMarkers(
  FX,
  group.by = "celltype",
  only.pos = TRUE
)

markers_TQ <- FindAllMarkers(
  TQ,
  group.by = "celltype",
  only.pos = TRUE
)


proto_file <- "./dataLib/athRoot/2019_DC_athroot_Tom_Denyer.xlsx"
sheets <- excel_sheets(proto_file)
sheets

proto <- read_excel(proto_file, sheet = sheets[3])
colnames(proto) <- gsub(" ", "_", colnames(proto))

proto_genes <- proto$Gene_ID[proto$Log2_FC >= 2]
length(proto_genes)

mFX <- markers_FX %>%
  group_by(cluster) %>%
  top_n(n = 1000, wt = avg_log2FC) %>%
  nest()

names(mFX$data) <- mFX$cluster

mTQ <- markers_TQ %>%
  group_by(cluster) %>%
  top_n(n = 1000, wt = avg_log2FC) %>%
  nest()

names(mTQ$data) <- mTQ$cluster


sapply(names(mFX$data), function(x) {
  c(
    FX = prop.table(table(mFX$data[[x]]$gene %in% proto_genes))[2] * 100,
    TQ = prop.table(table(mTQ$data[[x]]$gene %in% proto_genes))[2] * 100,
    samples = prop.table(table(rownames(obj)[sample(length(rownames(obj)), 1000)] %in% proto_genes))[2] * 100
  )
})


sapply(names(mFX$data), function(x) {
  c(
    FXs = prop.table(
      table(mFX$data[[x]]$gene %in% rownames(obj)[sample(length(rownames(obj)), length(proto_genes))])
    )[2] * 100,
    TQs = prop.table(
      table(mTQ$data[[x]]$gene %in% rownames(obj)[sample(length(rownames(obj)), length(proto_genes))])
    )[2] * 100
  )
})
