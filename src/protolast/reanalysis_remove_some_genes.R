library(magrittr)
library(Seurat)
library(readxl)
library(dplyr)

# read in the protolasted data -----------
dp <- "./dataLib/athLeaf/PC_Wolf_B_Frommer_protoplasted.xlsx"
sheet <- excel_sheets(dp)
protolast <- read_excel(dp, sheet = sheet[4])
colnames(protolast) <- gsub(" ", "_", colnames(protolast))
protolast <- protolast %>% subset(Log2_FC > 1.58 & padj < 0.05)
dim(protolast)


# scRNAseq data ------------
fl <- list.files("./dataLib/athLeaf", full.names = TRUE, pattern = "rds$")
fl <- fl[c(1, 2)]
n <- gsub(".rds", "", basename(fl))
n[2] <- "JC"


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

obj$orig.ident[obj$orig.ident == "SeuratProject"] <- "JC"
table(obj$orig.ident)

obj <- NormalizeData(obj)

# regress out and remove the protolasted genes ------------
obj <- AddModuleScore(
  obj,
  features = list(protolast$Gene_ID[protolast$Gene_ID %in% rownames(obj)]),
  name = "protolast"
)

colnames(obj@meta.data)

obj %<>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = round(0.1 * nrow(obj))
  ) %>%
  ScaleData(vars.to.regress = "protolast1")


table(VariableFeatures(obj) %in% protolast$Gene_ID)

VariableFeatures(obj) <- VariableFeatures(obj)[!VariableFeatures(obj) %in% protolast$Gene_ID]

obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 100, verbose = FALSE)

obj <- IntegrateLayers(
  object = obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE,
  theta = 4
)

dims <- 1:15
mm <- "harmony"

obj <- FindNeighbors(obj, dims = dims, reduction = mm)

obj %<>%
  FindClusters(resolution = 0.6, cluster.name = "harmony_clusters") %>%
  RunUMAP(dims = dims, n.neighbors = 50, reduction = mm)


saveRDS(obj, "./dataLib/athLeaf/LC_public_remove_protolasted_harmony.Rds")
