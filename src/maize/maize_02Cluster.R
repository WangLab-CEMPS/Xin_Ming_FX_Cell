library(Seurat)
library(ggplot2)
library(magrittr)
library(reticulate)


obj <- readRDS("./dataLib/maize/maize.rds")
obj


# ----
obj %<>%
  NormalizeData() %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = 2000
  ) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)

ElbowPlot(obj, reduction = "pca", ndims = 100)


dims <- 1:20
obj <- FindNeighbors(obj, dims = dims)

obj %<>%
  FindClusters(resolution = 0.7) %>%
  RunUMAP(dims = dims, metric = "correlation", n.neighbors = 30, spread = 2)



# -----------
p_umap <- DimPlot(
  obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE, pt.size = .6,
  cols = pigRutils::select_colors("col29")
) &
  tidydr::theme_dr() &
  theme(
    panel.grid = element_blank(),
    legend.position = "top", title = element_blank()
  ) &
  guides(
    color = guide_legend(
      keywidth = 0.4,
      nrow = 2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

p_umap

ggsave(filename = "./Plots/maize/maize_umap.pdf", p_umap, width = 7, height = 7.5)

saveRDS(obj, "./dataLib/maize/maize.rds")


# ------------------
anno <- read.delim("./dataLib/maize/genes_all.txt", sep = "\t", fill = TRUE, header = TRUE)
mm <- read.table("./dataLib/maize/marker.txt")[[1]]

anno <- anno[anno$v5.Gene.Model.ID %in% mm, ]
anno <- anno[, 1:2] %>% unique()

mm <- mm[mm %in% rownames(obj)]

for (i in mm) {
  p <- FeaturePlot(obj, i)
  png(paste0("./tmp/", i, ".png"), height = 400, width = 410)
  print(p)
  dev.off()
}

DotPlot(obj, features = mm)

# Endodermis      Zm00001eb083140 Zm00001eb422410 Zm00001eb372230 Zm00001eb249210 Zm00001eb333290
# Endo early      Zm00001eb408720
# Initials        Zm00001eb370670 Zm00001eb035350 Zm00001eb188180
# Phloem & Sieve  Zm00001eb266840 Zm00001eb002390 Zm00001eb173060 Zm00001eb315650 Zm00001eb342070
# Pericycle       Zm00001eb356620
# Pericycle       Zm00001eb397500 Zm00001eb302250
# Stele           Zm00001eb175860
# Xylem           Zm00001eb267050
## Stele 7        Zm00001eb338510
# Cortex          Zm00001eb000450 Zm00001eb330530 Zm00001eb249760 Zm00001eb429560
# Cortex Initials Zm00001eb069630
# CortexII(M)     Zm00001eb187430
# Endodermal Initials Zm00001eb388430
## Epidermis      Zm00001eb255380
# QC              Zm00001eb200450
# Root cap        Zm00001eb091140 Zm00001eb140000 Zm00001eb091170

obj[["celltype"]] <- obj$seurat_clusters
obj[["celltype"]] <- as.character(obj$celltype)

obj$celltype[obj$celltype %in% c("6", "18")] <- "Endodermis"
obj$celltype[obj$celltype %in% c("11")] <- "Initials CC"
obj$celltype[obj$celltype %in% c("13", "14", "2", "15")] <- "Phloem & Sieve"
obj$celltype[obj$celltype %in% c("8")] <- "Pericycle"
obj$celltype[obj$celltype %in% c("0", "10")] <- "Stele(M)"
obj$celltype[obj$celltype %in% c("7")] <- "Stele(E)"
obj$celltype[obj$celltype %in% c("9")] <- "Xylem"
obj$celltype[obj$celltype %in% c("1", "5", "4")] <- "Cortex"
obj$celltype[obj$celltype %in% c("16")] <- "Cortex Initials"
obj$celltype[obj$celltype %in% c("3")] <- "CortexII(m)"
obj$celltype[obj$celltype %in% c("19", "17")] <- "Root Cap"
obj$celltype[obj$celltype %in% c("12")] <- "QC"


DimPlot(obj, group.by = "celltype")

saveRDS(obj, "./dataLib/maize/maize.rds")
