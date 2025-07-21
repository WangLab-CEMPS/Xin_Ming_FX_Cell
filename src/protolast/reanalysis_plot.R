library(Seurat)
library(ggplot2)
library(patchwork)


obj <- readRDS("./dataLib/athLeaf/LC_public.Rds")
obj_normal <- readRDS("./dataLib/athLeaf/LC_public_harmony.Rds")
obj_protolast <- readRDS("./dataLib/athLeaf/LC_public_remove_protolasted_harmony.Rds")

cc <- read.csv("~/ref/Arabidopsis_thaliana/Arabidopsis_cell_cycle_genes_fx_update.csv")
table(cc$phase)

obj_normal <- CellCycleScoring(
  obj_normal,
  s.features = cc$gene[cc$phase == "S"],
  g2m.features = cc$gene[cc$phase == "G2M"]
)

obj_protolast <- CellCycleScoring(
  obj_protolast,
  s.features = cc$gene[cc$phase == "S"],
  g2m.features = cc$gene[cc$phase == "G2M"]
)


p1 <- DimPlot(
  obj,
  pt.size = 1.7,
  group.by = "orig.ident", reduction = "umap",
  raster = TRUE
) + ggtitle("Normal pipeline without harmony")

p2 <- DimPlot(
  obj_normal,
  pt.size = 1.7,
  group.by = "orig.ident", reduction = "umap",
  raster = TRUE
) + ggtitle("Normal pipeline with harmony")

p3 <- DimPlot(
  obj_protolast,
  pt.size = 1.7,
  group.by = "orig.ident",
  reduction = "umap",
  raster = TRUE
) + ggtitle("harmony after regressing out Protolasted genes and remove Protolasted genes in VariableFeatures")


pdf("./Plots/protolast_reanalysis.pdf", width = 29, height = 9.4)
(p1 + p2 + p3) & NoAxes()
dev.off()
