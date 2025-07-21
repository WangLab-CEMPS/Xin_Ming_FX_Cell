library(Seurat)
library(dplyr)
library(ComplexHeatmap)

# read in RDS
obj <- readRDS("./dataLib/osRoot/ggm_dgm_gdm.Rds")
colnames(obj@meta.data)
table(obj$orig.ident)
table(obj$celltype)

# pseudo-bulk by tech and celltype
Assays(obj)

bulk <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = c("tech", "celltype"),
  return.seurat = TRUE
)
counts <- GetAssayData(bulk, assay = "RNA", layer = "counts")
dim(counts)

# get the correlation matrix
cor_matrix <- cor(as.matrix(counts))

col_fun <- circlize::colorRamp2(
  c(min(cor_matrix), mean(cor_matrix), max(cor_matrix)),
  c("#268c3a", "#ffffff", "#db1474")
)

gg <- factor(
  stringr::str_split(colnames(cor_matrix), "_", simplify = TRUE)[, 2]
)

od <- order(stringr::str_split(colnames(cor_matrix), "_", simplify = TRUE)[, 2])
od2 <- order(stringr::str_split(rownames(cor_matrix), "_", simplify = TRUE)[, 2])
pd1 <- cor_matrix[od2, od]

p1 <- Heatmap(
  pd1,
  name = "Correlation",
  col = col_fun,
  row_split = factor(rep(1:11, each = 3)),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_dend = FALSE,
  show_column_names = FALSE,
  row_title = NULL
)

p2 <- Heatmap(
  cor_matrix,
  name = "Correlation",
  col = col_fun,
  row_split = 6,
  column_split = 6,
  show_column_dend = FALSE,
  show_column_names = FALSE,
  column_title = NULL,
  row_title = NULL
)

pdf("./figures/FX_FXcryo_cryoFX_pseudobulk_cor.pdf", width = 9.6, height = 6)
p1
p2
dev.off()
