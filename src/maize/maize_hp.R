library(dplyr)
library(readxl)
library(Seurat)
library(circlize)
library(ComplexHeatmap)

# load seurat obj and labeled markers
obj <- readRDS("./dataLib/maize/maize.rds")

mm <- obj@meta.data %>%
  select(seurat_clusters, celltype) %>%
  unique()
rownames(mm) <- paste0("C", mm$seurat_clusters)

cell_order <- c(1, 5, 4, 16, 3, 6, 18, 11, 8, 13, 14, 2, 15, 12, 19, 17, 7, 10, 0, 9)
cell_order <- paste0("C", cell_order)
mm <- mm[cell_order, ]

# load all marker list
marker_file <- "./dataLib/maize/maize_markers.xlsx"

markers <- purrr::map(
  cell_order, ~ read_excel(marker_file, sheet = .)
)
names(markers) <- cell_order

tt <- lapply(cell_order, function(x) {
  markers[[x]] %>%
    mutate(cluster = x) %>%
    filter(avg_log2FC >= 1.32) %>%
    top_n(n = 200, wt = avg_log2FC) %>%
    dplyr::pull(gene)
})

names(tt) <- cell_order

features <- c(
  "Zm00001eb372230", "Zm00001eb422410", "Zm00001eb249210",
  "Zm00001eb370670", "Zm00001eb035350", "Zm00001eb188180",
  "Zm00001eb266840", "Zm00001eb315650", "Zm00001eb002390", "Zm00001eb342070",
  "Zm00001eb397500", "Zm00001eb302250",
  "Zm00001eb175860", "Zm00001eb338510",
  "Zm00001eb267050",
  "Zm00001eb249760", "Zm00001eb429560",
  "Zm00001eb069630",
  "Zm00001eb187430",
  "Zm00001eb200450",
  "Zm00001eb091140", "Zm00001eb091170", "Zm00001eb140000"
)


alias <- read.csv("./dataLib/maize/genes_all.csv") %>%
  subset(v5GeneModelID %in% features) %>%
  select(v5GeneModelID, GeneSymbol) %>%
  unique()
rownames(alias) <- alias$v5GeneModelID

# adjust gene loc
all_genes <- c(unlist(tt), features) %>% unique()

# get AverageExpression use log normalized data
obj[["cluster"]] <- paste0("C", obj$seurat_clusters)
avg <- AggregateExpression(obj, group.by = "cluster", slot = "data")
avg <- as.matrix(avg$RNA)
avg <- avg[all_genes, ]

# outlier
q05 <- quantile(avg, .05)
q95 <- quantile(avg, .95)
avg[avg < q05] <- q05
avg[avg > q95] <- q95

# scale for heatmap
pd <- t(scale(t(avg[, cell_order])))
pd[pd >= 3.5] <- 3.5
pd[pd <= -3.5] <- -3.5

col_fun <- colorRamp2(c(-3.5, 0, 3.5), c("#16bbd8", "#ecf5f6", "#ed4300"))

obj$celltype <- gsub("CortexII\\(m\\)", "Cortex Initials", obj$celltype)
# readin colors info
cell_colors <- data.frame(
  colors = c(
    "#b7889d", "#ecd0d0", # "#c6dbef",
    "#a9e097", "#04a5fc",
    "#66c2a5", "#74c476", "#8dcaec", "#f01778", "#cab9f5", "#df8633", "#2c6917"
  ),
  celltype = levels(factor(obj$celltype))
)
rownames(cell_colors) <- cell_colors$celltype

mm$celltype <- gsub("CortexII\\(m\\)", "Cortex Initials", mm$celltype)
Group <- factor(mm$celltype)

top_annotation <- HeatmapAnnotation(
  cluster = anno_block(
    gp = gpar(fill = cell_colors$colors),
    labels_gp = gpar(col = "black", fontsize = 9),
    labels = cell_colors$celltype,
    width = unit(1, "cm")
  ),
  border = TRUE
)

alias$GeneSymbol %in% c("Zm00001d002191")

hp <- Heatmap(
  pd,
  col = col_fun,
  column_split = Group,
  column_order = cell_order,
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_title = "Scaled Expression of Cell-Type-Specific Markers",
  top_annotation = top_annotation,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  use_raster = TRUE,
  raster_resize_mat = TRUE,
  heatmap_legend_param = list(
    at = c(-3.5, -1.5, 0, 1.5, 3.5),
    direction = "horizontal",
    title_position = "topcenter",
    legend_width = unit(3.4, "cm"),
    title = "Normalized\nExpression"
  )
) + rowAnnotation(
  link = anno_mark(
    at = which(rownames(avg) %in% alias$v5GeneModelID),
    labels = alias[rownames(avg)[rownames(avg) %in% alias$v5GeneModelID], "GeneSymbol"],
    which = "rows"
  )
)

pdf("./figures/maize/maize_hp.pdf", width = 12, height = 6)
hp
dev.off()
