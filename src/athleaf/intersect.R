library(clusterProfiler)
library(org.At.tair.db)
library(ComplexHeatmap)
library(ggplot2)
library(Seurat)
library(dplyr)

# Read reference data
alias <- read.csv("~/ref/Athaliana/Araport11/gene_aliases_20220331_tidy.csv")

# Define cell types
cell_types <- c(
  "epi" = "Epidermis", "meso" = "Mesophyll", "bs" = "Bundle-Sheet",
  "ph" = "Phloem", "hyd" = "Hydathode"
)

# Function to read DE gene data
read_de_data <- function(cell_type_name) {
  file_path <- paste0("./dataLib/athLeaf/de/LC_LPZ_", cell_type_name, "_de_sig.csv")
  data <- read.csv(file_path)
  data
}

# Read DE data for all cell types
de_data <- lapply(cell_types, read_de_data)

# Print dimensions of each dataset
lapply(de_data, dim)

# Function to separate up and down regulated genes
get_regulated_genes <- function(de_data) {
  up <- de_data$geneId[de_data$avg_log2FC > 0]
  down <- de_data$geneId[de_data$avg_log2FC < 0]
  list(up = up, down = down)
}

# Get up and down regulated genes for each cell type
regulated_genes <- lapply(de_data, get_regulated_genes)

# Create lists for upset plots
up_list <- lapply(regulated_genes, function(x) x$up)
down_list <- lapply(regulated_genes, function(x) x$down)

# Name the lists with cell type short names
names(up_list) <- names(cell_types)
names(down_list) <- names(cell_types)

# Print counts
sapply(up_list, length)
sapply(down_list, length)

all_intersect <- Reduce(intersect, down_list)
write.csv(all_intersect, "./dataLib/athLeaf/intersect_all_42.csv", row.names = FALSE)

# Create combination matrices for upset plots
m_up <- make_comb_mat(up_list)
m_down <- make_comb_mat(down_list)

# Set colors based on combination degree
comb_colors <- c("#fe1f1f", "#ffa600", "#00eaff", "#1d1df9", "black")

# Create upset plot for up-regulated genes
p1 <- UpSet(
  m_up,
  set_order = c("meso", "epi", "bs", "ph", "hyd"),
  column_title = "Up-regulated (Distinct Mode)",
  comb_col = comb_colors[comb_degree(m_up)],
  top_annotation = upset_top_annotation(m_up, add_numbers = TRUE, annotation_name_rot = 90),
  left_annotation = upset_left_annotation(m_up, add_numbers = TRUE, show_annotation_name = FALSE)
)

# Create upset plot for down-regulated genes
p2 <- UpSet(
  m_down,
  set_order = c("meso", "epi", "bs", "ph", "hyd"),
  column_title = "Down-regulated (Distinct Mode)",
  comb_col = comb_colors[comb_degree(m_down)],
  top_annotation = upset_top_annotation(m_down, add_numbers = TRUE, show_annotation_name = FALSE),
  right_annotation = upset_right_annotation(m_down, add_numbers = TRUE, show_annotation_name = FALSE)
)

# Save plots to PDF
pdf("./figures/leaf/upset_up.pdf", width = 5, height = 3)
p1
dev.off()

pdf("./figures/leaf/upset_down.pdf", width = 7, height = 3)
p2
dev.off()

# Function to find unique genes for each cell type
get_unique_genes <- function(target_list_name, all_lists) {
  Reduce(
    setdiff,
    c(
      list(all_lists[[target_list_name]]),
      all_lists[names(all_lists) != target_list_name]
    )
  )
}

# Find unique up-regulated genes for each cell type
unique_up <- lapply(names(up_list), function(cell_type) {
  get_unique_genes(cell_type, up_list)
})
names(unique_up) <- names(cell_types)

# Find unique down-regulated genes for each cell type
unique_down <- lapply(names(down_list), function(cell_type) {
  get_unique_genes(cell_type, down_list)
})
names(unique_down) <- names(cell_types)

# Print counts of unique genes
sapply(unique_up, length)
sapply(unique_down, length)


# load Seurat obj
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public_anno.Rds")
obj <- JoinLayers(obj)
obj$celltype <- gsub("Mesophyll\\*", "Mesophyll", obj$celltype)
obj <- subset(obj, subset = orig.ident %in% c("LC", "LPZ"))
sel_ct <- c("Phloem", "Bundle-Sheet", "Epidermis", "Hydathode", "Mesophyll")
obj <- subset(obj, subset = celltype %in% sel_ct)
meta <- obj@meta.data

# Get scaled data
Layers(obj)
cts <- GetAssayData(obj, layer = "scale.data")

# Get up and down regulated genes
gup <- unique_up %>% unlist()
gup <- data.frame(gene = gup, ct = gsub("[0-9]+", "", names(gup)))
dim(gup)

gdown <- unique_down %>% unlist()
gdown <- data.frame(gene = gdown, ct = gsub("[0-9]+", "", names(gdown)))
dim(gdown)

# Get up and down Matrix
meta <- meta %>%
  as_tibble(rownames = "cell") %>%
  group_by(orig.ident, celltype) %>%
  arrange(orig.ident, celltype)

cts_up <- cts[gup$gene, meta$cell]
cts_down <- cts[gdown$gene, meta$cell]

# up
pd <- cts_up
pd <- MinMax(pd, min = -2.5, max = 2.5)

col_fun <- circlize::colorRamp2(
  c(min(pd), 0.5, max(pd)) %>% round(., 1),
  c("#bbdddd", "#dbdada", "#ff0000")
)

# Define custom colors for cell types and groups
celltype_colors <- c(
  "Phloem" = "#e4c5c6",
  "Bundle-Sheet" = "#90b9c7",
  "Epidermis" = "#f9c59c",
  "Hydathode" = "#f17a51",
  "Mesophyll" = "#60bba2"
)

group_colors <- c(
  "LC" = "#c575ff",
  "LPZ" = "#7784ff"
)

ha <- HeatmapAnnotation(
  CellType = meta$celltype,
  Groups = meta$orig.ident,
  col = list(
    CellType = celltype_colors,
    Groups = group_colors
  )
)

pdf("./figures/leaf/heatmap_up_distinct.pdf", width = 15, height = 9)
Heatmap(
  pd,
  name = "Scaled Expression",
  row_split = gup$ct,
  column_split = meta$celltype,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = NULL,
  top_annotation = ha,
  col = col_fun,
  use_raster = TRUE
)
dev.off()

# down
pd <- cts_down
pd <- MinMax(pd, min = -2.5, max = 2.5)

col_fun <- circlize::colorRamp2(
  c(min(pd), 0.5, max(pd)) %>% round(., 1),
  c("#bbdddd", "#dbdada", "#ff0000")
)

# Define custom colors for cell types and groups
celltype_colors <- c(
  "Phloem" = "#e4c5c6",
  "Bundle-Sheet" = "#90b9c7",
  "Epidermis" = "#f9c59c",
  "Hydathode" = "#f17a51",
  "Mesophyll" = "#60bba2"
)

group_colors <- c(
  "LC" = "#c575ff",
  "LPZ" = "#7784ff"
)

ha <- HeatmapAnnotation(
  CellType = meta$celltype,
  Groups = meta$orig.ident,
  col = list(
    CellType = celltype_colors,
    Groups = group_colors
  )
)

pdf("./figures/leaf/heatmap_down_distinct.pdf", width = 15, height = 9)
Heatmap(
  pd,
  name = "Scaled Expression",
  row_split = gdown$ct,
  column_split = meta$celltype,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = NULL,
  top_annotation = ha,
  col = col_fun,
  use_raster = TRUE
)
dev.off()

# heatmap with single celltype
pup_list <- lapply(unique(gup$ct), function(ct) {
  print(ct)

  cellt <- cell_types[ct]
  sub_meta <- meta[meta$celltype == cellt, ]
  sub_cts <- cts_up[gup$gene[gup$ct == ct], sub_meta$cell]
  sub_cts <- MinMax(sub_cts, min = -2.5, max = 2.5)

  ctc <- celltype_colors[cellt]

  col_fun <- circlize::colorRamp2(
    c(min(sub_cts), 0.5, max(sub_cts)) %>% round(., 1),
    c("#c6e1e1", "#fff2f2", "#ff0000")
  )
  ha <- HeatmapAnnotation(
    CellType = sub_meta$celltype,
    Groups = factor(sub_meta$orig.ident),
    col = list(
      CellType = ctc,
      Groups = group_colors
    )
  )
  p <- Heatmap(
    sub_cts,
    name = "Scaled Expression",
    column_split = factor(sub_meta$orig.ident),
    show_row_names = FALSE,
    show_column_dend = FALSE,
    show_column_names = FALSE,
    column_title = NULL,
    top_annotation = ha,
    col = col_fun,
    use_raster = TRUE
  )
  p
})

names(pup_list) <- unique(gup$ct)
table(meta$celltype)
table(gup$ct)

pdf("./figures/leaf/heatmap_up_distinct_Hydathode.pdf", width = 5.5, height = 1.8)
print(pup_list[["hyd"]])
dev.off()

pdf("./figures/leaf/heatmap_up_distinct_Bundle-Sheet.pdf", width = 6.3, height = 2.6)
print(pup_list[["bs"]])
dev.off()

pdf("./figures/leaf/heatmap_up_distinct_Phloem.pdf", width = 6, height = 2.6)
print(pup_list[["ph"]])
dev.off()

pdf("./figures/leaf/heatmap_up_distinct_Epidermis.pdf", width = 9, height = 5)
print(pup_list[["epi"]])
dev.off()

pdf("./figures/leaf/heatmap_up_distinct_Mesophyll.pdf", width = 10, height = 2.8)
print(pup_list[["meso"]])
dev.off()

pdown_list <- lapply(unique(gdown$ct), function(ct) {
  print(ct)

  cellt <- cell_types[ct]
  sub_meta <- meta[meta$celltype == cellt, ]
  sub_cts <- cts_down[gdown$gene[gdown$ct == ct], sub_meta$cell]
  sub_cts <- MinMax(sub_cts, min = -2.5, max = 2.5)

  ctc <- celltype_colors[cellt]

  col_fun <- circlize::colorRamp2(
    c(min(sub_cts), 0.5, max(sub_cts)) %>% round(., 1),
    c("#c6e1e1", "#fff2f2", "#ff0000")
  )
  ha <- HeatmapAnnotation(
    CellType = sub_meta$celltype,
    Groups = factor(sub_meta$orig.ident),
    col = list(
      CellType = ctc,
      Groups = group_colors
    )
  )
  p <- Heatmap(
    sub_cts,
    name = "Scaled Expression",
    column_split = factor(sub_meta$orig.ident),
    show_row_names = FALSE,
    show_column_dend = FALSE,
    show_column_names = FALSE,
    column_title = NULL,
    top_annotation = ha,
    col = col_fun,
    use_raster = TRUE
  )
  p
})

names(pdown_list) <- unique(gdown$ct)
table(meta$celltype)
table(gdown$ct)

pdf("./figures/leaf/heatmap_down_distinct_Hydathode.pdf", width = 5.5, height = 1.8)
print(pdown_list[["hyd"]])
dev.off()

pdf("./figures/leaf/heatmap_down_distinct_Bundle-Sheet.pdf", width = 6.3, height = 3.6)
print(pdown_list[["bs"]])
dev.off()

pdf("./figures/leaf/heatmap_down_distinct_Phloem.pdf", width = 6, height = 2.8)
print(pdown_list[["ph"]])
dev.off()

pdf("./figures/leaf/heatmap_down_distinct_Epidermis.pdf", width = 9, height = 4)
print(pdown_list[["epi"]])
dev.off()

pdf("./figures/leaf/heatmap_down_distinct_Mesophyll.pdf", width = 10, height = 6)
print(pdown_list[["meso"]])
dev.off()



for (ct in names(unique_up)) {
  print(ct)
  ff <- unique_up[[ct]]
  small_obj <- subset(obj, subset = celltype == cell_types[ct])

  if (!dir.exists(paste0("./Plots/athleaf/up_down/up/", ct))) {
    dir.create(paste0("./Plots/athleaf/up_down/up/", ct))
  }

  for (g in ff) {
    print(g)
    p <- FeaturePlot(
      small_obj,
      features = g,
      pt.size = 0.1,
      order = TRUE,
      split.by = "orig.ident"
    )
    ggsave(
      paste0("./Plots/athleaf/up_down/up/", ct, "/", g, ".png"),
      p,
      width = 6, height = 3
    )
  }
}


for (ct in names(unique_down)) {
  print(ct)
  ff <- unique_down[[ct]]
  small_obj <- subset(obj, subset = celltype == cell_types[ct])

  if (!dir.exists(paste0("./Plots/athleaf/up_down/down/", ct))) {
    dir.create(paste0("./Plots/athleaf/up_down/down/", ct))
  }

  for (g in ff) {
    print(g)
    p <- FeaturePlot(
      small_obj,
      pt.size = 0.1,
      features = g,
      order = TRUE,
      split.by = "orig.ident"
    )
    ggsave(
      paste0("./Plots/athleaf/up_down/down/", ct, "/", g, ".png"),
      p,
      width = 6, height = 3
    )
  }
}


openxlsx::write.xlsx(unique_up, "./dataLib/athLeaf/LC_LPZ_wound_distinct_up.xlsx")
openxlsx::write.xlsx(unique_down, "./dataLib/athLeaf/LC_LPZ_wound_distinct_down.xlsx")
