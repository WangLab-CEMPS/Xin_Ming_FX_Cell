library(dplyr)
library(tidyr)
library(readxl)
library(eulerr)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(patchwork)
library(ComplexHeatmap)


obj <- readRDS("./dataLib/athRoot/FX_TQ_v5.Rds")

fp <- "./dataLib/athRoot/FX_TQ_protolasted_only_pos_TRUE.xlsx"
sheets <- excel_sheets(fp)
sheets <- sheets[1:12]

de_data <- lapply(sheets, function(x) {
  read_excel(fp, sheet = x)
})
names(de_data) <- sheets
sapply(de_data, nrow)
head(de_data[[1]])

tq_fx_up <- lapply(de_data, function(x) {
  x$gene[x$avg_log2FC >= 2]
})
sapply(tq_fx_up, length)
length(unique(unlist(tq_fx_up)))

tq_fx_down <- lapply(de_data, function(x) {
  x$gene[x$avg_log2FC <= -2]
})
sapply(tq_fx_down, length)
length(unique(unlist(tq_fx_down)))

# proto genes from Tom Denyer et al. (2019) ----
proto_file <- "./dataLib/athRoot/2019_DC_athroot_Tom_Denyer.xlsx"
sheets <- excel_sheets(proto_file)
sheets

proto <- read_excel(proto_file, sheet = sheets[3])
colnames(proto) <- gsub(" ", "_", colnames(proto))
rm(sheets)

proto <- proto[proto$Gene_ID %in% rownames(obj), ]

proto_genes_up <- proto$Gene_ID[proto$Log2_FC >= 2]
length(proto_genes_up)

proto_genes_down <- proto$Gene_ID[proto$Log2_FC <= -2]
length(proto_genes_down)

# marker genes from FindAllMarkers ----
markers <- read.csv("./dataLib/athRoot/FX_TQ_FindAllMarkers_celltype.csv", row.names = 1)
rownames(markers) <- NULL
head(markers)

obj <- AddModuleScore(obj, features = list(tq_fx_up[["Root-Cap"]]), name = "tmp")
FeaturePlot(obj, features = "tmp1", split.by = "tech")

unique(markers$cluster)
marker <- markers[!markers$cluster %in% c("Epidermis", "Unknown"), ]
marker <- marker[marker$p_val_adj < 0.05 & marker$avg_log2FC > 2, ]
dim(marker)

pp1 <- plot(
  euler(
    list(
      "tq_fx_up" = unique(unlist(tq_fx_up)),
      "bulk_proto_up" = unique(proto_genes_up)
    )
  ),
  quantities = list(type = c("percent", "counts")),
  edges = list(col = "grey90", lex = 3),
  fills = c("#ff8b87", "#ffe46e")
)

pp2 <- plot(
  euler(
    list(
      "tq_fx_down" = unique(unlist(tq_fx_down)),
      "bulk_proto_down" = unique(proto_genes_down)
    )
  ),
  quantities = list(type = c("percent", "counts")),
  edges = list(col = "grey90", lex = 3),
  fills = c("#ff8b87", "#ffe46e")
)

pdf("./figures/athroot/FX_TQ_proto_veen.pdf", width = 8, height = 4)
print(pp1)
print(pp2)
dev.off()


# --------
marker_genes <- marker[marker$cluster %in% c("Endodermis"), ]
marker_genes <- marker_genes[marker_genes$gene %in% unique(tq_fx_up[["Endodermis"]]), ]
dim(marker_genes)

marker_genes %>% arrange(desc(pct.1)) %>% top_n(10)

mm <- c("AT1G66700", "AT5G07310", "AT5G64230", "AT5G62140", "AT5G59540", "AT5G38005")
mm %in% unique(unlist(tq_fx_up))
mm %in% unique(proto_genes_up)

FeaturePlot(obj, features = "AT5G38005", split.by = "tech")

ppp <- FeaturePlot(obj, features = mm)

meta <- obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(cell, tech)

meta$tech <- ifelse(meta$tech == "tq", "TQ", "FX-Cell")

p_list <- lapply(seq.int(ppp), function(x) {
  xx <- ppp[[x]]
  pd <- xx$data %>%
    as_tibble(rownames = "cell") %>%
    left_join(meta)
  t <- xx$labels$title
  print(t)
  ggplot(pd, aes(x = UMAP_1, y = UMAP_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    facet_wrap(~tech) +
    ggtitle(t)
})

pp <- wrap_plots(p_list, ncol = 1) &
  scale_color_gradientn(colours = c("#c8c6c3", "#f6ee6c", "#f9a432", "#eb212f", "#88181d")) &
  theme_linedraw() &
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

pdf("./figures/athroot/FX_TQ_proto_feature.pdf", width = 8.6, height = 27)
pp
dev.off()

# upset plot ----
select_celltype <- c("QC", "Cortex", "Endodermis", "Root-Cap", "Stele")

m_up <- make_comb_mat(tq_fx_up[select_celltype])
m_down <- make_comb_mat(tq_fx_down[select_celltype])



# Set colors based on combination degree
comb_colors <- c("#fe1f1f", "#ffa600", "#00eaff", "#1d1df9", "black")

# Create upset plot for up-regulated genes
p1 <- UpSet(
  m_up,
  set_order = c("Root-Cap", "Stele", "QC", "Endodermis", "Cortex"),
  column_title = "Up-regulated (Distinct Mode)",
  comb_col = comb_colors[comb_degree(m_up)],
  top_annotation = upset_top_annotation(m_up, add_numbers = TRUE, annotation_name_rot = 90),
  left_annotation = upset_left_annotation(m_up, add_numbers = TRUE, show_annotation_name = FALSE)
)

# Create upset plot for down-regulated genes
p2 <- UpSet(
  m_down,
  set_order = c("Root-Cap", "Stele", "QC", "Endodermis", "Cortex"),
  column_title = "Down-regulated (Distinct Mode)",
  comb_col = comb_colors[comb_degree(m_down)],
  top_annotation = upset_top_annotation(m_down, add_numbers = TRUE, show_annotation_name = FALSE),
  right_annotation = upset_right_annotation(m_down, add_numbers = TRUE, show_annotation_name = FALSE)
)

# Save plots to PDF
pdf("./figures/athroot/upset_up.pdf", width = 8, height = 3)
p1
dev.off()

pdf("./figures/athroot/upset_down.pdf", width = 7, height = 3)
p2
dev.off()
