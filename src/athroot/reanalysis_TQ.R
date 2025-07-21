library(magrittr)
library(ggplot2)
library(ggrastr)
library(readxl)
library(Seurat)
library(dplyr)

# read data --------------
obj <- readRDS("./dataLib/athRoot/FX_TQ_v5.Rds")
obj <- subset(obj, tech == "tq")

# normal pipeline -------
obj_norm <- CreateSeuratObject(
  counts = GetAssayData(obj, layer = "counts"),
  meta.data = obj@meta.data
)

obj_norm %<>%
  NormalizeData() %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = round(0.1 * length(rownames(.)))
  ) %>%
  ScaleData(features = rownames(.), var.to.regress = "nCount_RNA") %>%
  RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)

vv <- VariableFeatures(obj_norm)

dims <- 1:30
obj_norm <- FindNeighbors(obj_norm, dims = dims)

obj_norm %<>%
  FindClusters(resolution = 0.6) %>%
  RunUMAP(dims = dims)

# remove proto genes (Bulk) which are in VariableFeatures -------
# read proto genes from Tom Denyer et al. 2019
proto_file <- "./dataLib/athRoot/2019_DC_athroot_Tom_Denyer.xlsx"
sheets <- excel_sheets(proto_file)
sheets

proto <- read_excel(proto_file, sheet = sheets[3])
colnames(proto) <- gsub(" ", "_", colnames(proto))

proto_genes_up <- proto$Gene_ID[proto$Log2_FC >= 2]
length(proto_genes_up)

proto_genes_down <- proto$Gene_ID[proto$Log2_FC <= -2]
length(proto_genes_down)

obj_rmproto <- CreateSeuratObject(
  counts = GetAssayData(obj, layer = "counts"),
  meta.data = obj@meta.data
)

obj_rmproto %<>%
  NormalizeData() %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = round(0.1 * length(rownames(.)))
  )

table(VariableFeatures(obj_rmproto) %in% c(proto_genes_up, proto_genes_down))
gg_bulk <- VariableFeatures(obj_rmproto)[VariableFeatures(obj_rmproto) %in% c(proto_genes_up, proto_genes_down)]
gg_bulk <- vv[vv %in% c(proto_genes_up, proto_genes_down)]


VariableFeatures(obj_rmproto) <- VariableFeatures(obj_rmproto)[
  !VariableFeatures(obj_rmproto) %in% c(proto_genes_up, proto_genes_down)
]

obj_rmproto %<>%
  ScaleData(features = rownames(.), var.to.regress = "nCount_RNA") %>%
  RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)


dims <- 1:30
obj_rmproto <- FindNeighbors(obj_rmproto, dims = dims)

obj_rmproto %<>%
  FindClusters(resolution = 0.6) %>%
  RunUMAP(dims = dims)

# remove proto genes (scRNA-seq) which are in VariableFeatures -------
# read proto genes from scRNA-seq
fp <- "./dataLib/athRoot/FX_TQ_protolasted_only_pos_TRUE.xlsx"
sheets <- excel_sheets(fp)
sheets <- sheets[1:12]

de_data <- lapply(sheets, function(x) {
  read_excel(fp, sheet = x) %>%
    subset(abs(avg_log2FC) >= 2) %>%
    pull(gene)
})

proto_genes_scrna <- unique(unlist(de_data))

obj_rmproto_scrna <- CreateSeuratObject(
  counts = GetAssayData(obj, layer = "counts"),
  meta.data = obj@meta.data
)

obj_rmproto_scrna %<>%
  NormalizeData() %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = round(0.12 * length(rownames(.)))
  )

table(VariableFeatures(obj_rmproto_scrna) %in% proto_genes_scrna)
gg_scrna <- VariableFeatures(obj_rmproto_scrna)[VariableFeatures(obj_rmproto_scrna) %in% proto_genes_scrna]
gg_scrna <- vv[vv %in% proto_genes_scrna]


VariableFeatures(obj_rmproto_scrna) <- VariableFeatures(obj_rmproto_scrna)[
  !VariableFeatures(obj_rmproto_scrna) %in% proto_genes_scrna
]

obj_rmproto_scrna %<>%
  ScaleData(features = rownames(.), var.to.regress = "nCount_RNA") %>%
  RunPCA(features = VariableFeatures(.), npcs = 100, verbose = FALSE)


dims <- 1:30
obj_rmproto_scrna <- FindNeighbors(obj_rmproto_scrna, dims = dims)

obj_rmproto_scrna %<>%
  FindClusters(resolution = 0.6) %>%
  RunUMAP(dims = dims)


# plot --------------
emb_norm <- Embeddings(obj_norm, reduction = "umap") %>%
  as.data.frame() %>%
  mutate(cell = rownames(.))
colnames(emb_norm) <- c("UMAP1_norm", "UMAP2_norm", "cell")

emb_rmproto <- Embeddings(obj_rmproto, reduction = "umap") %>%
  as.data.frame() %>%
  mutate(cell = rownames(.))
colnames(emb_rmproto) <- c("UMAP1_rmproto", "UMAP2_rmproto", "cell")

emb_rmproto_scrna <- Embeddings(obj_rmproto_scrna, reduction = "umap") %>%
  as.data.frame() %>%
  mutate(cell = rownames(.))
colnames(emb_rmproto_scrna) <- c("UMAP1_rmproto_scrna", "UMAP2_rmproto_scrna", "cell")

meta <- obj@meta.data
meta$cell <- rownames(meta)

pd <- left_join(emb_norm, emb_rmproto, by = "cell") %>%
  left_join(emb_rmproto_scrna, by = "cell") %>%
  left_join(meta, by = "cell")
head(pd)


p1 <- ggplot(pd, aes(x = UMAP1_norm, y = UMAP2_norm, color = celltype)) +
  geom_point_rast(size = 0.3) +
  ggtitle("Normal pipeline")

p2 <- ggplot(pd, aes(x = UMAP1_rmproto, y = UMAP2_rmproto, color = celltype)) +
  geom_point_rast(size = 0.3) +
  ggtitle("Remove proto genes (bulk) which are in VariableFeatures ")

p3 <- ggplot(pd, aes(x = UMAP1_rmproto_scrna, y = UMAP2_rmproto_scrna, color = celltype)) +
  geom_point_rast(size = 0.3) +
  ggtitle("Remove proto genes (scRNA-seq) which are in VariableFeatures")

colors <- pigRutils::select_colors("col_33")[14:28]
colors[13] <- "grey"

pdf("./figures/athroot/reanalysis_TQ_uamp.pdf", width = 20, height = 5.2)
(p1 + p2 + p3) & scale_color_manual(values = colors) & tidydr::theme_dr() +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  ) &
  guides(
    color = guide_legend(
      keywidth = 0.4,
      ncol = 1,
      override.aes = list(size = 2, stroke = 2)
    )
  )
dev.off()


# FindAllMarkers --------------
Idents(obj_norm) <- "seurat_clusters"
markers_norm <- FindAllMarkers(obj_norm, logfc.threshold = 2)

Idents(obj_rmproto) <- "seurat_clusters"
markers_rmproto <- FindAllMarkers(obj_rmproto, logfc.threshold = 2)

Idents(obj_rmproto_scrna) <- "seurat_clusters"
markers_rmproto_scrna <- FindAllMarkers(obj_rmproto_scrna, logfc.threshold = 2)


markers_norm <- markers_norm[markers_norm$p_val_adj < 0.05, ]
markers_rmproto <- markers_rmproto[markers_rmproto$p_val_adj < 0.05, ]
markers_rmproto_scrna <- markers_rmproto_scrna[markers_rmproto_scrna$p_val_adj < 0.05, ]



library(eulerr)

eu <- euler(
  list(
    "Normal pipeline" = unique(markers_norm$gene),
    "Remove proto (bulk)" = unique(markers_rmproto$gene),
    "Remove proto (scRNA-seq)" = unique(markers_rmproto_scrna$gene)
  )
)


p_venn <- plot(eu,
  quantities = list(type = c("percent")),
  edges = list(col = "white", lex = 1),
  fills = c("#87cdff", "#ffe96e", "#ff8b87"),
  legend = list(side = "right")
)

pdf("./figures/athroot/reanalysis_TQ_markers_venn.pdf", width = 8, height = 6)
print(p_veen)
dev.off()


eu2 <- euler(
  list(
    "Remove proto (bulk)" = gg_bulk,
    "Remove proto (scRNA)" = gg_scrna
  )
)

p_venn2 <- plot(eu2,
  quantities = list(type = c("counts")),
  edges = list(col = "white", lex = 1),
  fills = c("#87cdff", "#ffc2cc"),
  legend = list(side = "right")
)

pdf("./figures/athroot/reanalysis_TQ_remove_vars_venn.pdf", width = 8, height = 6)
print(p_venn2)
dev.off()


pse <- read.csv("./dataLib/athRoot/TQ_vs_FX_protolasted.csv", row.names = 1)
head(pse)



eu3 <- euler(
  list(
    "Remove proto (bulk)" = gg_bulk,
    "Protolasted" = rownames(pse)[rownames(pse) %in% VariableFeatures(obj_rmproto_scrna)]
  )
)

p_venn3 <- plot(eu3,
  quantities = list(type = c("counts")),
  edges = list(col = "white", lex = 1),
  fills = c("#87cdff", "#ffc2cc"),
  legend = list(side = "right")
)


