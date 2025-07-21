library(patchwork)
library(ggrastr)
library(ggplot2)
library(Seurat)
library(dplyr)
library(pigRutils)

# LC UMPA & TSNE -----------------
# load data
obj <- readRDS("./dataLib/athLeaf/LC.rds")
Reductions(obj)
colnames(obj@meta.data)

# prepare plot data
pdumap <- SeuratObject::Embeddings(obj, reduction = "umap") %>%
  as_tibble(rownames = "cell") %>%
  mutate(c = obj$seurat_clusters)

pdtsne <- SeuratObject::Embeddings(obj, reduction = "tsne") %>%
  as_tibble(rownames = "cell")

pd <- inner_join(pdumap, pdtsne, by = "cell")
head(pd)

pumap <- ggplot() +
  ggrastr::geom_point_rast(data = pd, aes(x = umap_1, y = umap_2, color = c), size = 0.1)

ptsne <- ggplot() +
  ggrastr::geom_point_rast(data = pd, aes(x = tSNE_1, y = tSNE_2, color = c), size = 0.1)

p <- (pumap + ptsne) &
  scale_color_manual(values = pigRutils::select_colors("col_33")) &
  tidydr::theme_dr() &
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) &
  guides(color = guide_legend(keywidth = 0.4, ncol = 1, override.aes = list(size = 2, stroke = 2)))

ggsave("./figures/leaf/LC_umap_tsne.pdf", p, width = 12, height = 6)

# LPZ UMPA & TSNE -----------------
# load data
obj <- readRDS("./dataLib/athLeaf/LPZ.rds")
Reductions(obj)
colnames(obj@meta.data)

# prepare plot data
pdumap <- SeuratObject::Embeddings(obj, reduction = "umap") %>%
  as_tibble(rownames = "cell") %>%
  mutate(c = obj$seurat_clusters)

pdtsne <- SeuratObject::Embeddings(obj, reduction = "tsne") %>%
  as_tibble(rownames = "cell")

pd <- inner_join(pdumap, pdtsne, by = "cell")
head(pd)

pumap <- ggplot(data = pd, aes(x = umap_1, y = umap_2, color = c)) +
  ggrastr::geom_point_rast(size = 0.1)

ptsne <- ggplot(data = pd, aes(x = tSNE_1, y = tSNE_2, color = c)) +
  ggrastr::geom_point_rast(size = 0.1)

p <- (pumap + ptsne) &
  scale_color_manual(values = pigRutils::select_colors("col_33")) &
  tidydr::theme_dr() &
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) &
  guides(color = guide_legend(keywidth = 0.4, ncol = 1, override.aes = list(size = 2, stroke = 2)))

ggsave("./figures/leaf/LPZ_umap_tsne.pdf", p, width = 13, height = 6)
# leaf: LC & LPZ & Public data -----------------
# load data
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public.Rds")
Reductions(obj)
colnames(obj@meta.data)


# Anno UMAP merge/split -----------------
# merge
table(obj$celltype)
col <- pigRutils::select_colors("col_33")
col <- col[c(1, 7, 5, 3, 6, 2, 14, 9, 10, 8, 12)]
embSCdim(
  obj, "harmony_umap",
  group_by = "celltype",
  colors = col,
  save = TRUE,
  width = 7.5, height = 5.6,
  save_path = "./figures/leaf/LCLPZ_pub_umap_with_anno.pdf"
)
# split
embSCdim(
  obj, "harmony_umap",
  group_by = "celltype",
  split_by = "orig.ident",
  colors = col,
  save = TRUE,
  width = 17.2, height = 5.5,
  save_path = "./figures/leaf/LCLPZ_pub_umap_with_anno_split.pdf"
)

# Cluster UMAP merge/split ------------------
# merge
embSCdim(
  obj, "harmony_umap",
  group_by = "celltype",
  colors = col,
  save = TRUE,
  width = 7, height = 6,
  save_path = "./figures/leaf/LCLPZ_pub_umap_with_Cluster.pdf"
)

# split
embSCdim(
  obj, "harmony_umap",
  group_by = "celltype",
  split_by = "orig.ident",
  colors = col,
  save = TRUE,
  width = 15.6, height = 5,
  save_path = "./figures/leaf/LCLPZ_pub_umap_with_Cluster_split.pdf"
)
# Sample UMAP merge/split ------------------
# prepare plot data
pd <- SeuratObject::Embeddings(obj, reduction = "harmony_umap") %>%
  as_tibble(rownames = "cell") %>%
  mutate(c = obj$seurat_clusters, cc = obj$celltype, sample = obj$orig.ident)

colnames(pd)[2:3] <- c("UMAP_1", "UMAP_2")

# split
pd$sample <- factor(pd$sample, levels = c("LC", "LPZ", "JC"))

p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
  ggrastr::geom_point_rast(size = 0.2) +
  scale_color_manual(values = c("#1dc1ab", "#56acfc", "#fc877f")) +
  facet_wrap(~sample) +
  tidydr::theme_dr() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  guides(
    color = guide_legend(
      keywidth = 0.4,
      ncol = 1,
      override.aes = list(size = 2, stroke = 2)
    )
  )

ggsave("./figures/leaf/LCLPZ_pub_umap_with_sample_split.pdf", width = 15.2, height = 5)

# merge
lc <- pd[pd$sample == "LC", ]
lpz <- pd[pd$sample == "LPZ", ]
jc <- pd[pd$sample == "JC", ]

pd_merge <- rbind(lpz, lc, jc)

p <- ggplot(data = pd_merge, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
  ggrastr::geom_point_rast(size = 0.1) +
  scale_color_manual(values = c("#1dc1ab", "#56acfc", "#fc877f")) +
  tidydr::theme_dr() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  guides(
    color = guide_legend(
      keywidth = 0.4,
      ncol = 1,
      override.aes = list(size = 2, stroke = 2)
    )
  )

ggsave("./figures/leaf/LCLPZ_pub_umap_with_sample.pdf", width = 8.1, height = 7)

# rish umap -----------------
mm <- c(
  "AT3G17609", "AT1G56060", "AT2G02990",
  "AT3G04300", "AT2G35290", "AT5G54300",
  "AT1G22810", "AT4G34410", "AT1G65390"
)

p <- FeaturePlot(obj, features = mm, reduction = "harmony_umap")

meta <- obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(cell, orig.ident)

p_list <- lapply(seq.int(p), function(x) {
  xx <- p[[x]]
  pd <- xx$data %>%
    as_tibble(rownames = "cell") %>%
    left_join(meta)
  t <- xx$labels$title
  print(t)
  ggplot(pd, aes(x = harmonyumap_1, y = harmonyumap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    facet_wrap(~orig.ident) +
    ggtitle(t)
})

pp <- wrap_plots(p_list, ncol = 2) &
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

pdf("./figures/leaf/leaf_rish.pdf", width = 17.2, height = 16)
pp
dev.off()

# go AUCell -----------------
library(AUCell)
library(scales)

go <- read.csv("./dataLib/athLeaf/GO_response_to_wounding.tsv", sep = "\t", header = FALSE)
go <- go[[1]]

obj <- JoinLayers(obj)
Layers(obj)
cts <- LayerData(obj, layer = "counts")

GeneSets <- list(
  wound = unique(go)
)

CellsAUC <- AUCell_run(cts, GeneSets)

pd <- Embeddings(obj, reduction = "harmony_umap") %>% as.data.frame()
pd$auc <- pigRutils::mf_minmax(getAUC(CellsAUC)["wound", ])
pd$group <- obj$orig.ident
head(pd)

p <- ggplot(pd, aes(x = harmonyumap_1, y = harmonyumap_2, color = auc)) +
  ggrastr::geom_point_rast(size = 0.3) +
  facet_wrap(~group) +
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) +
  theme_linedraw() +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  ggtitle("AUCell Score: Response to wounding")

ggsave("./figures/leaf/goID_wound_AUCell_split.pdf", p, height = 6.2, width = 18)

# wound response gene umap plot -----------
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public.Rds")

genes <- c(
  "AT5G64750", "AT4G34410", "AT3G59220", "AT3G28210", "AT3G25780",
  "AT3G23250", "AT3G02840", "AT2G38530", "AT2G29500", "AT2G02990",
  "AT1G53540", "AT1G07400", "AT1G02920", "AT1G02400"
)

p <- FeaturePlot(obj, features = genes, reduction = "harmony_umap")

meta <- obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(cell, orig.ident)

sapply(seq.int(p), function(x) {
  xx <- p[[x]]
  pd <- xx$data %>%
    as_tibble(rownames = "cell") %>%
    left_join(meta)
  t <- xx$labels$title
  print(t)
  p <- ggplot(pd, aes(x = harmonyumap_1, y = harmonyumap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    facet_wrap(~orig.ident) +
    scale_color_gradientn(colours = c("#c8c6c3", "#f6ee6c", "#f9a432", "#eb212f", "#88181d")) +
    ggtitle(t) +
    theme_linedraw() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )

  ggsave(paste0("./Plots/athleaf/wound/", t, ".pdf"), p, height = 5, width = 13.3)
})
