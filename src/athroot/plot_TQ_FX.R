library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(patchwork)

# read data --------------
obj <- readRDS("./dataLib/athRoot/FX_TQ_v5.Rds")
obj$seurat_clusters <- as.character(obj$seurat_clusters)
obj$seurat_clusters[obj$celltype == "Epidermis"] <- "22"
obj$seurat_clusters <- factor(obj$seurat_clusters, levels = 0:22)

pd <- SeuratObject::Embeddings(obj, reduction = "umap") %>%
  as_tibble(rownames = "cell") %>%
  mutate(c = obj$seurat_clusters, ct = obj$celltype, sample = obj$orig.ident, tech = obj$tech)

colnames(pd)[2:3] <- c("UMAP_1", "UMAP_2")

# 注释 UMAP -----------------
colors <- pigRutils::select_colors("col_33")[14:28]
colors[13] <- "grey"
p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = ct)) +
  ggrastr::geom_point_rast(size = 0.2) +
  scale_color_manual(values = colors) +
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

ggsave("./figures/athroot/umap_anno.pdf", width = 7.2, height = 5)

p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = ct)) +
  ggrastr::geom_point_rast(size = 0.2) +
  scale_color_manual(values = colors) +
  facet_wrap(~tech) +
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

ggsave("./figures/athroot/umap_anno_split.pdf", width = 10.8, height = 5)

# merge/split -----------------

p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
  ggrastr::geom_point_rast(size = 0.2) +
  scale_color_manual(values = pigRutils::select_colors("col_33")[8:32]) +
  facet_wrap(~tech) +
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
      ncol = 2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

ggsave("./figures/athroot/umap_cluster_split.pdf", width = 10.1, height = 5)

p <- ggplot(data = pd, aes(x = UMAP_1, y = UMAP_2, color = c)) +
  ggrastr::geom_point_rast(size = 0.3) +
  scale_color_manual(values = pigRutils::select_colors("col_33")[8:32]) +
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
      ncol = 2,
      override.aes = list(size = 2, stroke = 2)
    )
  )

ggsave("./figures/athroot/umap_cluster_split_sample.pdf", width = 15.1, height = 5)

# markers umap-----------------
mm <- c(
  "AT1G33280", "AT1G79580",
  "AT2G26760", "AT3G25980",
  "AT3G20840", "AT1G51190",
  "AT5G49270", "AT1G16440",
  "AT1G29025", "AT1G62510",
  "AT4G11290", "AT5G57620",
  "AT4G14020", "AT4G19840",
  "AT1G68810", "AT1G20850",
  "AT1G77690", "AT4G14650",
  "AT2G28790", "AT1G73590",
  "AT1G79840", "AT5G40330",
  "AT1G68740", "AT1G32450",
  "AT3G23430", "AT3G45700"
)

p <- FeaturePlot(obj, mm, order = TRUE)

p_list <- lapply(p, function(x) {
  dd <- x$data %>%
    as_tibble(rownames = "cell") %>%
    left_join(pd)
  t <- x$labels$title

  ggplot(dd, aes(x = UMAP_1, y = UMAP_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    facet_wrap(~sample) +
    ggtitle(t)
})


pp <- wrap_plots(p_list, ncol = 3) &
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

pdf("./figures/athroot/markers_umap_split.pdf", width = 25, height = 30)
pp
dev.off()
