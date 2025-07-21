library(Seurat)
library(ggplot2)
library(ggrastr)
library(dplyr)
library(patchwork)

obj <- readRDS("./dataLib/tillering_node_rhizome/rhizome_S0.rds")
DimPlot(obj, group.by = "celltype")

mm <- c(
  "Bochr02G072640", "Bochr03G119130", "Bochr01G001790",
  "Bochr03G117850", "Bochr12G355820", "Bochr01G043670",
  "Bochr10G305000", "Bochr06G205260", "Bochr01G014440"
)

# mm <- c(
#   "Os01g0642200", "Os02g0232900", "Os04g0496300",
#   "Os01g0507700", "Os05g0550300", "Os05g0160300",
#   "Os09g0458900", "Os03g0135100", "Os03g0807500",
#   "Os06g0218600", "Os09g0451400", "Os12g0476200"
# )

p <- FeaturePlot(obj, mm, order = TRUE)


pd <- SeuratObject::Embeddings(obj, reduction = "umap") %>%
  as_tibble(rownames = "cell") %>%
  mutate(c = obj$seurat_clusters, ct = obj$celltype, sample = obj$orig.ident)

p_list <- lapply(p, function(x) {
  dd <- x$data %>%
    as_tibble(rownames = "cell") %>%
    left_join(pd)
  t <- x$labels$title

  ggplot(dd, aes(x = umap_1, y = umap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
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

pdf("./figures/tillering_node_rhizome_fig4s4/markers_umap_tillering.pdf", width = 12, height = 16.6)
pp
dev.off()
