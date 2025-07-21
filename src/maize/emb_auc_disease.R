library(AUCell)
library(Seurat)
library(scales)
library(dplyr)
library(ggplot2)
library(ggrastr)

obj <- readRDS("./dataLib/maize/maize_pub.Rds")

genes4 <- read.table("./dataLib/maize/supp.txt")[[1]]
anno <- read.delim("./dataLib/maize/genes_all.txt")
anno <- anno[anno$v4.Gene.Model.ID %in% genes4, ]

sel_genes <- unique(anno$v5.Gene.Model.ID)
sel_genes <- sel_genes[sel_genes %in% rownames(obj)]


obj <- JoinLayers(obj)
Layers(obj)
cts <- LayerData(obj, layer = "counts")

geneSets <- list(
  dis = sel_genes
)

cells_AUC <- AUCell_run(cts, geneSets)

pd <- Embeddings(obj, reduction = "harmony_umap") %>% as.data.frame()
pd$auc <- pigRutils::mf_minmax(getAUC(cells_AUC)["dis", ])
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
  )


ggsave("./Plots/maize/maize_disease_auc.pdf", height = 6.2, width = 18)


p <- FeaturePlot(obj, sel_genes)

meta <- obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(cell, orig.ident)

p_list <- lapply(p, function(x) {
  pd <- x$data %>%
    as_tibble(rownames = "cell") %>%
    left_join(meta)
  t <- x$labels$title

  p <- ggplot(pd, aes(x = harmonyumap_1, y = harmonyumap_2, color = !!sym(t))) +
    geom_point_rast(size = 0.1) +
    facet_wrap(~orig.ident) +
    ggtitle(t)
  ggsave(paste0("./tmp/", t, ".png"), height = 4, width = 8)
})
