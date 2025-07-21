library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(patchwork)

# load data
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public_anno.Rds")
obj <- JoinLayers(obj)
obj <- subset(obj, subset = orig.ident %in% c("LC", "LPZ"))
obj$celltype <- gsub("Mesophyll\\*", "Mesophyll", obj$celltype)
sel_ct <- c("Phloem", "Bundle-Sheet", "Epidermis", "Hydathode", "Mesophyll")
obj <- subset(obj, subset = celltype %in% sel_ct)

cell_types <- c(
  "epi" = "Epidermis",
  "meso" = "Mesophyll",
  "bs" = "Bundle-Sheet",
  "ph" = "Phloem",
  "hyd" = "Hydathode"
)

epi <- c("AT1G68870", "AT3G48100", "AT5G10250", "AT3G04110", "AT5G58670")
bs <- c("AT3G25710", "AT4G36470", "AT5G61480", "AT1G17745", "AT2G32150", "AT4G27070")
hyd <- c("AT1G28330", "AT2G28950", "AT2G36830", "AT5G48930")
meso <- c("AT1G17140", "AT4G14130", "AT5G01240", "AT1G30040", "AT2G25735", "AT3G19270")
ph <- c("AT1G68740", "AT3G54720", "AT1G78850", "AT2G39010")

g_list <- list(epi = epi, bs = bs, hyd = hyd, meso = meso, ph = ph)

for (i in seq_along(g_list)) {
  ct <- names(g_list)[i]
  print(ct)
  sub_obj <- subset(obj, subset = celltype == cell_types[ct])
  meta <- sub_obj@meta.data %>%
    as_tibble(rownames = "cell") %>%
    select(cell, orig.ident)
  p <- FeaturePlot(sub_obj, features = g_list[[ct]], reduction = "harmony_umap", order = TRUE)

  sapply(seq.int(p), function(x) {
    xx <- p[[x]]
    pd <- xx$data %>%
      as_tibble(rownames = "cell") %>%
      left_join(meta)
    t <- xx$labels$title
    print(t)
    p <- ggplot(pd, aes(x = harmonyumap_1, y = harmonyumap_2, color = !!sym(t))) +
      geom_point_rast(size = 0.2, alpha = 0.9) +
      facet_wrap(~orig.ident) +
      scale_color_gradientn(colours = c("#c8c6c3", "#f6ee6c", "#ffa42c", "#eb212f", "#a60008")) +
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

    ggsave(paste0("./figures/leaf/select_ct_umap/", t, ".pdf"), p, height = 4, width = 7)
  })
}


small_obj <- subset(obj, subset = celltype == "Phloem")
meta <- small_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(cell, orig.ident)

p <- FeaturePlot(small_obj, features = "AT2G39010", reduction = "harmony_umap", order = TRUE)

pd <- p$data %>%
  as_tibble(rownames = "cell") %>%
  left_join(meta)
t <- p$labels$title
print(t)
p <- ggplot(pd, aes(x = harmonyumap_1, y = harmonyumap_2, color = !!sym(t))) +
  geom_point_rast(size = 2, alpha = 0.9) +
  facet_wrap(~orig.ident) +
  scale_color_gradientn(colours = c("#c8c6c3", "#f6ee6c", "#ffa42c", "#eb212f", "#a60008")) +
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

ggsave("./figures/leaf/select_ct_umap/AT2G39010_phloem.pdf", p, height = 4, width = 7)
