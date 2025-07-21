library(AUCell)
library(Seurat)
library(scales)
library(dplyr)
library(ggplot2)
library(pigRutils)

go <- read.csv("./dataLib/athLeaf/GO_response_to_wounding.tsv", sep = "\t", header = FALSE)
go <- go[[1]]

obj <- readRDS("./dataLib/athLeaf/scsnRNA_harmony.Rds")
obj <- JoinLayers(obj)
Layers(obj)
cts <- LayerData(obj, layer = "counts")

GeneSets <- list(
  wound = unique(go)
)

CellsAUC <- AUCell_run(cts, GeneSets)

pd <- Embeddings(obj, reduction = "harmony_umap") %>% as.data.frame()
pd$auc <- pigRutils::mf_minmax(getAUC(CellsAUC)["wound", ])
pd$cluster <- obj$seurat_clusters
pd$celltype <- obj$celltype
pd$group <- obj$orig.ident

head(pd)

p <- ggplot(pd, aes(x = harmonyumap_1, y = harmonyumap_2, color = auc)) +
  ggrastr::geom_point_rast(size = 0.2) +
  facet_wrap(~group, nrow = 1) +
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

ggsave("./figures/leaf/snRNA/goID_wound_AUCell_umap.pdf", p, height = 4.5, width = 20)

table(pd$celltype)
pd$celltype <- gsub("Mesophyll\\*", "Mesophyll", pd$celltype)
table(pd$celltype)

pd$group <- factor(pd$group, levels = c("Public_data", "scLC", "scLPZ", "snLC", "snLPZ"))

head(pd)
p3 <- ggplot(pd, aes(x = cluster, y = auc)) +
  geom_violin(aes(fill = group), drop = FALSE) +
  geom_boxplot(width = 0.1) +
  geom_hline(aes(yintercept = round(mean(auc), 1)), linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~group, nrow = 1) +
  ggtitle("AUCell Score: Response to wounding") +
  ylab("AUCell Score") +
  xlab("Cluster") +
  ggprism::theme_prism(base_size = 11) +
  # scale_fill_manual(values = c("#f3827e", "#22bba5", "#59a1d7")) +
  coord_flip() +
  theme(
    legend.title = element_blank(), axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.text = element_text(face = "plain")
  )

pdf("./figures/leaf/snRNA/goID_wound_AUCell_violin.pdf", width = 18, height = 5)
p3
dev.off()

cc <- c("#e5dcca", "#e79ade", "#f786a8", "#60cae4", "#2394c8")

median(pd$auc)

p4 <- ggplot(pd, aes(x = celltype, y = auc)) +
  geom_violin(aes(fill = group), drop = FALSE) +
  geom_boxplot(width = 0.1) +
  geom_hline(aes(yintercept = 0.26), linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~group, nrow = 1) +
  ggtitle("AUCell Score: Response to wounding") +
  ylab("AUCell Score") +
  xlab("Cluster") +
  ggprism::theme_prism(base_size = 11) +
  scale_fill_manual(values = cc) +
  coord_flip() +
  theme(
    legend.title = element_blank(), axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.text = element_text(face = "plain")
  )

pdf("./figures/leaf/snRNA/goID_wound_AUCell_violin_celltype.pdf", width = 18, height = 5)
p4
dev.off()


pp <- p <- embSCdim(
  obj, "harmony_umap",
  group_by = "orig.ident",
  colors = cc,
  save = TRUE,
  width = 7.2, height = 5.5,
  save_path = "./figures/leaf/snRNA/scsnRNA_harmony_umap_with_samples.pdf"
)
