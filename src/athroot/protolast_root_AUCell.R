library(dplyr)
library(AUCell)
library(readxl)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(scales)
library(patchwork)
library(ggprism)


obj <- readRDS("./dataLib/athRoot/FX_TQ_v5.Rds")

proto_file <- "./dataLib/athRoot/2019_DC_athroot_Tom_Denyer.xlsx"
sheets <- excel_sheets(proto_file)
sheets

proto <- read_excel(proto_file, sheet = sheets[3])
colnames(proto) <- gsub(" ", "_", colnames(proto))

proto_genes_up <- proto$Gene_ID[proto$Log2_FC >= 2]
length(proto_genes_up)

proto_genes_down <- proto$Gene_ID[proto$Log2_FC <= -2]
length(proto_genes_down)

cts <- GetAssayData(obj, layer = "counts")

cells_AUC <- AUCell_run(
  cts,
  list(
    pp_up = proto_genes_up,
    pp_down = proto_genes_down,
    pp_all = c(proto_genes_up, proto_genes_down)
  )
)
obj$pp_up <- getAUC(cells_AUC)["pp_up", ]
obj$pp_down <- getAUC(cells_AUC)["pp_down", ]
obj$pp_all <- getAUC(cells_AUC)["pp_all", ]

FeaturePlot(obj, features = "pp_up", split.by = "tech")
FeaturePlot(obj, features = "pp_down", split.by = "tech")
FeaturePlot(obj, features = "pp_all", split.by = "tech")

saveRDS(obj, "./dataLib/athRoot/FX_TQ_v5.Rds")

# plot ----
# umap split by tech
pd <- Embeddings(obj, reduction = "umap") %>%
  as.data.frame() %>%
  mutate(cells = rownames(.)) %>%
  left_join(obj@meta.data %>% mutate(cells = rownames(.)), by = "cells")

pd$tech <- ifelse(pd$tech == "tq", "TQ", "FX-Cell")

head(pd)

p1 <- ggplot(pd, aes(x = UMAP_1, y = UMAP_2, color = pp_all)) +
  ggrastr::geom_point_rast(size = 0.2) +
  facet_wrap(~tech) +
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) +
  theme_linedraw(base_size = 11) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  ggtitle("AUCell Score of Protoplastization Affects Gene Sets, from Tom Denyer et al. (DC, 2019)")

ggsave("./figures/athroot/FX_TQ_protolast_AUCell.pdf", p1, width = 10, height = 5)

pd %>%
  group_by(tech) %>%
  summarise(mean(pp_all))

quantile(pd$pp_all, probs = seq(0, 1, 0.1))

pd$celltype <- factor(
  pd$celltype,
  levels = pd %>%
    subset(tech == "TQ") %>%
    group_by(celltype) %>%
    summarise(mean_pp_all = mean(pp_all)) %>%
    arrange(desc(mean_pp_all)) %>%
    pull(celltype)
)

# violin plot of celltype split by tech
p2 <- ggplot(pd, aes(x = celltype, y = pp_all)) +
  geom_violin(aes(fill = tech), width = 1, position = position_dodge(width = 0.7)) +
  geom_boxplot(aes(fill = tech), width = 0.1, position = position_dodge(width = 0.7)) +
  geom_hline(aes(yintercept = 0.18), linetype = "dashed", color = "#ecc61a", linewidth = 0.2) +
  geom_hline(aes(yintercept = 0.05), linetype = "dashed", color = "#ff8b87", linewidth = 0.2) +
  # facet_wrap(~tech) +
  ylab("AUCell Score") +
  xlab("CellType") +
  theme_prism(base_size = 11) +
  scale_fill_manual(values = c("#ff8b87", "#ecc61a")) +
  scale_x_discrete(
    guide = guide_prism_bracket()
  ) +
  theme(
    legend.title = element_blank(), axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.text = element_text(face = "plain")
  )

ggsave("./figures/athroot/FX_TQ_protolast_violin.pdf", p2, width = 16, height = 2)


pdf("./figures/athroot/FX_TQ_protolast_AUCell_umap_violin.pdf", width = 17, height = 11)
(p1 / p2) + plot_layout(heights = c(5, 2))
dev.off()
