library(ggplot2)
library(dplyr)
library(ggrepel)

alias <- read.csv("~/ref/Athaliana/Araport11/gene_aliases_20220331_tidy.csv")
deg <- read.csv("./dataLib/athLeaf/LPZ_Mesophyll_de.csv")
head(deg)

label_genes <- c(
  "AT5G43080" = "CYCA3;1",
  "AT1G02970" = "ATWEE1",
  "AT1G04020" = "ATBARD1",
  "AT5G08020" = "ATRPA70B"
)

label_genes2 <- c(
  "AT1G02930" = "ATGST1",
  "AT1G72520" = "ATLOX4",
  "AT1G73080" = "ATPEPR1",
  "AT1G02400" = "ATGA2OX4",
  "AT1G02920" = "ATGST11",
  "AT3G23250" = "ATMYB15",
  "AT3G28210" = "PMZ"
)

deg$p_val_adj[deg$p_val_adj == 0] <- 10^(-323)

label_dat <- deg %>%
  filter(geneId %in% names(c(label_genes, label_genes2))) %>%
  left_join(alias, by = c("geneId" = "name"))
label_dat <- label_dat[, -c(8, 9)]

deg <- deg %>%
  mutate(
    Change = case_when(
      avg_log2FC > 2 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < -2 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "Not"
    )
  )

p <- ggplot() +
  ggrastr::geom_point_rast(
    data = deg,
    aes(
      x = avg_log2FC,
      y = -log10(p_val_adj),
      color = Change,
      size = abs(avg_log2FC),
      alpha = 0.7
    )
  ) +
  scale_size(range = c(0.5, 4)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_label_repel(
    data = label_dat,
    aes(x = avg_log2FC, y = -log10(p_val_adj), label = symbol),
    size = 4,
    fontface = 2,
    segment.size = 1,
    box.padding = 0.7
  ) +
  geom_point(
    data = label_dat,
    aes(
      x = avg_log2FC,
      y = -log10(p_val_adj),
      color = "Labeled Gene",
      size = abs(avg_log2FC)
    )
  ) +
  scale_y_continuous(expand = c(0, 8), limits = c(0, 340)) +
  theme_classic() +
  theme(
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  scale_color_manual(
    values = c("#2f80d6", "#e05328", "#888888", "black"),
    breaks = c("Down", "Up", "Not", "Labeled Gene")
  ) +
  labs(x = "Log2FC", y = "-Log10(Adjusted Pvalue)") +
  guides(
    size = guide_legend(title = "|Log2FC|"),
    color = guide_legend(title = "Sig"),
    alpha = "none"
  )

pdf("./figures/leaf/LPZ_Mesophyll_de_volcano.pdf", width = 7, height = 9)
p
dev.off()
