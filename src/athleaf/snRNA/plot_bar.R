library(Seurat)
library(ggplot2)
library(dplyr)


snlc <- readRDS("./dataLib/athLeaf/snRNA/snLC.rds")
snlpz <- readRDS("./dataLib/athLeaf/snRNA/snLPZ.rds")
lc <- readRDS("./dataLib/athLeaf/LC.rds")
lpz <- readRDS("./dataLib/athLeaf/LPZ.rds")
# stat n_Feature_RNA & n_Count

pd <- data.frame(
  sample = c("snLC", "snLPZ", "scLC", "scLPZ"),
  group = c("snRNA", "snRNA", "scRNA", "scRNA"),
  type = c("LC", "LPZ", "LC", "LPZ"),
  cell_num = c(
    ncol(snlc), ncol(snlpz), ncol(lc), ncol(lpz)
  ),
  mean_nFeature_RNA = c(
    mean(snlc$nFeature_RNA), mean(snlpz$nFeature_RNA), mean(lc$nFeature_RNA), mean(lpz$nFeature_RNA)
  ),
  mean_nCount = c(
    mean(snlc$nCount_RNA), mean(snlpz$nCount_RNA), mean(lc$nCount_RNA), mean(lpz$nCount_RNA)
  ),
  median_nFeature_RNA = c(
    median(snlc$nFeature_RNA), median(snlpz$nFeature_RNA), median(lc$nFeature_RNA), median(lpz$nFeature_RNA)
  ),
  median_nCount = c(
    median(snlc$nCount_RNA), median(snlpz$nCount_RNA), median(lc$nCount_RNA), median(lpz$nCount_RNA)
  )
)

write.csv(pd, "./dataLib/athLeaf/snRNA/snRNA_stat.csv")


p <- ggplot(pd, aes(x = sample, y = mean_nFeature_RNA, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(mean_nFeature_RNA, 0)),
    vjust = 1.6, size = 11
  ) +
  scale_fill_manual(values = c("#2cb5a2", "#5fa1d8")) +
  ggprism::theme_prism(base_size = 11) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Mean nFeature_RNA")

p

ggsave("./figures/athleaf/snscRNA_mean_nFeature_RNA.pdf", p, width = 7, height = 7)

p2 <- ggplot(pd, aes(x = sample, y = median_nFeature_RNA, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(median_nFeature_RNA, 0)),
    vjust = 1.6, size = 11
  ) +
  scale_fill_manual(values = c("#2cb5a2", "#5fa1d8")) +
  ggprism::theme_prism(base_size = 11) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Median nFeature_RNA")

p2

ggsave("./figures/athleaf/snscRNA_median_nFeature_RNA.pdf", p2, width = 7, height = 7)
