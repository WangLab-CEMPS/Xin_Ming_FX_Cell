library(ggplot2)
library(dplyr)
library(tidyr)
library(Seurat)

# read data --------------
root1 <- readRDS("./dataLib/athRoot/atroot1.rds")
root2 <- readRDS("./dataLib/athRoot/atroot2.rds")
tq <- readRDS("./dataLib/athRoot/TQ_root.rds")
# stat --------------
pd <- data.frame(
  samples = c("root1", "root2", "tq"),
  mean = c(mean(root1$nFeature_RNA), mean(root2$nFeature_RNA), mean(tq$nFeature_RNA)),
  median = c(median(root1$nFeature_RNA), median(root2$nFeature_RNA), median(tq$nFeature_RNA))
)

pd
# plot --------------
p1 <- ggplot(pd, aes(x = samples, y = mean)) +
  geom_bar(stat = "identity", fill = "#338dac") +
  ggprism::theme_prism(base_size = 11) +
  ggtitle("mean genes per cell")

p2 <- ggplot(pd, aes(x = samples, y = median)) +
  geom_bar(stat = "identity", fill = "#338dac") +
  ggprism::theme_prism(base_size = 11) +
  ggtitle("median genes per cell")

(p1 + p2) & scale_y_continuous(expand = c(0, 0)) & theme(axis.title = element_blank())

ggsave("./figures/athroot/stat_genes_per_cell.pdf", width = 7, height = 5)
