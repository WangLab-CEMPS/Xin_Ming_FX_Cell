library(dplyr)
library(Seurat)
library(ggplot2)

fl <- c(
  "./rawdata/ref_method/refmethod1", # 1
  "./rawdata/ref_method/refmethod2", # 1
  "./rawdata/ref_method/snRNA1", # 2
  "./rawdata/ref_method/snRNA2", # 2
  "./rawdata/athLeaf/scRNA/LC",
  "./rawdata/athLeaf/scRNA/LPZ",
  "./rawdata/osRoot/ggm1", # 3
  "./rawdata/osRoot/ggm2", # 3
  "./rawdata/osRoot/gdm1", # 4
  "./rawdata/osRoot/gdm2", # 4
  "./rawdata/osRoot/dgm1", # 5
  "./rawdata/osRoot/dgm2" # 5
)

nn <- basename(fl)

stat <- lapply(fl, function(x) {
  oo <- Read10X(x) %>% CreateSeuratObject(., min.cells = 0, min.features = 0)

  cc <- c(
    unique(oo$orig.ident),
    length(colnames(oo)),
    length(rownames(oo)),
    mean(oo$nFeature_RNA),
    mean(oo$nCount_RNA),
    median(oo$nFeature_RNA),
    median(oo$nCount_RNA)
  )

  names(cc) <- c("samples", "nCells", "nGenes", "mean_nFeature_RNA", "mean_nCount_RNA", "median_nFeature_RNA", "median_nCount_RNA")
  cc
}) %>%
  do.call(rbind, .) %>%
  as.data.frame()

stat$samples <- nn

stat

write.csv(stat, "./dataLib/refmethod_snRNA_athleaf_osroot_stat_raw.csv", quote = FALSE, row.names = FALSE)

# plot cell number & mean nFeature_RNA
pd <- stat
pd <- pd[!pd$samples %in% c("LC", "LPZ"), ]
pd$group <- gsub("[1-2]$", "", pd$samples)
pd$group <- factor(pd$group, levels = c("refmethod", "snRNA", "ggm", "gdm", "dgm"))
pd <- pd[order(pd$group), ]

pd_ncells <- select(pd, group, samples, nCells)
pd_nFeature <- select(pd, group, samples, mean_nFeature_RNA, median_nFeature_RNA)


p1 <- ggplot(pd_ncells, aes(fill = samples, y = nCells, x = group)) +
  geom_bar(position = "dodge", stat = "identity", color = "grey85") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = rep(c("#D1313A", "#E77E67", "#F2A47E", "#eee1d9", "#E3B194"), each = 2)) +
  ggprism::theme_prism() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
  )


p2 <- ggplot(pd_nFeature, aes(fill = samples, y = mean_nFeature_RNA, x = group)) +
  geom_bar(position = "dodge", stat = "identity", color = "grey85") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = rep(c("#D1313A", "#E77E67", "#F2A47E", "#eee1d9", "#E3B194"), each = 2)) +
  ggprism::theme_prism() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
  )


p3 <- ggplot(pd_nFeature, aes(fill = samples, y = median_nFeature_RNA, x = group)) +
  geom_bar(position = "dodge", stat = "identity", color = "grey85") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = rep(c("#D1313A", "#E77E67", "#F2A47E", "#eee1d9", "#E3B194"), each = 2)) +
  ggprism::theme_prism() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
  )

ggsave("./Plots/nCells_raw_bar.pdf", p1, width = 8, height = 4)
ggsave("./Plots/mean_nFeature_raw_bar.pdf", p2, width = 8, height = 4)
ggsave("./Plots/median_nFeature_raw_bar.pdf", p3, width = 8, height = 4)

