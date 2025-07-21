library(stringr)
library(DESeq2)
library(dplyr)


cts <- read.table("./athLeaf/bulk/LC_LPC.txt", header = TRUE, row.names = 1)
cts <- cts[, -c(1:5)]
head(cts)

colnames(cts) <- gsub("\\.sorted\\.bam$", "", colnames(cts))
colnames(cts)

coldat <- data.frame(
  sample = colnames(cts),
  group = str_split(colnames(cts), "_", simplify = TRUE)[, 1]
)

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldat,
  design = ~group
)

corre <- cor(cts)

col_fun <- circlize::colorRamp2(
  c(min(corre), median(corre), max(corre)),
  c("#3bb0f0", "#ebebeb", "#f95426")
)

library(ComplexHeatmap)

p1 <- Heatmap(
  corre,
  col = col_fun,
  name = "Correlation",
  show_column_names = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  rect_gp = gpar(col = "white")
)


png("./Plots/LCLPZ_cor.png", width = 450, height = 270, res = 100)
p1
dev.off()

dds <- estimateSizeFactors(dds)
vsd <- vst(dds)

library(pigRutils)

df_vsd <- doPCA(indata = assay(vsd), coldata = coldat, pcsToUse = 1:3, intgroup = "group")
percentVar_vsd <- round(100 * attr(df_vsd, "percentVar"))

df <- df_vsd
percentVar <- percentVar_vsd

df2 <- df %>%
  group_by(group) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))
df2$group <- stringr::str_split(df2$group, "[.]", simplify = TRUE)[, 2]

library(ggplot2)

p2 <- ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group), stat = "unique", size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = df2,
    aes(label = group),
    size = 2.5,
    bg.r = 0.25,
    color = "#000103",
    bg.color = "#f5f0f0",
    seed = 1,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(1, "lines")
  ) +
  labs(title = "PCA analysis", caption = "Map to Tak1 V6") +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  coord_fixed()

ggsave(plot = p2, filename = "./Plots/LCLPZ_pca.png", width = 6, height = 3)


dds <- DESeq(dds)

count_norm <- counts(dds, normalize = TRUE)

res <- results(dds) %>% as_tibble(rownames = "geneId") %>% subset(abs(log2FoldChange) > 1)
res

drop_gene <- res[["geneId"]]

write.table(drop_gene, file = "dataLib/athLeaf/drop_gene.txt", row.names = FALSE)
