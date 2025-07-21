library(eulerr)
library(dplyr)
library(readxl)

# read in data
# pseudo bulk vs L24
pseudo_bulk_L24 <- read.table("./dataLib/proto/pseudo_bulk_vs_bulk24_DE.txt", sep = "\t")
pseudo_bulk_L24 <- subset(pseudo_bulk_L24, padj < 0.05 & abs(log2FoldChange) >= 2)
dim(pseudo_bulk_L24)


# sc
sc <- read.csv("./dataLib/athLeaf/sc_protolasted_leaf.csv")
sc <- subset(sc, p_val_adj < 0.05 & avg_log2FC >= 2)
dim(sc)

sc_2 <- "./dataLib/athLeaf/JC_vs_LC_protolasted.xlsx"
sheet <- excel_sheets(sc_2)
sc_2 <- lapply(sheet, function(x) {
  df <- read_excel(sc_2, sheet = x)
  df <- subset(df, p_val_adj < 0.05 & avg_log2FC >= 2)
  df
})
names(sc_2) <- sheet

sc_2 <- do.call(rbind, sc_2)
dim(sc_2)


# bulk
bulk17 <- read.csv("./dataLib/athLeaf/bulk_proto/L17_vs_P17_DE.csv")
bulk17 <- subset(bulk17, padj < 0.05 & log2FoldChange < -2)
dim(bulk17)

bulk24 <- read.csv("./dataLib/athLeaf/bulk_proto/L24_vs_P24_DE.csv")
bulk24 <- subset(bulk24, padj < 0.05 & log2FoldChange < -2)
dim(bulk24)

bulk_L24_vs_L17 <- read.csv("./dataLib/athLeaf/bulk_proto/L24d_VS_L17d_DE.csv")
bulk_L24_vs_L17 <- subset(bulk_L24_vs_L17, padj < 0.05 & abs(log2FoldChange) >= 2)
dim(bulk_L24_vs_L17)

bulk_P24_vs_P17 <- read.csv("./dataLib/athLeaf/bulk_proto/P24_VS_P17_DE.csv")
bulk_P24_vs_P17 <- subset(bulk_P24_vs_P17, padj < 0.05 & abs(log2FoldChange) >= 2)
dim(bulk_P24_vs_P17)

bulk_P17_vs_L24d <- read.csv("./dataLib/athLeaf/bulk_proto/P17_VS_L24d_DE.csv")
bulk_P17_vs_L24d <- subset(bulk_P17_vs_L24d, padj < 0.05 & log2FoldChange >= 2)
dim(bulk_P17_vs_L24d)

dp <- "./dataLib/athLeaf/PC_Wolf_B_Frommer_protoplasted.xlsx"
sheet <- excel_sheets(dp)
bulk_pub <- read_excel(dp, sheet = sheet[4])
colnames(bulk_pub) <- gsub(" ", "_", colnames(bulk_pub))
bulk_pub <- subset(bulk_pub, padj < 0.05 & Log2_FC >= 2)
dim(bulk_pub)

set_remove <- c(bulk_L24_vs_L17$geneId) %>% unique()
sc_genes <- sc$X[!sc$X %in% set_remove] %>% unique()
bulk_genes <- bulk17$geneId[!bulk17$geneId %in% set_remove] %>% unique()
bulk_p <- bulk_pub$Gene_ID[!bulk_pub$Gene_ID %in% set_remove] %>% unique()

# venn diagram
gl <- list(
  sc = sc$X,
  bulk = bulk_P17_vs_L24d$geneId
)


pdf("./figures/proto/leaf_sc_bulk_P17L24_venn.pdf", width = 8, height = 6)

plot(
  euler(gl),
  quantities = list(type = c("counts")),
  edges = list(col = "white", lex = 1),
  fills = c("#87cdff", "#ffe96e"),
  legend = list(side = "right")
)

dev.off()


gl <- list(
  sc = sc$X,
  bulk_public = bulk_p
)

pdf("./figures/proto/leaf_sc_bulk_pub_venn.pdf", width = 8, height = 6)

plot(
  euler(gl),
  quantities = list(type = c("counts")),
  edges = list(col = "white", lex = 1),
  fills = c("#87cdff", "#e9e8e2"),
  legend = list(side = "right")
)

dev.off()



gl <- list(
  sc = sc_genes,
  bulk_public = bulk_pub$Gene_ID
)


pdf("./figures/proto/leaf_sc_bulk_pub_venn_filtered.pdf", width = 8, height = 6)

plot(
  euler(gl),
  quantities = list(type = c("counts")),
  edges = list(col = "white", lex = 1),
  fills = c("#87cdff", "#e9e8e2"),
  legend = list(side = "right")
)

dev.off()


write.table(sc_genes, "./dataLib/proto/sc_genes.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(bulk_pub$Gene_ID, "./dataLib/proto/bulk_genes.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
