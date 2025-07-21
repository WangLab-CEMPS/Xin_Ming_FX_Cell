library(Seurat)
library(eulerr)
library(dplyr)
library(readxl)
library(clusterProfiler)
library(org.At.tair.db)
library(ggplot2)

# read data
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public_anno.Rds")
detected_genes <- rownames(obj)
rm(obj)

dp <- "./dataLib/athLeaf/PC_Wolf_B_Frommer_protoplasted.xlsx"
sheet <- excel_sheets(dp)
bulk_protolast <- read_excel(dp, sheet = sheet[4])
colnames(bulk_protolast) <- gsub(" ", "_", colnames(bulk_protolast))
head(bulk_protolast)

sc_protolast <- read.csv("./dataLib/athLeaf/sc_protolasted_leaf.csv", row.names = 1)
head(sc_protolast)

# filter
sc_protolast <- sc_protolast %>%
  subset(abs(avg_log2FC) > 2 & p_val_adj < 0.05)

dim(sc_protolast)

bulk_protolast <- bulk_protolast[bulk_protolast$Gene_ID %in% detected_genes, ] %>%
  subset(abs(Log2_FC) > 2 & padj < 0.05)

dim(bulk_protolast)
k <- 2
bulk_protolast <- bulk_protolast[rowSums(bulk_protolast[, c(4, 6, 8, 10)] >= 1) >= k, ]
dim(bulk_protolast)

# overlap
intersect(rownames(sc_protolast), bulk_protolast$Gene_ID) %>% length()

# venn diagram
gl <- list(
  SC_protolasted = rownames(sc_protolast),
  Bulk_protolasted = bulk_protolast$Gene_ID
)

fit <- euler(gl)

p <- plot(
  fit,
  edges = list(col = "lightgray", lex = 3),
  labels = list(font = 1.2, cex = 1.2),
  quantities = list(col = "black", font = 1.5),
  fill = c("#c9d9f4", "#e7ffd2")
)

pdf("./figures/venn_sc_bulk_protolasted_new.pdf", width = 7, height = 7)
p
dev.off()
