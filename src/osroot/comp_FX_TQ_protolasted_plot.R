library(eulerr)
library(readxl)
library(dplyr)

fp <- "./dataLib/osRoot/TQ_FX_protolasted_celltype.xlsx"
sheets <- excel_sheets(fp)

sc <- lapply(sheets, function(x) {
  read_excel(fp, sheet = x) %>% subset(avg_log2FC >= 1.58) %>% pull(gene)
})

sc <- unlist(sc) %>% unique()

sc <- read.csv("./dataLib/osRoot/TQ_vs_FX_protolasted.csv")
dim(sc)
sc <- sc$X[sc$avg_log2FC >= 1.58] %>% unique()

bulk <- read.csv("./dataLib/osRoot/bulk_proto_post_vs_pre_results.csv")
bulk <- bulk$gene_id[bulk$padj < 0.05 & bulk$log2FoldChange >= 2] %>% unique()


pp1 <- plot(
  euler(
    list(
      "sc" = sc,
      "bulk" = bulk
    )
  ),
  quantities = list(type = c("percent", "counts")),
  edges = list(col = "grey90", lex = 3),
  fills = c("#ff8b87", "#ffe46e")
)
pp1
