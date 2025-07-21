library(readxl)
library(dplyr)
library(eulerr)
library(Seurat)

micro <- read.csv("./dataLib/birnbaum.csv")
head(micro)
micro <- micro[micro$PROTOPLASTING.EFFECT %in% c("LOST", "INDUCED"), ]
dim(micro)
micro <- micro[micro$Gene.ID %in% rownames(obj), ]


fp <- "./dataLib/athRoot/FX_TQ_protolasted_only_pos_TRUE.xlsx"
sheets <- excel_sheets(fp)
sheets <- sheets[1:12]

de_data <- lapply(sheets, function(x) {
  read_excel(fp, sheet = x) %>%
    subset(avg_log2FC >= 2) %>%
    pull(gene)
})

de_data <- unique(unlist(de_data))


eu <- euler(
  list(
    "Benfy_INDUCED" = unique(micro$Gene.ID[micro$PROTOPLASTING.EFFECT == "INDUCED"]),
    "FX-Cell" = de_data
  )
)

eu <- euler(
  c("Benfy_INDUCED" = 95, "FX-cell" = 1500, "Benfy_INDUCED&FX-cell" = 246)
)

plot(
  eu,
  # quantities = list(type = c("percent", "counts")),
  edges = list(col = "grey90", lex = 3),
  fills = c("#ff8b87", "#ffe46e"),
  legend = list(side = "right")
)


obj <- readRDS("./dataLib/athRoot/FX_TQ_v5.Rds")

table(obj$tech)

Idents(obj) <- "tech"

deg <- FindMarkers(
  obj,
  ident.1 = "tq",
  ident.2 = "atroot",
)

deg <- subset(deg, p_val_adj < 0.05 & abs(avg_log2FC) >= 2)
head(deg)

up <- rownames(deg[deg$avg_log2FC > 0, ])
down <- rownames(deg[deg$avg_log2FC < 0, ])


eu <- list(
  "Benfy_INDUCED" = unique(micro$Gene.ID[micro$PROTOPLASTING.EFFECT == "INDUCED"]),
  "FX-Cell" = up
)

plot(euler(eu))

table(micro$Gene.ID %in% rownames(obj))
