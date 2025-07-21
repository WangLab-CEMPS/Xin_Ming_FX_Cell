library(Seurat)
library(magrittr)
library(dplyr)
library(eulerr)

# load data
obj <- readRDS("./dataLib/athLeaf/scsnRNA_harmony.Rds")
obj <- JoinLayers(obj)

table(obj$orig.ident)
obj_sn <- subset(obj, subset = orig.ident %in% c("snLC", "snLPZ"))
obj_sc <- subset(obj, subset = orig.ident == c("scLC", "scLPZ"))



deg_sn <- FindMarkers(
  object = obj_sn,
  group.by = "orig.ident",
  ident.1 = "snLC",
  ident.2 = "snLPZ"
)

sub_deg_sn <- subset(deg_sn, p_val < 0.05 & p_val_adj < 0.05 & abs(avg_log2FC) > 1) %>%
  mutate(sig = ifelse(avg_log2FC > 0, "up", "down"))

dim(sub_deg_sn)

deg_sc <- FindMarkers(
  object = obj_sc,
  group.by = "orig.ident",
  ident.1 = "scLC",
  ident.2 = "scLPZ"
)

sub_deg_sc <- subset(deg_sc, p_val < 0.05 & p_val_adj < 0.05 & abs(avg_log2FC) > 1) %>%
  mutate(sig = ifelse(avg_log2FC > 0, "up", "down"))

dim(sub_deg_sc)


p1 <- plot(
  euler(
    list(
      "sc_all" = rownames(sub_deg_sc),
      "sn_all" = rownames(sub_deg_sn)
    )
  ),
  quantities = TRUE,
  labels = list(fontsize = 12),
  edges = list(col = "white", lex = 2),
  fill = c("#95d3f7", "#fad88f"),
)

p2 <- plot(
  euler(
    list(
      "sc_LC_up" = rownames(sub_deg_sc[sub_deg_sc$sig == "up", ]),
      "sn_LC_up" = rownames(sub_deg_sn[sub_deg_sn$sig == "up", ])
    )
  ),
  quantities = TRUE,
  labels = list(fontsize = 12),
  edges = list(col = "white", lex = 2),
  fill = c("#95d3f7", "#fad88f"),
)

p3 <- plot(
  euler(
    list(
      "sc_LPZ_up" = rownames(sub_deg_sc[sub_deg_sc$sig == "down", ]),
      "sn_LPZ_up" = rownames(sub_deg_sn[sub_deg_sn$sig == "down", ])
    )
  ),
  quantities = TRUE,
  labels = list(fontsize = 12),
  edges = list(col = "white", lex = 2),
  fill = c("#95d3f7", "#fad88f")
)


pdf("./Plots/venn_sn_sc_1.pdf", width = 5, height = 4)
print(p1)
print(p2)
print(p3)
dev.off()