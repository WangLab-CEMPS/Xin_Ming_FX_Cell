library(eulerr)
library(dplyr)
library(tidyr)

scfl <- list.files("./dataLib/athLeaf/de", pattern = "de.csv$", full.names = TRUE)
snfl <- list.files("./dataLib/athLeaf/snRNA/de", pattern = "de.csv$", full.names = TRUE)


sc <- lapply(scfl, function(x) {
  read.csv(x) %>%
    subset(p_val < 0.05 & p_val_adj < 0.05 & avg_log2FC < -1.58) %>%
    pull(geneId)
}) %>%
  unlist() %>%
  unique()


sn <- lapply(snfl, function(x) {
  read.csv(x) %>%
    subset(p_val < 0.05 & p_val_adj < 0.05 & avg_log2FC < -1.58) %>%
    pull(geneId)
}) %>%
  unlist() %>%
  unique()

p1 <- plot(
  euler(
    list(
      sc = sc,
      sn = sn
    )
  ),
  quantities = TRUE,
  labels = list(fontsize = 12),
  edges = list(col = "white", lex = 2),
  fill = c("#95d3f7", "#fad88f")
)

pdf("./Plots/venn_LPZ_up.pdf", width = 4, height = 3)
print(p1)
dev.off()
