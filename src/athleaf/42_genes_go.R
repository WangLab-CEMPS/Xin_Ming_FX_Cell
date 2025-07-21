library(org.At.tair.db)
library(clusterProfiler)
library(dplyr)

genes <- read.csv("./dataLib/athLeaf/intersect_all_42.csv")[[1]]

go <- enrichGO(
  gene = genes,
  OrgDb = org.At.tair.db,
  ont = "BP",
  keyType = "TAIR"
)

go <- clusterProfiler::simplify(
  go,
  cutoff = 0.7,
  by = "pvalue",
  select_fun = min
)

pdf("./figures/go_42_genes.pdf", width = 7, height = 4)
barplot(go, showCategory = 15, label_format = 80)
dev.off()

as.data.frame(go) %>% write.csv("./dataLib/athLeaf/gores_42_genes.csv", row.names = FALSE)
