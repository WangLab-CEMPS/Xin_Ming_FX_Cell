library(Seurat)
library(tidyr)
library(dplyr)

obj <- readRDS("./dataLib/tillering_node_rhizome/rhizome_S0.rds")

markers <- FindAllMarkers(
  object = obj,
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 1.32
)

table(markers$cluster)

annoos <- read.csv("./dataLib/ls_os_anno.csv")
head(annoos)
annoat <- read.csv("./dataLib/ol_to_at.csv", row.names = 1)
head(annoat)

markers <- left_join(markers, annoos, by = c("gene" = "qseqid"))
markers <- left_join(markers, annoat, by = c("gene" = "ol_gene"))
head(markers)

go <- markers %>%
  nest(data = everything(), .by = cluster) %>%
  rowwise() %>%
  mutate(tair = list(data$ath_gene))

names(go$tair) <- paste0("Cluster", go$cluster)
names(go$data) <- paste0("Cluster", go$cluster)

openxlsx::write.xlsx(go$data, "./dataLib/tillering_node_rhizome/rhizome_S0_markers.xlsx")

# GO ---------
library(pigRutils)
library(org.At.tair.db)
library(clusterProfiler)

orgdb <- org.At.tair.db

go_res <- batch_GO(go$tair, orgdb = orgdb, keyType = "TAIR")

export_GO(go_res, "./dataLib/tillering_node_rhizome/rhizome_S0_markers_GO.xlsx")
