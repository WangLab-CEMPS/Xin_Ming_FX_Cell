library(Seurat)
library(dplyr)

obj <- readRDS("./dataLib/osRoot/ggm/ggm2_TQ2_strata.Rds")

table(obj$orig.ident)

obj_ggm <- subset(obj, subset = orig.ident == "ggm2")
obj_tq <- subset(obj, subset = orig.ident == "osRoot2")

# remove cluster cells number < 50
obj_ggm <- subset(
  obj_ggm,
  subset = seurat_clusters %in% names(table(obj_ggm$seurat_clusters))[table(obj_ggm$seurat_clusters) >= 50]
)
obj_tq <- subset(
  obj_tq,
  subset = seurat_clusters %in% names(table(obj_tq$seurat_clusters))[table(obj_tq$seurat_clusters) >= 50]
)

markers_ggm <- FindAllMarkers(
  object = obj_ggm,
  only.pos = TRUE,
)
markers_ggm <- markers_ggm[markers_ggm$pct.1 >= 0.2, ]
table(markers_ggm$cluster)

markers_tq <- FindAllMarkers(
  object = obj_tq,
  only.pos = TRUE,
)
markers_tq <- markers_tq[markers_tq$pct.1 >= 0.2, ]
table(markers_tq$cluster)

# Convert markers to list of gene sets by cluster
markers_ggm_by_cluster <- split(markers_ggm$gene, markers_ggm$cluster)
markers_tq_by_cluster <- split(markers_tq$gene, markers_tq$cluster)

# Function to calculate Overlap coefficient
calculate_overlap <- function(set1, set2) {
  intersection_size <- length(intersect(set1, set2))
  min_size <- min(length(set1), length(set2))
  return(intersection_size / min_size)
}

# Calculate coefficients for shared clusters
shared_clusters <- intersect(names(markers_ggm_by_cluster), names(markers_tq_by_cluster))
overlap_results <- data.frame(
  cluster = shared_clusters,
  overlap_coefficient = sapply(shared_clusters, function(cluster) {
    calculate_overlap(markers_ggm_by_cluster[[cluster]], markers_tq_by_cluster[[cluster]])
  })
)

# Add set sizes information
overlap_results$intersection_size <- sapply(shared_clusters, function(cluster) {
  length(intersect(markers_ggm_by_cluster[[cluster]], markers_tq_by_cluster[[cluster]]))
})
overlap_results$union_size <- sapply(shared_clusters, function(cluster) {
  length(union(markers_ggm_by_cluster[[cluster]], markers_tq_by_cluster[[cluster]]))
})
overlap_results$ggm_set_size <- sapply(shared_clusters, function(cluster) {
  length(markers_ggm_by_cluster[[cluster]])
})
overlap_results$tq_set_size <- sapply(shared_clusters, function(cluster) {
  length(markers_tq_by_cluster[[cluster]])
})

print("\nDetailed results including set sizes and all similarity coefficients:")
print(overlap_results)

print("\nInterpretation of coefficients:")
print("- Overlap coefficient (S): Similarity relative to the smaller set")

# GO -----
library(clusterProfiler)
library(org.At.tair.db)
library(pigRutils)
library(tidyr)

alias <- read.csv("~/ref/Oryza_sativa/annotation/IRGSP_OsAt_Orthologous.tsv", sep = "\t")
# remove duplicates according to RAP_ID
alias <- distinct(alias, RAP_ID, .keep_all = TRUE)
alias <- alias[!is.na(alias$TAIR), ]
dim(alias)

m_ggm <- markers_ggm %>%
  nest(data = everything(), .by = cluster) %>%
  rowwise() %>%
  mutate(tair = list(unique(
    alias$TAIR[match(data$gene, alias$RAP_ID)]
  )))

names(m_ggm$tair) <- m_ggm$cluster

m_tq <- markers_tq %>%
  nest(data = everything(), .by = cluster) %>%
  rowwise() %>%
  mutate(tair = list(unique(
    alias$TAIR[match(data$gene, alias$RAP_ID)]
  )))

names(m_tq$tair) <- m_tq$cluster


# do compareCluster of each shared cluster
comp_go_by_cluster <- lapply(shared_clusters, function(cluster) {
  comp_go <- compareCluster(
    geneClusters = list(ggm = m_ggm$tair[[cluster]], tq = m_tq$tair[[cluster]]),
    fun = "enrichGO",
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "BP"
  )
  return(comp_go)
})

names(comp_go_by_cluster) <- shared_clusters

overlap_coefficient <- lapply(comp_go_by_cluster, function(x) {
  tmp <- as.data.frame(x) %>%
    nest(data = everything(), .by = Cluster)
  names(tmp$data) <- tmp$Cluster
  overlap_coefficient <- calculate_overlap(tmp$data[["ggm"]]$ID, tmp$data[["tq"]]$ID)
  return(overlap_coefficient)
}) %>%
  unlist()


overlap_results$go_overlap_coefficient <- overlap_coefficient


# barplot
library(ggplot2)
library(patchwork)

pd <- overlap_results[, c("overlap_coefficient", "go_overlap_coefficient")]
pd <- round(pd, 1)
pd$cluster <- rownames(pd)
pd$cluster <- factor(pd$cluster, levels = pd$cluster)

p1 <- ggplot(pd, aes(y = cluster, x = overlap_coefficient)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0.4, linetype = "dashed")

p2 <- ggplot(pd, aes(y = cluster, x = go_overlap_coefficient)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0.4, linetype = "dashed")


pdf("./figures/ggm2_tq2_clusters_overlap_coefficient.pdf", width = 5, height = 5)
(p1 + p2) &
  scale_x_continuous(expand = c(0, 0), breaks = c(0.2, 0.4, 0.6, 0.8)) &
  ggprism::theme_prism(base_size = 12) &
  theme()
dev.off()
