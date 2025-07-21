library(dplyr)
library(magrittr)
library(sampling)
library(Seurat)
set.seed(131415)

obj <- readRDS("./dataLib/osRoot/osRoots_TQ_v5.rds")
obj <- subset(obj, subset = orig.ident != "10X_osRoot1_singlet")

ratio <- 0.51


table(obj$seurat_clusters)

cluster_dat <- data.frame(
  Cluster = obj$seurat_clusters,
  cell = names(obj$seurat_clusters)
)

cluster_dat <- cluster_dat[gtools::mixedorder(cluster_dat$Cluster), ]


size <- cluster_dat %>%
  group_by(Cluster) %>%
  summarise(size = n() * ratio) %>%
  .$size

stra_sample <- strata(
  cluster_dat,
  stratanames = "Cluster",
  size = round(size),
  method = "srswor"
)

dat <- getdata(cluster_dat, stra_sample)
# 保存分层抽样数据
write.csv(dat, "./dataLib/osRoot/root_tq_strata_sampling_51.csv")

obj <- subset(obj, cells = dat$cell)
# 构建Seurat对象
obj %<>% GetAssayData(layer = "counts") %>% CreateSeuratObject()

saveRDS(obj, "./dataLib/osRoot/root_tq_strata.rds")
