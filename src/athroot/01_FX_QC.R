library(Seurat)
library(magrittr)
library(pigRutils)

samples <- c(
  atroot1 = "./rawdata/atRoot/AtR1/AtR1F_Result/outs/filtered_feature_bc_matrix",
  atroot2 = "./rawdata/atRoot/AtR2/AtR2F_Result/outs/filtered_feature_bc_matrix"
)

sample <- samples[2]

(
  obj <- Read10X(data.dir = sample) %>%
    CreateSeuratObject(.,
      min.cells = 2,
      min.features = 200,
      project = names(sample)
    )
)

obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^ATMG")
obj[["percent.chlo"]] <- PercentageFeatureSet(obj, pattern = "^ATCG")

p <- qc_plot(obj, max_nFeature = 6000, min_nFeature = 500, max_nCount = 40000)

ggplot2::ggsave(paste0("./Plots/", names(sample), "_qc.pdf"), p, width = 10, height = 8)

# 定义过滤条件
filter_conditions <- obj$nCount_RNA < 40000 &
  obj$nCount_RNA > 599 &
  obj$nFeature_RNA < 6000 &
  obj$nFeature_RNA > 599 &
  obj$percent.mito < 10 &
  obj$percent.chlo < 5

# 查看过滤统计
table(filter_conditions)

# 应用过滤条件
obj <- subset(
  obj,
  cells = which(filter_conditions)
)

mean(obj$nFeature_RNA)
median(obj$nFeature_RNA)

saveRDS(obj, paste0("./dataLib/athRoot/", names(sample), ".rds"))
