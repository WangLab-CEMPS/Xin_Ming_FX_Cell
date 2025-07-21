library(Seurat)
library(dplyr)

# 设置
wd <- "/data/wmc_data/wkdir/MingX"
setwd(wd)
output_dir <- "/data/wmc_data/wkdir/MingX/dataLib/proto"
set.seed(123)

# 读取和预处理数据
obj <- readRDS("./dataLib/athLeaf/LC_LPZ_public_anno.Rds")
obj <- subset(obj, subset = orig.ident %in% c("LC"))
obj <- JoinLayers(obj)
obj$celltype <- gsub("Mesophyll\\*", "Mesophyll", obj$celltype)

# 分层抽样
cell_data <- data.frame(cell_id = colnames(obj), celltype = obj$celltype)
celltype_counts <- table(cell_data$celltype)
split_groups <- rep("g roup2", nrow(cell_data))

# 对每个细胞类型进行50%抽样分配给group1
for (celltype in names(celltype_counts)) {
  indices <- which(cell_data$celltype == celltype)
  count <- length(indices)

  if (count >= 3) {
    sample_size <- max(1, min(floor(count * 0.5), count - 1))
  } else if (count == 2) {
    sample_size <- 1
  } else {
    sample_size <- ifelse(runif(1) > 0.5, 1, 0)
  }

  if (sample_size > 0) {
    sampled <- sample(indices, sample_size)
    split_groups[sampled] <- "group1"
  }
}

# 创建分组对象
obj@meta.data$split_group <- split_groups
obj_group1 <- subset(obj, subset = split_group == "group1")
obj_group2 <- subset(obj, subset = split_group == "group2")

# 保存细胞列表
write.table(Cells(obj_group1),
  file = file.path(output_dir, "LC_split_group1.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(Cells(obj_group2),
  file = file.path(output_dir, "LC_split_group2.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Pseudo-bulk分析
pseudo_bulk <- AggregateExpression(obj,
  group.by = "split_group",
  assays = "RNA", slot = "counts", return.seurat = FALSE
)

# 保存结果
write.csv(pseudo_bulk$RNA[, "group1", drop = FALSE],
  file = file.path(output_dir, "LC_pseudo_bulk_group1.csv")
)
write.csv(pseudo_bulk$RNA[, "group2", drop = FALSE],
  file = file.path(output_dir, "LC_pseudo_bulk_group2.csv")
)

# 简要输出
cat("分层抽样完成\n")
cat("组1细胞数:", ncol(obj_group1), "组2细胞数:", ncol(obj_group2), "\n")
cat("结果已保存到:", output_dir, "\n")
