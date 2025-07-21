library(Seurat)
library(DESeq2)
library(dplyr)

tq <- readRDS("./dataLib/osRoot/ggm/osRoots_TQ_v5.rds")
tq <- AddMetaData(tq, metadata = tq$orig.ident, col.name = "samples")
tq$samples <- ifelse(tq$samples == "10X_osRoot1_singlet", "TQ1", "TQ2")
table(tq$samples)

fx <- readRDS("./dataLib/osRoot/ggm/ggm12_normal.Rds")
fx <- AddMetaData(fx, metadata = fx$orig.ident, col.name = "samples")
fx$samples <- ifelse(fx$samples == "ggm1", "FX1", "FX2")
table(fx$samples)

# 使用AddMetaData添加组别信息
tq <- AddMetaData(tq, metadata = "TQ", col.name = "group")
fx <- AddMetaData(fx, metadata = "FX", col.name = "group")

# 合并数据集
combined <- merge(tq, fx, add.cell.ids = c("TQ", "FX"))

# 使用AddMetaData创建样本和组别的组合标识
combined <- AddMetaData(combined,
  metadata = paste(combined$samples, combined$group, sep = "_"),
  col.name = "sample_group"
)

# 设置细胞标识为样本名
Idents(combined) <- "samples"

print("正在使用Seurat进行pseudo-bulk聚合...")

# 使用Seurat的AggregateExpression函数进行pseudo-bulk聚合
pseudo_bulk <- AggregateExpression(combined,
  group.by = "samples",
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

# 提取计数矩阵
pseudo_bulk_counts <- pseudo_bulk$RNA

print(paste("聚合后样本数量:", ncol(pseudo_bulk_counts)))
print(paste("基因数量:", nrow(pseudo_bulk_counts)))

# 创建样本元数据，使用Seurat函数提取信息
sample_meta <- combined@meta.data %>%
  select(samples, group) %>%
  distinct() %>%
  arrange(samples)

# 确保样本顺序一致
sample_meta <- sample_meta[match(colnames(pseudo_bulk_counts), sample_meta$samples), ]

# 过滤低表达基因（在至少一半样本中表达量>10）
keep_genes <- rowSums(pseudo_bulk_counts > 10) >= ncol(pseudo_bulk_counts) / 2
pseudo_bulk_counts <- pseudo_bulk_counts[keep_genes, ]

print(paste("保留", nrow(pseudo_bulk_counts), "个基因进行分析"))

# 准备DESeq2分析
coldata <- sample_meta
rownames(coldata) <- coldata$samples
coldata$group <- factor(coldata$group, levels = c("FX", "TQ")) # FX作为对照组

# 确保计数矩阵和元数据样本顺序一致
pseudo_bulk_counts <- pseudo_bulk_counts[, rownames(coldata)]

# 创建DESeq2对象
dds <- DESeqDataSetFromMatrix(
  countData = pseudo_bulk_counts,
  colData = coldata,
  design = ~group
)

# 运行差异分析
print("正在运行DESeq2差异分析...")
dds <- DESeq(dds)

# 提取结果
res <- results(dds, contrast = c("group", "TQ", "FX")) # TQ vs FX
res_df <- as.data.frame(res)

# 添加基因名称
res_df$gene <- rownames(res_df)

# 移除包含NA的行
res_df <- res_df[complete.cases(res_df), ]


# normalized expression
print("正在提取标准化表达量...")

# 从DESeq2对象提取标准化表达量
normalized_counts <- counts(dds, normalized = TRUE)

# 转换为数据框并添加基因名
normalized_df <- as.data.frame(normalized_counts)
normalized_df$gene <- rownames(normalized_df)

# 为标准化表达量列添加前缀，便于识别
colnames(normalized_df)[1:(ncol(normalized_df) - 1)] <- paste0(
  "norm_", colnames(normalized_df)[1:(ncol(normalized_df) - 1)]
)

# 将差异分析结果与标准化表达量合并
print("正在合并差异分析结果和标准化表达量...")
merged_results <- merge(res_df, normalized_df, by = "gene", all.x = TRUE)


# 保存结果
merged_file <- "./dataLib/osRoot/pseudobulk_TQ_vs_FX_proto.csv"

write.csv(merged_results, file = merged_file, row.names = FALSE)
