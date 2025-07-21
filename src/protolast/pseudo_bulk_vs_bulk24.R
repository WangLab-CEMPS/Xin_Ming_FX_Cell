library(DESeq2)
library(dplyr)
library(pigRutils)

# read in data ----
bulk_raw <- read.table("./dataLib/athLeaf/bulk_proto/raw_counts.txt", sep = "\t")
head(bulk_raw)

# read in pseudo bulk
g1 <- read.csv("./dataLib/proto/LC_pseudo_bulk_group1.csv")
head(g1)
g2 <- read.csv("./dataLib/proto/LC_pseudo_bulk_group2.csv")
head(g2)

pseudo <- inner_join(g1, g2, by = "X")
head(pseudo)

# merge data
cts <- bulk_raw %>%
  mutate(X = rownames(.)) %>%
  inner_join(pseudo, by = "X")
rownames(cts) <- cts$X
cts <- cts %>% select(-X)
head(cts)


# run DESeq2
(
  coldata <- data.frame(
    samples = colnames(cts),
    group = c(gsub("\\.[1-4]", "", colnames(cts)[1:15]), "pseudo", "pseudo")
  )
)


dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~group
)

dds <- DESeq(dds)

# compare ionfo

com <- CompInfo(colData = coldata, condition = "group")

apply(
  com, 1, one_one,
  Type = "group",
  cutoff_padj = 0.05,
  od = "./dataLib/proto/new",
  seq = "RNA"
)
