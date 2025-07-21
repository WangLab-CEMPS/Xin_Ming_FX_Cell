library(ggplot2)
library(Seurat)
library(dplyr)

obj <- readRDS("./dataLib/osRoot/ggm/ggm2_TQ2_strata.Rds")

mm <- read.csv("./dataLib/osRoot/root_tq_strata_sampling_51.csv")
tq <- readRDS("./dataLib/osRoot/osRoots_TQ_v5.rds")
tq <- subset(tq, cells = mm$cell)
tq_meta <- tq@meta.data
tq_meta <- tibble(
  nn = tq_meta$seurat_clusters,
  cells = rownames(tq_meta)
)

mm <- obj@meta.data %>%
  as_tibble(rownames = "cells") %>%
  left_join(tq_meta) %>%
  as.data.frame()

rownames(mm) <- mm$cells

obj@meta.data <- mm

colnames(obj@meta.data)

p <- DimPlot(obj, reduction = "harmony_umap", group.by = "nn", label = TRUE) &
  tidydr::theme_dr() &
  theme(panel.grid = element_blank(), plot.title = element_blank())

ggsave("./Plots/osroot/ggm2_TQ2_harmony_umap_strata_labeledbyTQSeurat.pdf", p, width = 5, height = 4)


ep <- c("Os10g0122600", "Os03g0150800", "Os02g0595900", "Os12g0637100")
rh <- c("Os10g0546100", "Os07g0499500", "Os01g0164300", "Os10g0453900", "Os10g0454200")
rc <- c("Os04g0552000", "Os01g0248000")
ex <- c("Os04g0452700", "Os03g0107300")
ms <- c(
  "Os05g0438700", "Os08g0490900", "Os03g0279200", "Os02g0805200",
  "Os02g0829100", "Os08g0512600", "Os01g0805600", "Os12g0555600",
  "Os04g0486500", "Os02g0810200", "Os06g0127800", "Os07g0105700"
)
co <- c(
  "Os01g0745700", "Os05g0597100", "Os03g0729500", "Os01g0131600",
  "Os05g0520300", "Os02g0662000", "Os04g0554600"
)

vc <- c("Os03g0820500", "Os04g0540900", "Os05g0420300")
x <- c("Os10g0467800", "Os09g0422500", "Os04g0536500", "Os03g0640800")
p <- c("Os06g0614000", "Os06g0676000", "Os01g0236300")
en <- c(
  "Os07g0134000", "Os02g0745100", "Os04g0684300", "Os12g0122000",
  "Os01g0679700", "Os05g0477300"
)
pe <- c("Os03g0216700", "Os10g0524300", "Os07g0525100", "Os04g0445000")

names(ep) <- rep("ep", length(ep))
names(rh) <- rep("rh", length(rh))
names(rc) <- rep("rc", length(rc))
names(ex) <- rep("ex", length(ex))
names(ms) <- rep("ms", length(ms))
names(vc) <- rep("vc", length(vc))
names(x) <- rep("Xylem", length(x))
names(p) <- rep("Phloem", length(p))
names(en) <- rep("en", length(en))
names(pe) <- rep("pe", length(pe))
names(co) <- rep("co", length(co))

ff <- c(ep, rh, rc, ex, vc, x, p, pe, en, ms, co)


p2 <- DotPlot(obj, features = ff, cols = c("#eeeeee", "#ff4800"), col.min = 0, col.max = 3) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 80, vjust = 0.5, hjust = 0.5, size = 10)
  )

lev <- c(
  "5", "18", "20", "14", "7",
  "11", "22", "10",
  "15", "16", "19", "6",
  "4", "1", "13",
  "21", "0", "2", "17",
  "9", "3", "12", "8"
)

p2$data$id <- factor(p2$data$id, levels = lev)
p2

ggsave("./Plots/osroot/ggm2_TQ2_harmony_dot_strata.pdf", p2, width = 12, height = 5)


## ------------
obj[["celltype"]] <- obj$seurat_clusters %>% as.character()

obj$celltype[obj$celltype %in% c("5", "18")] <- "epidermis/root hair"
obj$celltype[obj$celltype %in% c("7", "8", "14")] <- "exodermis"
obj$celltype[obj$celltype %in% c("12")] <- "sclerenchyma"
obj$celltype[obj$celltype %in% c("3", "9")] <- "cortex"
obj$celltype[obj$celltype %in% c("0", "13", "21")] <- "endodermis"
obj$celltype[obj$celltype %in% c("1")] <- "pericycle"
obj$celltype[obj$celltype %in% c("2", "17")] <- "meristem"
obj$celltype[obj$celltype %in% c("10", "22")] <- "vascular cylinder"
obj$celltype[obj$celltype %in% c("15", "16", "4")] <- "xylem"
obj$celltype[obj$celltype %in% c("6", "19")] <- "phloem"
obj$celltype[obj$celltype %in% c("11")] <- "putative root cap junction"
obj$celltype[obj$celltype %in% c("20")] <- "root cap"



p3 <- DimPlot(
  obj,
  reduction = "harmony_umap",
  group.by = "celltype",
  cols = pigRutils::select_colors("col_33")
) & tidydr::theme_dr() &
  theme(panel.grid = element_blank(), plot.title = element_blank())

ggsave("./Plots/osroot/ggm2_TQ2_harmony_strata_Anno.pdf", p3, width = 8, height = 6)


saveRDS(obj, "./dataLib/osRoot/ggm2_TQ2_strata.Rds")
