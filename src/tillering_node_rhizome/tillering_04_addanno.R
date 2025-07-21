library(ggplot2)
library(Seurat)
library(dplyr)

obj <- readRDS("./dataLib/tillering_node_rhizome/tillering.rds")
DimPlot(obj, label = TRUE)

obj[["celltype"]] <- obj$seurat_clusters
obj[["celltype"]] <- as.character(obj$celltype)

obj$celltype[obj$celltype %in% c("2")] <- "Cell-Cycle"
obj$celltype[obj$celltype %in% c("0")] <- "Meristem-like"
obj$celltype[obj$celltype %in% c("10")] <- "Xylem"
obj$celltype[obj$celltype %in% c("7")] <- "Xylem & Cambium"
obj$celltype[obj$celltype %in% c("13")] <- "Cambium"
obj$celltype[obj$celltype %in% c("8", "9")] <- "Phloem"
obj$celltype[obj$celltype %in% c("4")] <- "Pericycle & Parenchyma"
obj$celltype[obj$celltype %in% c("5")] <- "Epidermis"
obj$celltype[obj$celltype %in% c("1", "3")] <- "Endodermis"
obj$celltype[obj$celltype %in% c("12")] <- "Exodermis"
obj$celltype[obj$celltype %in% c("11")] <- "Root-cap"
obj$celltype[obj$celltype %in% c("6")] <- "Unkown"


DimPlot(obj, group.by = "celltype")

saveRDS(obj, "./dataLib/tillering_node_rhizome/tillering.rds")

# cellcycle: "Os12g0555600"
cyc <- c("Os03g0279200", "Os08g0490900", "Os05g0438700")
# ep: "Os01g0822900", "Os01g0201600", "Os03g0181500", "Os03g0729500"
ep <- c("Os02g0776900", "Os02g0813600", "Os06g0653000")
vc <- c("Os04g0540900", "Os05g0420300", "Os03g0329700", "Os06g0216700")
# xylem: "Os03g0640800",
xylem <- c("Os10g0467800", "Os09g0422500", "Os10g0532000", "Os02g0237100")
# phloem: "Os06g0614000", "Os06g0676000"
phloem <- c("Os06g0717200", "Os01g0971000", "Os10g0543800", "Os04g0105200")
# cortex: "Os01g0745700", "Os05g0597100"
co <- c("Os02g0662000", "Os01g0914300")
en <- c("Os02g0745100", "Os12g0122000")
ex <- c("Os04g0452700", "Os03g0107300")
pe <- c("Os03g0216700", "Os10g0524300", "Os02g0809800", "Os04g0556000")
rc <- c("Os04g0552000", "Os01g0248000")
cam <- c("Os02g0614300", "Os04g0649400")
mer <- c("Os04g0653600", "Os03g0764900")
osh1 <- c("Os03g0727000")
osh15 <- c("Os07g0129700")
unk <- c("Os10g0532300", "Os01g0232000")

names(cyc) <- rep("Cyc", length(cyc))
names(ep) <- rep("Epi", length(ep))
names(xylem) <- rep("Xylem", length(xylem))
names(phloem) <- rep("Phloem", length(phloem))
names(en) <- rep("En", length(en))
names(ex) <- rep("Ex", length(ex))
names(pe) <- rep("Pe", length(pe))
names(rc) <- rep("RC", length(rc))
names(cam) <- rep("Cam", length(cam))
names(mer) <- rep("Mer", length(mer))
names(unk) <- rep("Unk", length(unk))

# dotplot
ff <- c(cyc, mer, xylem, cam, phloem, pe, ep, en, ex, rc, unk)
p <- DotPlot(obj, features = ff, cols = c("#eeeeee", "#ff4800"), col.min = 0, col.max = 3) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 80, vjust = 0.5, hjust = 0.5, size = 10)
  )

lev <- c(
  "2", "0", "10", "7", "13", "8", "9",
  "4", "5", "1", "3", "12", "11", "6"
)

p$data$id <- factor(p$data$id, levels = lev)
p

ggsave("./figures/tillering_node_rhizome/tillering_dot.pdf", p, width = 10, height = 5)
