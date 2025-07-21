library(ggplot2)
library(Seurat)
library(dplyr)

obj <- readRDS("./dataLib/tillering_node_rhizome/rhizome_S0.rds")

DimPlot(obj, label = TRUE) + NoAxes() + NoLegend()

obj[["celltype"]] <- obj$seurat_clusters
obj[["celltype"]] <- as.character(obj$celltype)

obj$celltype[obj$celltype %in% c("18", "14")] <- "Cell-Cycle"
obj$celltype[obj$celltype %in% c("13", "11")] <- "Xylem"
obj$celltype[obj$celltype %in% c("23")] <- "Tracheary-Element"
obj$celltype[obj$celltype %in% c("5", "20", "16", "21")] <- "Phloem & Companion-Cell"
obj$celltype[obj$celltype %in% c("19", "9")] <- "Cortex"
obj$celltype[obj$celltype %in% c("15", "17", "8")] <- "Endodermis & Exodermis & Pericycle"
obj$celltype[obj$celltype %in% c("6")] <- "Cambium"
obj$celltype[obj$celltype %in% c("4", "7", "22")] <- "Procambium & LRC(Primordia)"
obj$celltype[obj$celltype %in% c("1", "2")] <- "QC-like & Root-meristem"
obj$celltype[obj$celltype %in% c("0", "10")] <- "Shoot-meristem-like"
obj$celltype[obj$celltype %in% c("3")] <- "Meristematic-Vasculature"
obj$celltype[obj$celltype %in% c("12")] <- "Stele-like"


DimPlot(obj, group.by = "celltype")

saveRDS(obj, "./dataLib/tillering_node_rhizome/rhizome_S0.rds")

# Cellcycle: "Bochr12G354610", "Bochr04G146540", "Bochr06G191150", "Bochr06G199430"
cellcycle <- c("Bochr03G100050", "Bochr08G266560", "Bochr05G179090")
# Xylem: "Bochr10G309090","Bochr08G263440", "Bochr09G281850", "Bochr01G022590", "Bochr04G156770"
xylem <- c("Bochr10G305000", "Bochr09G282600")
protoxylem <- c("Bochr08G263530")
TrachearyElement <- c("Bochr01G048000", "Bochr05G161730")

# Phloem and Companion-Cell: "Bochr12G351210", "Bochr02G067580", "Bochr06G217330"
# Phloem and Companion-Cell: "Bochr01G004010", "Bochr02G053670", "Bochr06G213530"
phloem_cc <- c("Bochr06G210850", "Bochr06G214280", "Bochr01G047960")
phloem_pole_pericycle <- c("Bochr01G012110")

# Cortex :"Bochr05G184180", "Bochr03G094070", "Bochr01G043680"
cortex <- c("Bochr01G032380", "Bochr02G076060")

# Endodermis: "Bochr07G241900", "Bochr08G251390", "Bochr02G081430", "Bochr03G122170"
# Endodermis: "Bochr08G266350", "Bochr06G190220"
# Exodermis: "Bochr03G088140"
# Pericycle: "Bochr03G096190", "Bochr03G089080", "Bochr03G101610", "Bochr02G085600"
# Pericycle: "Bochr10G308710"
endodermis <- c("Bochr02G082030", "Bochr04G156440")
exodermis <- c("Bochr04G144330", "Bochr04G130210")
pericycle <- c("Bochr09G286170", "Bochr01G013200")

# Cambium: "Bochr04G147290", "Bochr03G113550", "Bochr04G157570", "Bochr04G136490"
cambium <- c("Bochr07G233300")

# Procambium and LRC(Primordia): "Bochr04G156440", "Bochr02G074070", "Bochr01G031640"
procambium <- c("Bochr03G094680")
lateral_root_primordia <- c("Bochr05G189490", "Bochr05G186270")

# QC-like and Root-meristem: "Bochr04G129990", "Bochr03G092290"
QC <- c("Bochr01G042780", "Bochr02G072640")
rc <- c("Bochr06G206520")

# Meristematic-Vasculature
meristematic_vasculature <- c("Bochr04G134390")

# Shoot-meristem-like: "Bochr03G115760",
smc <- c("Bochr03G119130", "Bochr03G119150", "Bochr07G220090")

# Stele-like: "Bochr05G187890" ,"Bochr05G177570", "Bochr02G070870", "Bochr03G099320", "Bochr02G083910"
stele_like <- c("Bochr02G080320", "Bochr12G343120", "Bochr08G266770")


cellcycle <- cellcycle
xylem <- c(xylem, protoxylem)
TrachearyElement <- TrachearyElement
phloem_cc <- c(phloem_cc, phloem_pole_pericycle)
cortex <- cortex
en_ex_pe <- c(endodermis, exodermis, pericycle)
cambium <- cambium
procambium_lrc <- c(procambium, lateral_root_primordia)
qc_like_root_mer <- c(QC, rc)
smc <- smc
meristematic_vasculature <- meristematic_vasculature
stele_like <- stele_like

names(cellcycle) <- rep("Cyc", length(cellcycle))
names(xylem) <- rep("Xylem", length(xylem))
names(TrachearyElement) <- rep("TE", length(TrachearyElement))
names(phloem_cc) <- rep("Phloem&CC", length(phloem_cc))
names(cortex) <- rep("Cortex", length(cortex))
names(en_ex_pe) <- rep("EnExPe", length(en_ex_pe))
names(cambium) <- rep("Cam", length(cambium))
names(procambium_lrc) <- rep("Proc&LRC", length(procambium_lrc))
names(qc_like_root_mer) <- rep("QC", length(qc_like_root_mer))
names(smc) <- rep("SMC-like", length(smc))
names(meristematic_vasculature) <- rep("MV", length(meristematic_vasculature))
names(stele_like) <- rep("Stele-like", length(stele_like))

# dotplot
ff <- c(
  cellcycle, xylem, TrachearyElement,
  phloem_cc, cortex, en_ex_pe, cambium,
  procambium_lrc, qc_like_root_mer, smc,
  meristematic_vasculature, stele_like
)

p <- DotPlot(obj, features = ff, cols = c("#eeeeee", "#ff4800"), col.min = 0, col.max = 3) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 80, vjust = 0.5, hjust = 0.5, size = 10)
  )

lev <- c(
  "14", "18", "13", "11", "23", "5", "20", "16", "21", "19",
  "9", "15", "17", "8", "6", "4", "7", "22", "1", "2", "0", "10", "3", "12"
)

p$data$id <- factor(p$data$id, levels = lev)

ggsave("./figures/tillering_node_rhizome/rhizome_S0_dot.pdf", p, width = 10, height = 5)

#
OsARF10 <- c("Bochr04G148850")
OsIAA24 <- c("Bochr07G223100")
#
IPA1 <- c("Bochr08G267750")
GRF6 <- c("Bochr03G119330")
HDT701 <- c("Bochr05G189490")
GL1 <- c("Bochr06G212770")
HDZIP13 <- c("Bochr03G113910")
GSTF6 <- c("Bochr10G309770")
WOX4 <- c("Bochr04G157570")
DSH1 <- c("Bochr06G197770")
PLT3 <- c("Bochr02G072640")
PLT2 <- c("Bochr06G213020")
PLT1 <- c("Bochr04G157850")
AP25 <- c("Bochr03G094070")
OsbHLH068 <- c("Bochr04G156440")

epi_ <- c("Bochr06G209350") #
