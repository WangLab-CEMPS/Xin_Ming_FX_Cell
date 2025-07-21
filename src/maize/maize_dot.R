library(Seurat)
library(dplyr)
library(ggplot2)

obj <- readRDS("./dataLib/maize/maize.rds")




# endo Zm00001eb083140 Zm00001eb422410 Zm00001eb249210
# cyc Zm00001eb370670 Zm00001eb035350 Zm00001eb188180
# phloem Zm00001eb266840 Zm00001eb315650 Zm00001eb002390 Zm00001eb342070
# Phloem & Sieve  Zm00001eb266840 Zm00001eb002390 Zm00001eb173060 Zm00001eb315650 Zm00001eb342070
# Pericycle       Zm00001eb356620
# Pericycle       Zm00001eb397500 Zm00001eb302250
# Stele           Zm00001eb175860
# Xylem           Zm00001eb267050
## Stele 7        Zm00001eb338510
# Cortex          Zm00001eb000450 Zm00001eb330530 Zm00001eb249760 Zm00001eb429560
# Cortex Initials Zm00001eb069630
# CortexII(M)     Zm00001eb187430
# Endodermal Initials Zm00001eb388430
## Epidermis      Zm00001eb255380, Zm00001eb190830
# QC              Zm00001eb200450
# Root cap        Zm00001eb091140 Zm00001eb140000 Zm00001eb091170 Zm00001eb118690 Zm00001eb196940
features <- c(
  "Zm00001eb372230", "Zm00001eb422410", "Zm00001eb249210",
  "Zm00001eb370670", "Zm00001eb035350", "Zm00001eb188180",
  "Zm00001eb397500", "Zm00001eb302250",
  "Zm00001eb175860", "Zm00001eb338510",
  "Zm00001eb267050",
  "Zm00001eb249760", "Zm00001eb429560",
  "Zm00001eb069630",
  "Zm00001eb187430",
  "Zm00001eb091140", "Zm00001eb091170", "Zm00001eb140000",
  "Zm00001eb200450"
)


p <- DotPlot(obj, features = features, cols = c("#eeeeee", "#ff4800"), col.min = 0, col.max = 3) +
  theme_bw() +
  coord_flip() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text = element_text(face = "bold")
  )

lev <- c(6, 18, 11, 13, 14, 15, 2, 8, 0, 10, 7, 9, 1, 5, 4, 16, 3, 17, 19, 12)

p$data$id <- factor(p$data$id, levels = lev)
p

ggsave("./figures/maize/maize_dot.pdf", p, width = 7, height = 6)
