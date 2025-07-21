library(dplyr)
library(tidyr)
library(readxl)
library(ComplexHeatmap)

# read go data
# Function to read GO data
read_go_data <- function(file_path) {
  sheets <- excel_sheets(file_path)
  go_data <- lapply(sheets, function(sheet) {
    df <- read_excel(file_path, sheet = sheet)
    df$ct <- sheet
    df
  })
  names(go_data) <- sheets
  go_data
}

# Read up and down regulated GO data
up_go <- read_go_data("./dataLib/athLeaf/go/LC_LPZ_wound_distinct_up_go_mx.xlsx")
down_go <- read_go_data("./dataLib/athLeaf/go/LC_LPZ_wound_distinct_down_go_mx.xlsx")

up_go <- do.call(rbind, up_go)
down_go <- do.call(rbind, down_go)

# Prepare data for heatmap
# Create heatmap from down_go data
plot_go_heatmap <- function(go_data, output_file, height = 4.6, width = 6.6) {
  # Prepare data matrix
  dat <- go_data %>%
    dplyr::select(ID, Description, pvalue, ct) %>%
    pivot_wider(names_from = ct, values_from = pvalue) %>%
    as.data.frame()

  # Set row names and handle missing values
  rownames(dat) <- dat$ID
  dat[is.na(dat)] <- 1

  # Define color scale
  min_p <- min(go_data$pvalue)
  col_fun <- circlize::colorRamp2(
    c(min_p, 0.05, 1),
    c("#005326", "#4cae78", "#f6fffa")
  )

  # Select data columns for heatmap
  data_cols <- setdiff(colnames(dat), c("ID", "Description"))

  # Create heatmap
  hp <- Heatmap(
    as.matrix(dat[, data_cols]),
    col = col_fun,
    name = "P-Value",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_side = "left",
    column_names_rot = 70,
    rect_gp = gpar(col = "grey"),
    heatmap_height = unit(6, "mm") * nrow(dat),
    heatmap_width = unit(13, "mm") * length(data_cols),
    heatmap_legend_param = list(
      at = c(min_p, 0.05, 1),
      legend_width = unit(6, "cm"),
      direction = "horizontal",
      title_position = "lefttop"
    )
  )

  # Add row annotations
  ra <- rowAnnotation(Symbol = anno_text(dat$Description,
    location = 0.01,
    just = "left"
  ))

  # Save plot
  pdf(output_file, height = height, width = width)
  draw(hp + ra, heatmap_legend_side = "top")
  dev.off()
}

# Plot heatmap
plot_go_heatmap(up_go, "./figures/leaf/LCLPZ_distinct_go_up_heatmap.pdf", height = 6, width = 6.6)
plot_go_heatmap(down_go, "./figures/leaf/LCLPZ_distinct_go_down_heatmap.pdf", height = 5, width = 7)


# add 42 genes go 
go_42 <- read.csv("./dataLib/athLeaf/gores_42_genes.csv")
head(go_42[,1:4], 10)

gg <- go_42[c(1:5, 8, 9), ] %>% mutate(ct = "share") %>% bind_rows(down_go)

plot_go_heatmap(gg, "./figures/leaf/LCLPZ_distinct_go_down_heatmap_42.pdf", height = 6, width = 8)

