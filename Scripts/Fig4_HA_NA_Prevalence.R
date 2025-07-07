#================================================
#### Packages
#================================================
library(openxlsx)
library(circlize)
library(ComplexHeatmap)
library(grid)

#================================================
#### HA/NA Prevalence Heatmap
#================================================
meta <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
data <- read.xlsx("HA_NA_Tree_Metadata.xlsx", rowNames = TRUE)

missing_rows <- setdiff(rownames(meta), rownames(data))
zero_rows <- matrix(0,
                    nrow = length(missing_rows),
                    ncol = ncol(data),
                    dimnames = list(missing_rows, colnames(data)))
zero_df <- as.data.frame(zero_rows)
data_full <- rbind(data, zero_df)
data_full <- data_full[rownames(meta), ]
meta_cols <- meta[, c("Location_Grouping", "Experiment", "Location")]
meta_cols <- meta_cols[rownames(data_full), ]
data_annotated <- cbind(data_full, meta_cols)

group_palette <- c(
  "Chicken_Throat" = "#009E73",
  "Chicken_Cloacal" = "#800000",       
  "Duck_Throat" = "#0072B2",
  "Duck_Cloacal" = "#964B00",          
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",          
  "Air_Outside" = "#228B22",           
  "Wash_Water" = "#F0E442",
  "Drinking_Water" = "#56B4E9",
  "Cage" = "#D55E00")

Heatmap_Mat <- data_annotated[, !(colnames(data_annotated) %in% c("Location", "Experiment", "Location_Grouping"))]
Heatmap_Mat <- as.matrix(Heatmap_Mat)
Split <- data_annotated$Experiment
Split <- factor(Split)
col_fuPA <- colorRamp2(c(0, 1), c("white", "red"))
Annotation_Meta <- data_annotated[, c("Location", "Location_Grouping")]

Colors_Top <- list(
  "Location" = c("Phnom_Penh" = "grey", "Takeo" = "black"),
  "Location_Grouping" = group_palette  
)

Row_Annotations <- rowAnnotation(
  df = Annotation_Meta,
  col = Colors_Top,
  show_annotation_name = FALSE
)

HA_NA <- Heatmap(
  Heatmap_Mat,
  row_split = Split,
  col = col_fuPA,
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  show_column_dend = FALSE,
  show_row_names = FALSE,  
  show_heatmap_legend = FALSE,
  column_names_rot = 45,
  rect_gp = gpar(col = "black", lwd = 1.5),
  column_dend_height = unit(2, "cm"),
  column_dend_gp = gpar(fontface = "bold", fontsize = 12),
  column_gap = unit(10, "mm"),
  row_names_gp = gpar(fontface = "bold", fontsize = 10),
  column_names_gp = gpar(fontface = "bold"),
  left_annotation = Row_Annotations,
)
print(HA_NA)

pdf("HA_NA_Prevalence.pdf", width = 4, height = 12) 
draw(HA_NA)
dev.off()
