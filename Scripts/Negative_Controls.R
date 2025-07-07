#================================================
#### Packages
#================================================
library(openxlsx)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(grid)

Species_RPM <- read.xlsx("Species_RPM_Negative_Controls.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Negative_Controls_Metadata.xlsx", rowNames = TRUE)

rank_scale <- function(x)
{
  x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  return(y);
}

Species_Rank <- apply(Species_RPM,2,rank_scale)
Species_Rank <- as.data.frame(Species_Rank)
Species_Mat <- as.matrix(t(Species_Rank))
Clustering <- hclust(dist(Species_Mat, method = "euclidean"), method = "ward.D")

Chicken_Throat_Metadata <- subset(Metadata, Location_Grouping=="Chicken_Throat")
Chicken_Cloacal_Metadata <- subset(Metadata, Location_Grouping=="Chicken_Cloacal")
Duck_Throat_Metadata <- subset(Metadata, Location_Grouping=="Duck_Throat")
Duck_Cloacal_Metadata <- subset(Metadata, Location_Grouping=="Duck_Cloacal")
Air_Slaughter_Metadata <- subset(Metadata, Location_Grouping=="Air_Slaughter_Area")
Air_Holding_Metadata <- subset(Metadata, Location_Grouping=="Air_Holding_Area")
Air_Outside_Metadata <- subset(Metadata, Location_Grouping=="Air_Outside")
Wash_Water_Metadata <- subset(Metadata, Location_Grouping=="Environmental_Wash_Water")
Drinking_Water_Metadata <- subset(Metadata, Location_Grouping=="Environmental_Drinking_Water")
Cage_Metadata <- subset(Metadata, Location_Grouping=="Environmental_Cage")
Neg_Metadata <- subset(Metadata, Location_Grouping=="Negative_Control")

Chicken_Throat_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Chicken_Throat_Metadata),] 
Chicken_Cloacal_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Chicken_Cloacal_Metadata),] 
Duck_Throat_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Duck_Throat_Metadata),]
Duck_Cloacal_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Duck_Cloacal_Metadata),]
Air_Slaughter_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Air_Slaughter_Metadata),]
Air_Holding_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Air_Holding_Metadata),]
Air_Outside_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Air_Outside_Metadata),]
Wash_Water_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Wash_Water_Metadata),]
Drinking_Water_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Drinking_Water_Metadata),]
Cage_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Cage_Metadata),]
Neg_Species <- Species_Rank[row.names(Species_Rank) %in% row.names(Neg_Metadata),]

Chicken_Throat_Mat <- as.matrix(Chicken_Throat_Species)
Chicken_Cloacal_Mat <- as.matrix(Chicken_Cloacal_Species)
Duck_Throat_Mat <- as.matrix(Duck_Throat_Species)
Duck_Cloacal_Mat <- as.matrix(Duck_Cloacal_Species)
Air_Slaughter_Mat <- as.matrix(Air_Slaughter_Species)
Air_Holding_Mat <- as.matrix(Air_Holding_Species)
Air_Outside_Mat <- as.matrix(Air_Outside_Species)
Wash_Water_Mat <- as.matrix(Wash_Water_Species)
Drinking_Water_Mat <- as.matrix(Drinking_Water_Species)
Cage_Mat <- as.matrix(Cage_Species)
Neg_Mat <- as.matrix(Neg_Species)

Chicken_Throat_Clust <- hclust(dist(Chicken_Throat_Mat, method = "euclidean"), method = "ward.D")
Chicken_Cloacal_Clust <- hclust(dist(Chicken_Cloacal_Mat, method = "euclidean"), method = "ward.D")
Duck_Throat_Clust <- hclust(dist(Duck_Throat_Mat, method = "euclidean"), method = "ward.D")
Duck_Cloacal_Clust <- hclust(dist(Duck_Cloacal_Mat, method = "euclidean"), method = "ward.D")
Air_Slaughter_Clust <- hclust(dist(Air_Slaughter_Mat, method = "euclidean"), method = "ward.D")
Air_Holding_Clust <- hclust(dist(Air_Holding_Mat, method = "euclidean"), method = "ward.D")
Air_Outside_Clust <- hclust(dist(Air_Outside_Mat, method = "euclidean"), method = "ward.D")
Wash_Water_Clust <- hclust(dist(Wash_Water_Mat, method = "euclidean"), method = "ward.D")
Drinking_Water_Clust <- hclust(dist(Drinking_Water_Mat, method = "euclidean"), method = "ward.D")
Cage_Clust <- hclust(dist(Cage_Mat, method = "euclidean"), method = "ward.D")
Neg_Clust <- hclust(dist(Neg_Mat, method = "euclidean"), method = "ward.D")

Clustering_List <- as.data.frame(Clustering$labels[Clustering$order])
#Color_Scheme = brewer.pal(8,"RdBu")
Color_Scheme = brewer.pal(9,"Reds")

Neg_Meta <- Neg_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Negative_Control" = "black"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Throat" = "yellow", "Cloacal" = "pink"))
Neg_Top_Annotations = HeatmapAnnotation(df = Neg_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Neg_Sample_Type <- Neg_Metadata[c(5)]
#Neg_Bottom_Annotations = HeatmapAnnotation(df = Neg_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Neg_Heatmap <- Heatmap(t(Neg_Mat),
                              name = "Rank(RPM)",
                              col = Color_Scheme,
                              rect_gp = gpar(col = "white", lwd = 0.1),
                              cluster_rows = Clustering,
                              cluster_columns = Neg_Clust,
                              show_row_names = FALSE,
                              show_column_names = FALSE,
                              show_column_dend = FALSE,
                              show_row_dend = FALSE,
                              #column_title = "Gallus Gallus",
                              column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                              column_title_side = "bottom",
                              column_dend_height = unit(2, "cm"),
                              row_dend_width = unit(10, "cm"),
                              show_heatmap_legend = FALSE,
                              top_annotation = Neg_Top_Annotations)
#bottom_annotation = Neg_Bottom_Annotations)
#heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 1)))

Neg_Heatmap

Chicken_Throat_Meta <- Chicken_Throat_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Chicken_Throat" = "#009E73"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Throat" = "yellow", "Cloacal" = "pink"))
Chicken_Throat_Top_Annotations = HeatmapAnnotation(df = Chicken_Throat_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Chicken_Sample_Type <- Chicken_Metadata[c(5)]
#Chicken_Bottom_Annotations = HeatmapAnnotation(df = Chicken_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Chicken_Throat_Heatmap <- Heatmap(t(Chicken_Throat_Mat),
                                  name = "Rank(RPM)",
                                  col = Color_Scheme,
                                  rect_gp = gpar(col = "white", lwd = 0.1),
                                  cluster_rows = Clustering,
                                  cluster_columns = Chicken_Throat_Clust,
                                  show_row_names = FALSE,
                                  show_column_names = FALSE,
                                  show_column_dend = FALSE,
                                  show_row_dend = FALSE,
                                  #column_title = "Gallus Gallus",
                                  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                                  column_title_side = "bottom",
                                  column_dend_height = unit(2, "cm"),
                                  row_dend_width = unit(10, "cm"),
                                  show_heatmap_legend = FALSE,
                                  top_annotation = Chicken_Throat_Top_Annotations)
#bottom_annotation = Chicken_Bottom_Annotations)
#heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 1)))

Chicken_Throat_Heatmap

Chicken_Cloacal_Meta <- Chicken_Cloacal_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Chicken_Cloacal" = "#800000"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Cloacal" = "yellow", "Cloacal" = "pink"))
Chicken_Cloacal_Top_Annotations = HeatmapAnnotation(df = Chicken_Cloacal_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Chicken_Sample_Type <- Chicken_Metadata[c(5)]
#Chicken_Bottom_Annotations = HeatmapAnnotation(df = Chicken_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Chicken_Cloacal_Heatmap <- Heatmap(t(Chicken_Cloacal_Mat),
                                   name = "Rank(RPM)",
                                   col = Color_Scheme,
                                   rect_gp = gpar(col = "white", lwd = 0.1),
                                   cluster_rows = Clustering,
                                   cluster_columns = Chicken_Cloacal_Clust,
                                   show_row_names = FALSE,
                                   show_column_names = FALSE,
                                   show_column_dend = FALSE,
                                   show_row_dend = FALSE,
                                   #column_title = "Gallus Gallus",
                                   column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                                   column_title_side = "bottom",
                                   column_dend_height = unit(2, "cm"),
                                   row_dend_width = unit(10, "cm"),
                                   show_heatmap_legend = FALSE,
                                   top_annotation = Chicken_Cloacal_Top_Annotations)
#bottom_annotation = Chicken_Bottom_Annotations)
#heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 1)))

Chicken_Cloacal_Heatmap

Duck_Throat_Meta <- Duck_Throat_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Duck_Throat" = "#0072B2"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Throat" = "yellow", "Cloacal" = "pink"))
Duck_Throat_Top_Annotations = HeatmapAnnotation(df = Duck_Throat_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Duck_Throat_Sample_Type <- Duck_Throat_Metadata[c(5)]
#Duck_Throat_Bottom_Annotations = HeatmapAnnotation(df = Duck_Throat_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Duck_Throat_Heatmap <- Heatmap(t(Duck_Throat_Mat),
                               name = "Rank(RPM)",
                               col = Color_Scheme,
                               rect_gp = gpar(col = "white", lwd = 0.1),
                               cluster_rows = Clustering,
                               cluster_columns = Duck_Throat_Clust,
                               show_row_names = FALSE,
                               show_column_names = FALSE,
                               show_column_dend = FALSE,
                               show_row_dend = FALSE,
                               #column_title = "Anas Platyrhynchos",
                               column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                               column_title_side = "bottom",
                               column_dend_height = unit(2, "cm"),
                               show_heatmap_legend = FALSE,
                               top_annotation = Duck_Throat_Top_Annotations)
#bottom_annotation = Duck_Throat_Bottom_Annotations)
#heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 3))

Duck_Throat_Heatmap

Duck_Cloacal_Meta <- Duck_Cloacal_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Duck_Cloacal" = "#964B00"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Cloacal" = "yellow", "Cloacal" = "pink"))
Duck_Cloacal_Top_Annotations = HeatmapAnnotation(df = Duck_Cloacal_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Duck_Cloacal_Sample_Type <- Duck_Cloacal_Metadata[c(5)]
#Duck_Cloacal_Bottom_Annotations = HeatmapAnnotation(df = Duck_Cloacal_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Duck_Cloacal_Heatmap <- Heatmap(t(Duck_Cloacal_Mat),
                                name = "Rank(RPM)",
                                col = Color_Scheme,
                                rect_gp = gpar(col = "white", lwd = 0.1),
                                cluster_rows = Clustering,
                                cluster_columns = Duck_Cloacal_Clust,
                                show_row_names = FALSE,
                                show_column_names = FALSE,
                                show_column_dend = FALSE,
                                show_row_dend = FALSE,
                                #column_title = "Anas Platyrhynchos",
                                column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                                column_title_side = "bottom",
                                column_dend_height = unit(2, "cm"),
                                show_heatmap_legend = FALSE,
                                top_annotation = Duck_Cloacal_Top_Annotations)
#bottom_annotation = Duck_Cloacal_Bottom_Annotations)
#heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 3))

Duck_Cloacal_Heatmap

Air_Slaughter_Meta <- Air_Slaughter_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Air_Slaughter_Area" = "#CC79A7"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Throat" = "yellow", "Cloacal" = "pink"))
Air_Slaughter_Top_Annotations = HeatmapAnnotation(df = Air_Slaughter_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Air_Slaughter_Sample_Type <- Air_Slaughter_Metadata[c(5)]
#Air_Slaughter_Bottom_Annotations = HeatmapAnnotation(df = Air_Slaughter_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Air_Slaughter_Heatmap <- Heatmap(t(Air_Slaughter_Mat),
                                 name = "Rank(RPM)",
                                 col = Color_Scheme,
                                 rect_gp = gpar(col = "white", lwd = 0.1),
                                 cluster_rows = Clustering,
                                 cluster_columns = Air_Slaughter_Clust,
                                 show_row_names = FALSE,
                                 show_column_names = FALSE,
                                 show_column_dend = FALSE,
                                 show_row_dend = FALSE,
                                 #column_title = "Air_Slaughter",
                                 column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                                 column_title_side = "bottom",
                                 column_dend_height = unit(2, "cm"),
                                 show_heatmap_legend = FALSE,
                                 top_annotation = Air_Slaughter_Top_Annotations)
#heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 3)))

Air_Slaughter_Heatmap

Air_Holding_Meta <- Air_Holding_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Air_Holding_Area" = "#6A5ACD"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Throat" = "yellow", "Cloacal" = "pink"))
Air_Holding_Top_Annotations = HeatmapAnnotation(df = Air_Holding_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Air_Holding_Sample_Type <- Air_Holding_Metadata[c(5)]
#Air_Holding_Bottom_Annotations = HeatmapAnnotation(df = Air_Holding_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Air_Holding_Heatmap <- Heatmap(t(Air_Holding_Mat),
                               name = "Rank(RPM)",
                               col = Color_Scheme,
                               rect_gp = gpar(col = "white", lwd = 0.1),
                               cluster_rows = Clustering,
                               cluster_columns = Air_Holding_Clust,
                               show_row_names = FALSE,
                               show_column_names = FALSE,
                               show_column_dend = FALSE,
                               show_row_dend = FALSE,
                               #column_title = "Air_Holding",
                               column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                               column_title_side = "bottom",
                               column_dend_height = unit(2, "cm"),
                               show_heatmap_legend = FALSE,
                               top_annotation = Air_Holding_Top_Annotations)
#heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 3)))

Air_Holding_Heatmap

Air_Outside_Meta <- Air_Outside_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Air_Outside" = "#228B22"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Throat" = "yellow", "Cloacal" = "pink"))
Air_Outside_Top_Annotations = HeatmapAnnotation(df = Air_Outside_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Air_Outside_Sample_Type <- Air_Outside_Metadata[c(5)]
#Air_Outside_Bottom_Annotations = HeatmapAnnotation(df = Air_Outside_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Air_Outside_Heatmap <- Heatmap(t(Air_Outside_Mat),
                               name = "Rank(RPM)",
                               col = Color_Scheme,
                               rect_gp = gpar(col = "white", lwd = 0.1),
                               cluster_rows = Clustering,
                               cluster_columns = Air_Outside_Clust,
                               show_row_names = FALSE,
                               show_column_names = FALSE,
                               show_column_dend = FALSE,
                               show_row_dend = FALSE,
                               #column_title = "Air_Outside",
                               column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                               column_title_side = "bottom",
                               column_dend_height = unit(2, "cm"),
                               show_heatmap_legend = FALSE,
                               top_annotation = Air_Outside_Top_Annotations)
#heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 3)))

Air_Outside_Heatmap

Wash_Water_Meta <- Wash_Water_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Environmental_Wash_Water" = "#F0E442"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Throat" = "yellow", "Cloacal" = "pink"))
Wash_Water_Top_Annotations = HeatmapAnnotation(df = Wash_Water_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Wash_Water_Sample_Type <- Wash_Water_Metadata[c(5)]
#Wash_Water_Bottom_Annotations = HeatmapAnnotation(df = Wash_Water_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Wash_Water_Heatmap <- Heatmap(t(Wash_Water_Mat),
                              name = "Rank(RPM)",
                              col = Color_Scheme,
                              rect_gp = gpar(col = "white", lwd = 0.1),
                              cluster_rows = Clustering,
                              cluster_columns = Wash_Water_Clust,
                              show_row_names = FALSE,
                              show_column_names = FALSE,
                              show_column_dend = FALSE,
                              show_row_dend = FALSE,
                              #column_title = "Wash Water",
                              column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                              column_title_side = "bottom",
                              column_dend_height = unit(2, "cm"),
                              show_heatmap_legend = FALSE,
                              top_annotation = Wash_Water_Top_Annotations)
#heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 3)))

Wash_Water_Heatmap

Drinking_Water_Meta <- Drinking_Water_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Environmental_Drinking_Water" = "#56B4E9"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Throat" = "yellow", "Cloacal" = "pink"))
Drinking_Water_Top_Annotations = HeatmapAnnotation(df = Drinking_Water_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Drinking_Water_Sample_Type <- Drinking_Water_Metadata[c(5)]
#Drinking_Water_Bottom_Annotations = HeatmapAnnotation(df = Drinking_Water_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Drinking_Water_Heatmap <- Heatmap(t(Drinking_Water_Mat),
                                  name = "Rank(RPM)",
                                  col = Color_Scheme,
                                  rect_gp = gpar(col = "white", lwd = 0.1),
                                  cluster_rows = Clustering,
                                  cluster_columns = Drinking_Water_Clust,
                                  show_row_names = FALSE,
                                  show_column_names = FALSE,
                                  show_column_dend = FALSE,
                                  show_row_dend = FALSE,
                                  #column_title = "Drinking Water",
                                  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                                  column_title_side = "bottom",
                                  column_dend_height = unit(2, "cm"),
                                  show_heatmap_legend = FALSE,
                                  top_annotation = Drinking_Water_Top_Annotations)
#heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 3)))

Drinking_Water_Heatmap

Cage_Meta <- Cage_Metadata[c(1)]
Colors_Top <- list("Location_Grouping" = c("Environmental_Cage" = "#D55E00"))
#Colors_Bottom <- list("Alternative_Grouping" = c("Throat" = "yellow", "Cloacal" = "pink"))
Cage_Top_Annotations = HeatmapAnnotation(df = Cage_Meta, show_annotation_name = FALSE, col = Colors_Top)
#Cage_Sample_Type <- Cage_Metadata[c(5)]
#Cage_Bottom_Annotations = HeatmapAnnotation(df = Cage_Sample_Type, show_annotation_name = FALSE, col = Colors_Bottom)

Cage_Heatmap <- Heatmap(t(Cage_Mat),
                        name = "Rank(RPM)",
                        col = Color_Scheme,
                        rect_gp = gpar(col = "white", lwd = 0.1),
                        cluster_rows = Clustering,
                        cluster_columns = Cage_Clust,
                        show_row_names = TRUE,
                        row_names_gp = gpar(fontsize = 6, fontface = "bold"),
                        show_column_names = FALSE,
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        #column_title = "Cage",
                        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                        column_title_side = "bottom",
                        column_dend_height = unit(2, "cm"),
                        show_heatmap_legend = TRUE,
                        top_annotation = Cage_Top_Annotations,
                        heatmap_legend_param = list(legend_gp = gpar(col = "black", lwd = 1.5)))

Cage_Heatmap

#================================================
#### Combine and Plot
#================================================
Heatmap <- Neg_Heatmap + Chicken_Throat_Heatmap + Chicken_Cloacal_Heatmap + 
  Duck_Throat_Heatmap + Duck_Cloacal_Heatmap + Air_Slaughter_Heatmap + Air_Holding_Heatmap + Air_Outside_Heatmap +
  Wash_Water_Heatmap + Drinking_Water_Heatmap + Cage_Heatmap

pdf("Negative_Controls_Heatmap.pdf", width = 16, height = 12)  
draw(Heatmap)
dev.off()
