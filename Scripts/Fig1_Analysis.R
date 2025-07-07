#================================================
#### Packages
#================================================
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(ade4)
library(vegan)
library(gplots)
library(gridExtra)
library(grid)
library(ape)
library(phyloseq)
library(ggforce)
library(plyr)
library(psych)
library(dplyr)
library(tidyverse)
library(data.table)
library(compositions)
library(devtools)
library(export)
library(emmeans)
library(lme4)
library(simr)
library(reshape)
library(ComplexHeatmap)
library(circlize)
library(PMCMRplus)
library(ggrepel)
library(ggdist)
library(janitor)
library(lefser)
library(cowplot)
library(openxlsx)
library(sf)
library(ggpattern)
library(ggridges)
library(lsa)
library(patchwork)
library(viridis)
library(ggrastr)
library(gghalves)
library(ggbeeswarm)
library(ggtree)
library(ape)
library(treeio)
library(tidytree)
library(ggplot2)
library(openxlsx)
library(ggnewscale)
library(ggstream)
library(RColorBrewer)
library(tidyr)
library(sf)
library(ggplot2)
library(dplyr)
library(gggenes)
library(IRanges)

#================================================
#### Colour Palette For Analysis
#================================================
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

#================================================
#### Alpha-Diversity - Poultry Viruses
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Species <- (Species > 0) * 1
Species <- t(Species)
Species <- as.data.frame(Species)

ASV_Alpha <- otu_table(Species, taxa_are_rows = TRUE)
Metadata <- sample_data(Metadata)
physeq <- phyloseq(ASV_Alpha, Metadata)

Alpha_Diversity <- estimate_richness(
  physeq, split = TRUE,
  measures = c("Observed", "Shannon", "Simpson", "Chao1")
)

Alpha_Diversity <- cbind(Alpha_Diversity, as.data.frame(sample_data(physeq)))
Alpha_Diversity$Location_Grouping <- as.factor(Alpha_Diversity$Location_Grouping)

Kruskall_Observed <- kwAllPairsDunnTest(Alpha_Diversity$Observed, Alpha_Diversity$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Observed$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Alpha_Diversity_Observed.xlsx", rowNames = TRUE)
Kruskall_Chao1 <- kwAllPairsDunnTest(Alpha_Diversity$Chao1, Alpha_Diversity$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chao1$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Alpha_Diversity_Chao1.xlsx", rowNames = TRUE)
Kruskall_Shannon <- kwAllPairsDunnTest(Alpha_Diversity$Shannon, Alpha_Diversity$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Shannon$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Alpha_Diversity_Shannon.xlsx", rowNames = TRUE)
Kruskall_Simpson <- kwAllPairsDunnTest(Alpha_Diversity$Simpson, Alpha_Diversity$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Simpson$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Alpha_Diversity_Simpson.xlsx", rowNames = TRUE)

alpha_long <- pivot_longer(
  Alpha_Diversity,
  cols = c("Observed", "Shannon", "Simpson", "Chao1"),
  names_to = "Metric",
  values_to = "Value"
)

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
  "Cage" = "#D55E00"
)

alpha_long$Location_Grouping <- factor(alpha_long$Location_Grouping, levels = names(group_palette))

plot_list <- alpha_long %>%
  split(.$Metric) %>%
  map(~ ggplot(.x, aes(x = Location_Grouping, y = Value, fill = Location_Grouping)) +
        geom_boxplot(outlier.shape = NA, color = "black") +
        geom_jitter(aes(color = Location_Grouping), width = 0.2, alpha = 0.6, size = 1.5) +
        scale_fill_manual(values = group_palette) +
        scale_color_manual(values = group_palette) +
        theme_minimal() +
        labs(
          title = unique(.x$Metric),
          x = NULL,
          y = "Alpha Diversity"
        ) +
        theme(
          axis.text.x = element_blank(),      
          axis.ticks.x = element_blank(),  
          axis.ticks.y = element_line(color = "black"),
          axis.line.x = element_line(color = "black"),  
          axis.line.y = element_line(color = "black"), 
          axis.text.y = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          plot.title = element_text(hjust = 0, face = "bold"),
          strip.text = element_text(face = "bold"),
          panel.border = element_rect(color = "grey50", fill = NA, size = 1),
          legend.position = "right"
        ))

final_alpha_plot <- plot_list$Chao1 + plot_list$Shannon + plot_list$Simpson +
  plot_layout(nrow = 3)
final_alpha_plot
ggsave(filename="Alpha_Diversity_Boxplots.pdf", plot = final_alpha_plot, width=200, height=280, units="mm")
ggsave(filename="Observed_Species_Boxplots.pdf", plot = plot_list$Observed, width=200, height=140, units="mm")

#================================================
#### Avian Species Detection Heatmap - Summary
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Binary_Species <- (Species > 0) * 1

avian_groups <- c("Chicken_Throat", "Chicken_Cloacal", "Duck_Throat", "Duck_Cloacal")
all_experiments <- unique(Metadata$Experiment)
all_viruses <- colnames(Binary_Species)

summary_matrix <- matrix(0, nrow = length(all_experiments), ncol = length(all_viruses))
rownames(summary_matrix) <- paste0("Experiment_", all_experiments)
colnames(summary_matrix) <- all_viruses

for (exp_id in all_experiments) {
  exp_samples <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  avian_samples <- rownames(Metadata[Metadata$Experiment == exp_id &
                                       Metadata$Location_Grouping %in% avian_groups, ])
  env_samples <- setdiff(exp_samples, avian_samples)
  
  if (length(avian_samples) == 0 || length(env_samples) == 0) next
  
  avian_data <- Binary_Species[avian_samples, , drop = FALSE]
  env_data <- Binary_Species[env_samples, , drop = FALSE]
  
  avian_detected <- colSums(avian_data) > 0
  env_detected <- colSums(env_data) > 0
  
  row_label <- paste0("Experiment_", exp_id)
  
  for (virus in all_viruses) {
    a <- avian_detected[virus]
    e <- env_detected[virus]
    
    if (a && e) {
      summary_matrix[row_label, virus] <- 1
    }
    if (a && !e) {
      summary_matrix[row_label, virus] <- 2
    }
    if (!a && e) {
      summary_matrix[row_label, virus] <- 3
    }
  }
}

Label_Colors <- colorRamp2(
  c(0, 1, 2, 3),
  c("white", "#F8766D", "#F6E09B", "#00BFC4")
)

status_legend <- Legend(
  title = "Detection",
  at = c(1, 2, 3),
  labels = c("Detected", "Undetected", "Unique"),
  legend_gp = gpar(fill = Label_Colors(c(1, 2, 3))),
  border = "black",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

experiment_order <- paste0("Experiment_", sort(as.numeric(gsub("Experiment_", "", rownames(summary_matrix)))))
summary_matrix <- summary_matrix[experiment_order, ]

summary_heatmap <- Heatmap(
  summary_matrix,
  col = Label_Colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontface = "bold", fontsize = 9),
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontface = "bold", fontsize = 7),
  rect_gp = gpar(col = "black", lwd = 0.4),
  column_title = "",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  show_heatmap_legend = FALSE
)

Summary_Final_Heatmap <- draw(
  summary_heatmap,
  padding = unit(c(10, 20, 10, 10), "mm"),
  heatmap_legend_list = list(status_legend),
  heatmap_legend_side = "right"
)

pdf("Summary_Heatmap_A4.pdf", width = 15, height = 4)  
draw(Summary_Final_Heatmap)
dev.off()

#================================================
#### Avian Species Detection - Line Graph Summary
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Binary_Species <- (Species > 0) * 1

avian_groups <- c("Chicken_Throat", "Chicken_Cloacal", "Duck_Throat", "Duck_Cloacal")
all_experiments <- sort(unique(Metadata$Experiment))
all_viruses <- colnames(Binary_Species)

summary_matrix <- matrix(0, nrow = length(all_experiments), ncol = length(all_viruses))
rownames(summary_matrix) <- paste0("Experiment_", all_experiments)
colnames(summary_matrix) <- all_viruses

for (exp_id in all_experiments) {
  exp_samples <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  avian_samples <- rownames(Metadata[Metadata$Experiment == exp_id &
                                       Metadata$Location_Grouping %in% avian_groups, ])
  env_samples <- setdiff(exp_samples, avian_samples)
  
  if (length(avian_samples) == 0 || length(env_samples) == 0) next
  
  avian_data <- Binary_Species[avian_samples, , drop = FALSE]
  env_data <- Binary_Species[env_samples, , drop = FALSE]
  
  avian_detected <- colSums(avian_data) > 0
  env_detected <- colSums(env_data) > 0
  
  row_label <- paste0("Experiment_", exp_id)
  
  for (virus in all_viruses) {
    a <- avian_detected[virus]
    e <- env_detected[virus]
    
    if (a && e) {
      summary_matrix[row_label, virus] <- 1  
    }
    if (a && !e) {
      summary_matrix[row_label, virus] <- 2  
    }
    if (!a && e) {
      summary_matrix[row_label, virus] <- 3  
    }
  }
}

summary_df <- as.data.frame(summary_matrix) %>%
  rownames_to_column("Experiment") %>%
  pivot_longer(-Experiment, names_to = "Virus", values_to = "Status")

count_df <- summary_df %>%
  group_by(Experiment, Status) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Status = factor(Status, levels = c(1, 2, 3),
                         labels = c("Detected", "Undetected", "Unique")))

wide_df <- count_df %>%
  pivot_wider(names_from = Status, values_from = Count, values_fill = 0) %>%
  mutate(
    Total = Detected + Undetected
  ) %>%
  filter(Total > 0) %>%
  mutate(
    Percent_Detected = (Detected / Total) * 100
  )

main_df <- wide_df %>%
  pivot_longer(cols = c("Detected", "Total"), names_to = "Status", values_to = "Count")

unique_df <- wide_df %>%
  select(Experiment, Unique)

percent_df <- wide_df %>%
  select(Experiment, Percent_Detected)

ordered_experiments <- paste0("Experiment_", sort(as.numeric(gsub("Experiment_", "", wide_df$Experiment))))
main_df$Experiment <- factor(main_df$Experiment, levels = ordered_experiments)
unique_df$Experiment <- factor(unique_df$Experiment, levels = ordered_experiments)
percent_df$Experiment <- factor(percent_df$Experiment, levels = ordered_experiments)

status_colors <- c(
  "Detected" = "#F8766D",
  "Total" = "#228B22",
  "Unique" = "#00BFC4"
)
main_df <- main_df %>%
  mutate(Status = factor(Status, levels = c("Total", "Detected")))

main_plot <- ggplot(main_df, aes(x = Experiment, y = Count, fill = Status, group = Status)) +
  geom_area(position = "identity", alpha = 0.8, size = 1.2, color = "black", lineend = "round") +
  geom_line(data = unique_df, aes(x = Experiment, y = Unique, group = 1),
            inherit.aes = FALSE, color = status_colors["Unique"], size = 1.5) +
  scale_fill_manual(values = status_colors[c("Total", "Detected")]) +
  labs(
    title = "Poultry Virus Detection Summary",
    x = "Experiment",
    y = "Virus Species Detected",
    fill = "Detection Type"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title.x = element_text(size = 12, color = "black", face = "bold"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    plot.title = element_text(hjust = 0, size = 12, face = "bold", color = "black"),
    legend.position = "right"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

percentage_plot <- ggplot(percent_df, aes(x = Experiment, y = Percent_Detected, group = 1)) +
  geom_area(fill = "#F8766D", alpha = 0.6, color = "black", size = 1.2) +
  labs(
    title = "Percentage of Detected Viruses",
    x = "Experiment",
    y = "% Detected"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", color = "black"),
    axis.title.y = element_text(size = 11, face = "bold", color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0, color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(size = 1, color = "black"),
    axis.line.y = element_line(size = 1, color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

main_plot
percentage_plot

ggsave(filename="Line_Area_Graph_Summary.pdf", plot = main_plot, width=200, height=140, units="mm")
ggsave(filename="Line_Area_Percentage_Summary.pdf", plot = percentage_plot, width=200, height=140, units="mm")

#================================================
#### Avian Species Detection Heatmap - Chicken Throat
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Chicken_Throat <- subset(Metadata, Location_Grouping == "Chicken_Throat")
Chicken_Throat_Samples <- rownames(Chicken_Throat)
Chicken_Cloacal <- subset(Metadata, Location_Grouping == "Chicken_Cloacal")
Chicken_Cloacal_Samples <- rownames(Chicken_Cloacal)
Duck_Cloacal <- subset(Metadata, Location_Grouping == "Duck_Cloacal")
Duck_Cloacal_Samples <- rownames(Duck_Cloacal)
Duck_Throat <- subset(Metadata, Location_Grouping == "Duck_Throat")
Duck_Throat_Samples <- rownames(Duck_Throat)

Chicken_Throat_Species <- Species[rownames(Species) %in% Chicken_Throat_Samples, ]
Chicken_Throat_Species <- Chicken_Throat_Species[, colSums(Chicken_Throat_Species != 0) > 0]
Chicken_Throat_Taxa <- colnames(Chicken_Throat_Species)

Chicken_Throat_Data <- Species[, colnames(Species) %in% Chicken_Throat_Taxa]
Chicken_Throat_Data <- Chicken_Throat_Data[!rownames(Chicken_Throat_Data) %in% Chicken_Cloacal_Samples, ]
Chicken_Throat_Data <- Chicken_Throat_Data[!rownames(Chicken_Throat_Data) %in% Duck_Cloacal_Samples, ]
Chicken_Throat_Data <- Chicken_Throat_Data[!rownames(Chicken_Throat_Data) %in% Duck_Throat_Samples, ]

Chicken_Throat_Meta <- Metadata[rownames(Metadata) %in% rownames(Chicken_Throat_Data), ]

Binary_Data <- (Chicken_Throat_Data > 0) * 1
meta <- Metadata[rownames(Binary_Data), ]

Labelled_Matrix <- matrix(NA, nrow = nrow(Binary_Data), ncol = ncol(Binary_Data))
rownames(Labelled_Matrix) <- rownames(Binary_Data)
colnames(Labelled_Matrix) <- colnames(Binary_Data)

for (exp_id in unique(meta$Experiment)) {
  exp_samples <- rownames(meta[meta$Experiment == exp_id, ])
  CT_samples <- rownames(meta[meta$Experiment == exp_id & meta$Location_Grouping == "Chicken_Throat", ])
  other_samples <- setdiff(exp_samples, CT_samples)
  
  if (length(CT_samples) == 0 || length(other_samples) == 0) next
  
  CT_block <- Binary_Data[CT_samples, , drop = FALSE]
  CT_consensus <- colSums(CT_block) > 0
  
  for (sample in other_samples) {
    for (taxon in colnames(Binary_Data)) {
      ct_val <- CT_consensus[taxon] * 1
      other_val <- Binary_Data[sample, taxon]
      
      if (ct_val == 1 && other_val == 1) {
        Labelled_Matrix[sample, taxon] <- 10
      } else if (ct_val == 1 && other_val == 0) {
        Labelled_Matrix[sample, taxon] <- 20
      } else if (ct_val == 0 && other_val == 1) {
        Labelled_Matrix[sample, taxon] <- 30
      } else {
        Labelled_Matrix[sample, taxon] <- 0
      }
    }
  }
}

Labelled_Matrix <- Labelled_Matrix[meta$Location_Grouping != "Chicken_Throat", ]
meta_filtered <- meta[rownames(Labelled_Matrix), ]

Label_Colors <- colorRamp2(
  c(0, 10, 20, 30),
  c("white", "#F8766D", "#F6E09B", "#00BFC4")
)

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
  "Cage" = "#D55E00"
)

location_colors <- c("Phnom_Penh" = "grey", "Takeo" = "black")

all_taxa <- colnames(Labelled_Matrix)
experiments <- unique(meta_filtered$Experiment)

heatmap_list <- lapply(experiments, function(exp_id) {
  samples <- rownames(meta_filtered[meta_filtered$Experiment == exp_id, ])
  sub_matrix <- t(Labelled_Matrix[samples, , drop = FALSE]) 
  
  missing_rows <- setdiff(all_taxa, rownames(sub_matrix))
  if (length(missing_rows) > 0) {
    for (v in missing_rows) {
      sub_matrix <- rbind(sub_matrix, setNames(rep(0, ncol(sub_matrix)), colnames(sub_matrix)))
      rownames(sub_matrix)[nrow(sub_matrix)] <- v
    }
  }
  sub_matrix <- sub_matrix[all_taxa, , drop = FALSE]
  
  anno_subset <- meta_filtered[samples, ]
  
  top_anno <- HeatmapAnnotation(
    "Sample Type" = anno_subset$Location_Grouping,
    col = list("Sample Type" = group_palette),
    show_annotation_name = FALSE
  )
  
  bottom_anno <- HeatmapAnnotation(
    Location = anno_subset$Location,
    col = list(Location = location_colors),
    show_annotation_name = FALSE,
    annotation_name_side = "right",
    which = "column"
  )
  
  Heatmap(
    sub_matrix,
    col = Label_Colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontface = "bold", fontsize = 6),
    show_column_names = FALSE,
    rect_gp = gpar(col = "black", lwd = 0.5),
    top_annotation = top_anno,
    bottom_annotation = bottom_anno,
    name = NULL,
    column_title = paste("Experiment", exp_id),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_title_side = "bottom",
    show_heatmap_legend = FALSE
  )
})

status_legend <- Legend(
  title = "Detection",
  at = c(10, 20, 30),
  labels = c("Detected", "Undetected", "ES"),
  legend_gp = gpar(fill = Label_Colors(c(10, 20, 30))),
  border = "black",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

Chicken_Throat_Final_Heatmap <- draw(
  Reduce(`+`, heatmap_list),
  padding = unit(c(10, 10, 10, 10), "mm"),
  heatmap_legend_list = list(status_legend),
  heatmap_legend_side = "right"
)

pdf("Chicken_Throat_Heatmap_A4.pdf", width = 10, height = 4)  # A4 size in inches
draw(Chicken_Throat_Final_Heatmap)
dev.off()

#================================================
#### Avian Species Detection Heatmap - Chicken Cloacal
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Chicken_Throat <- subset(Metadata, Location_Grouping == "Chicken_Throat")
Chicken_Throat_Samples <- rownames(Chicken_Throat)
Chicken_Cloacal <- subset(Metadata, Location_Grouping == "Chicken_Cloacal")
Chicken_Cloacal_Samples <- rownames(Chicken_Cloacal)
Duck_Cloacal <- subset(Metadata, Location_Grouping == "Duck_Cloacal")
Duck_Cloacal_Samples <- rownames(Duck_Cloacal)
Duck_Throat <- subset(Metadata, Location_Grouping == "Duck_Throat")
Duck_Throat_Samples <- rownames(Duck_Throat)

Ref_Group <- "Chicken_Cloacal"
Ref_Samples <- Chicken_Cloacal_Samples

Ref_Species <- Species[Ref_Samples, ]
Ref_Species <- Ref_Species[, colSums(Ref_Species != 0) > 0]
Ref_Taxa <- colnames(Ref_Species)

Ref_Data <- Species[, colnames(Species) %in% Ref_Taxa]
Ref_Data <- Ref_Data[!rownames(Ref_Data) %in% Chicken_Throat_Samples, ]
Ref_Data <- Ref_Data[!rownames(Ref_Data) %in% Duck_Cloacal_Samples, ]
Ref_Data <- Ref_Data[!rownames(Ref_Data) %in% Duck_Throat_Samples, ]

Ref_Meta <- Metadata[rownames(Metadata) %in% rownames(Ref_Data), ]

Binary_Data <- (Ref_Data > 0) * 1
meta <- Metadata[rownames(Binary_Data), ]

Labelled_Matrix <- matrix(NA, nrow = nrow(Binary_Data), ncol = ncol(Binary_Data))
rownames(Labelled_Matrix) <- rownames(Binary_Data)
colnames(Labelled_Matrix) <- colnames(Binary_Data)

for (exp_id in unique(meta$Experiment)) {
  exp_samples <- rownames(meta[meta$Experiment == exp_id, ])
  Ref_Exp_Samples <- rownames(meta[meta$Experiment == exp_id & meta$Location_Grouping == Ref_Group, ])
  other_samples <- setdiff(exp_samples, Ref_Exp_Samples)
  
  if (length(Ref_Exp_Samples) == 0 || length(other_samples) == 0) next
  
  Ref_Block <- Binary_Data[Ref_Exp_Samples, , drop = FALSE]
  Ref_Consensus <- colSums(Ref_Block) > 0
  
  for (sample in other_samples) {
    for (taxon in colnames(Binary_Data)) {
      ref_val <- Ref_Consensus[taxon] * 1
      other_val <- Binary_Data[sample, taxon]
      
      if (ref_val == 1 && other_val == 1) {
        Labelled_Matrix[sample, taxon] <- 10
      } else if (ref_val == 1 && other_val == 0) {
        Labelled_Matrix[sample, taxon] <- 20
      } else if (ref_val == 0 && other_val == 1) {
        Labelled_Matrix[sample, taxon] <- 30
      } else {
        Labelled_Matrix[sample, taxon] <- 0
      }
    }
  }
}

Labelled_Matrix <- Labelled_Matrix[meta$Location_Grouping != Ref_Group, ]
meta_filtered <- meta[rownames(Labelled_Matrix), ]

Label_Colors <- colorRamp2(
  c(0, 10, 20, 30),
  c("white", "#F8766D", "#F6E09B", "#00BFC4")
)

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
  "Cage" = "#D55E00"
)

location_colors <- c("Phnom_Penh" = "grey", "Takeo" = "black")

all_taxa <- colnames(Labelled_Matrix)
experiments <- unique(meta_filtered$Experiment)

heatmap_list <- lapply(experiments, function(exp_id) {
  samples <- rownames(meta_filtered[meta_filtered$Experiment == exp_id, ])
  sub_matrix <- t(Labelled_Matrix[samples, , drop = FALSE]) 
  
  missing_rows <- setdiff(all_taxa, rownames(sub_matrix))
  if (length(missing_rows) > 0) {
    for (v in missing_rows) {
      sub_matrix <- rbind(sub_matrix, setNames(rep(0, ncol(sub_matrix)), colnames(sub_matrix)))
      rownames(sub_matrix)[nrow(sub_matrix)] <- v
    }
  }
  sub_matrix <- sub_matrix[all_taxa, , drop = FALSE]
  
  anno_subset <- meta_filtered[samples, ]
  
  top_anno <- HeatmapAnnotation(
    "Sample Type" = anno_subset$Location_Grouping,
    col = list("Sample Type" = group_palette),
    show_annotation_name = FALSE
  )
  
  bottom_anno <- HeatmapAnnotation(
    Location = anno_subset$Location,
    col = list(Location = location_colors),
    show_annotation_name = FALSE,
    annotation_name_side = "right",
    which = "column"
  )
  
  Heatmap(
    sub_matrix,
    col = Label_Colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontface = "bold", fontsize = 6),
    show_column_names = FALSE,
    rect_gp = gpar(col = "black", lwd = 0.5),
    top_annotation = top_anno,
    bottom_annotation = bottom_anno,
    name = NULL,
    column_title = paste("Experiment", exp_id),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_title_side = "bottom",
    show_heatmap_legend = FALSE
  )
})

status_legend <- Legend(
  title = "Detection",
  at = c(10, 20, 30),
  labels = c("Detected", "Undetected", "ES"),
  legend_gp = gpar(fill = Label_Colors(c(10, 20, 30))),
  border = "black",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

Chicken_Cloacal_Final_Heatmap <- draw(
  Reduce(`+`, heatmap_list),
  padding = unit(c(10, 10, 10, 10), "mm"),
  heatmap_legend_list = list(status_legend),
  heatmap_legend_side = "right"
)

pdf("Chicken_Cloacal_Heatmap_A4.pdf", width = 10, height = 6)  # A4 size in inches
draw(Chicken_Cloacal_Final_Heatmap)
dev.off()

#================================================
#### Avian Species Detection Heatmap - Duck Throat
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Chicken_Throat <- subset(Metadata, Location_Grouping == "Chicken_Throat")
Chicken_Throat_Samples <- rownames(Chicken_Throat)
Chicken_Cloacal <- subset(Metadata, Location_Grouping == "Chicken_Cloacal")
Chicken_Cloacal_Samples <- rownames(Chicken_Cloacal)
Duck_Cloacal <- subset(Metadata, Location_Grouping == "Duck_Cloacal")
Duck_Cloacal_Samples <- rownames(Duck_Cloacal)
Duck_Throat <- subset(Metadata, Location_Grouping == "Duck_Throat")
Duck_Throat_Samples <- rownames(Duck_Throat)

Ref_Group <- "Duck_Throat"
Ref_Samples <- Duck_Throat_Samples

Ref_Species <- Species[Ref_Samples, ]
Ref_Species <- Ref_Species[, colSums(Ref_Species != 0) > 0]
Ref_Taxa <- colnames(Ref_Species)

Ref_Data <- Species[, colnames(Species) %in% Ref_Taxa]
Ref_Data <- Ref_Data[!rownames(Ref_Data) %in% Chicken_Throat_Samples, ]
Ref_Data <- Ref_Data[!rownames(Ref_Data) %in% Chicken_Cloacal_Samples, ]
Ref_Data <- Ref_Data[!rownames(Ref_Data) %in% Duck_Cloacal_Samples, ]

Ref_Meta <- Metadata[rownames(Metadata) %in% rownames(Ref_Data), ]

Binary_Data <- (Ref_Data > 0) * 1
meta <- Metadata[rownames(Binary_Data), ]

Labelled_Matrix <- matrix(NA, nrow = nrow(Binary_Data), ncol = ncol(Binary_Data))
rownames(Labelled_Matrix) <- rownames(Binary_Data)
colnames(Labelled_Matrix) <- colnames(Binary_Data)

for (exp_id in unique(meta$Experiment)) {
  exp_samples <- rownames(meta[meta$Experiment == exp_id, ])
  Ref_Exp_Samples <- rownames(meta[meta$Experiment == exp_id & meta$Location_Grouping == Ref_Group, ])
  other_samples <- setdiff(exp_samples, Ref_Exp_Samples)
  
  if (length(Ref_Exp_Samples) == 0 || length(other_samples) == 0) next
  
  Ref_Block <- Binary_Data[Ref_Exp_Samples, , drop = FALSE]
  Ref_Consensus <- colSums(Ref_Block) > 0
  
  for (sample in other_samples) {
    for (taxon in colnames(Binary_Data)) {
      ref_val <- Ref_Consensus[taxon] * 1
      other_val <- Binary_Data[sample, taxon]
      
      if (ref_val == 1 && other_val == 1) {
        Labelled_Matrix[sample, taxon] <- 10
      } else if (ref_val == 1 && other_val == 0) {
        Labelled_Matrix[sample, taxon] <- 20
      } else if (ref_val == 0 && other_val == 1) {
        Labelled_Matrix[sample, taxon] <- 30
      } else {
        Labelled_Matrix[sample, taxon] <- 0
      }
    }
  }
}

Labelled_Matrix <- Labelled_Matrix[meta$Location_Grouping != Ref_Group, ]
meta_filtered <- meta[rownames(Labelled_Matrix), ]

Label_Colors <- colorRamp2(
  c(0, 10, 20, 30),
  c("white", "#F8766D", "#F6E09B", "#00BFC4")
)

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
  "Cage" = "#D55E00"
)

location_colors <- c("Phnom_Penh" = "grey", "Takeo" = "black")

all_taxa <- colnames(Labelled_Matrix)
experiments <- unique(meta_filtered$Experiment)

heatmap_list <- lapply(experiments, function(exp_id) {
  samples <- rownames(meta_filtered[meta_filtered$Experiment == exp_id, ])
  sub_matrix <- t(Labelled_Matrix[samples, , drop = FALSE]) 
  
  missing_rows <- setdiff(all_taxa, rownames(sub_matrix))
  if (length(missing_rows) > 0) {
    for (v in missing_rows) {
      sub_matrix <- rbind(sub_matrix, setNames(rep(0, ncol(sub_matrix)), colnames(sub_matrix)))
      rownames(sub_matrix)[nrow(sub_matrix)] <- v
    }
  }
  sub_matrix <- sub_matrix[all_taxa, , drop = FALSE]
  
  anno_subset <- meta_filtered[samples, ]
  
  top_anno <- HeatmapAnnotation(
    "Sample Type" = anno_subset$Location_Grouping,
    col = list("Sample Type" = group_palette),
    show_annotation_name = FALSE
  )
  
  bottom_anno <- HeatmapAnnotation(
    Location = anno_subset$Location,
    col = list(Location = location_colors),
    show_annotation_name = FALSE,
    annotation_name_side = "right",
    which = "column"
  )
  
  Heatmap(
    sub_matrix,
    col = Label_Colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontface = "bold", fontsize = 6),
    show_column_names = FALSE,
    rect_gp = gpar(col = "black", lwd = 0.5),
    top_annotation = top_anno,
    bottom_annotation = bottom_anno,
    name = NULL,
    column_title = paste("Experiment", exp_id),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_title_side = "bottom",
    show_heatmap_legend = FALSE
  )
})

status_legend <- Legend(
  title = "Detection",
  at = c(10, 20, 30),
  labels = c("Detected", "Undetected", "ES"),
  legend_gp = gpar(fill = Label_Colors(c(10, 20, 30))),
  border = "black",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

Duck_Throat_Final_Heatmap <- draw(
  Reduce(`+`, heatmap_list),
  padding = unit(c(10, 10, 10, 10), "mm"),
  heatmap_legend_list = list(status_legend),
  heatmap_legend_side = "right"
)

pdf("Duck_Throat_Heatmap_A4.pdf", width = 10, height = 4)  # A4 size in inches
draw(Duck_Throat_Final_Heatmap)
dev.off()

#================================================
#### Avian Species Detection Heatmap - Duck Cloacal
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Chicken_Throat <- subset(Metadata, Location_Grouping == "Chicken_Throat")
Chicken_Throat_Samples <- rownames(Chicken_Throat)
Chicken_Cloacal <- subset(Metadata, Location_Grouping == "Chicken_Cloacal")
Chicken_Cloacal_Samples <- rownames(Chicken_Cloacal)
Duck_Cloacal <- subset(Metadata, Location_Grouping == "Duck_Cloacal")
Duck_Cloacal_Samples <- rownames(Duck_Cloacal)
Duck_Throat <- subset(Metadata, Location_Grouping == "Duck_Throat")
Duck_Throat_Samples <- rownames(Duck_Throat)

Ref_Group <- "Duck_Cloacal"
Ref_Samples <- Duck_Cloacal_Samples

Ref_Species <- Species[Ref_Samples, ]
Ref_Species <- Ref_Species[, colSums(Ref_Species != 0) > 0]
Ref_Taxa <- colnames(Ref_Species)

Ref_Data <- Species[, colnames(Species) %in% Ref_Taxa]
Ref_Data <- Ref_Data[!rownames(Ref_Data) %in% Chicken_Throat_Samples, ]
Ref_Data <- Ref_Data[!rownames(Ref_Data) %in% Chicken_Cloacal_Samples, ]
Ref_Data <- Ref_Data[!rownames(Ref_Data) %in% Duck_Throat_Samples, ]

Ref_Meta <- Metadata[rownames(Metadata) %in% rownames(Ref_Data), ]

Binary_Data <- (Ref_Data > 0) * 1
meta <- Metadata[rownames(Binary_Data), ]

Labelled_Matrix <- matrix(NA, nrow = nrow(Binary_Data), ncol = ncol(Binary_Data))
rownames(Labelled_Matrix) <- rownames(Binary_Data)
colnames(Labelled_Matrix) <- colnames(Binary_Data)

for (exp_id in unique(meta$Experiment)) {
  exp_samples <- rownames(meta[meta$Experiment == exp_id, ])
  Ref_Exp_Samples <- rownames(meta[meta$Experiment == exp_id & meta$Location_Grouping == Ref_Group, ])
  other_samples <- setdiff(exp_samples, Ref_Exp_Samples)
  
  if (length(Ref_Exp_Samples) == 0 || length(other_samples) == 0) next
  
  Ref_Block <- Binary_Data[Ref_Exp_Samples, , drop = FALSE]
  Ref_Consensus <- colSums(Ref_Block) > 0
  
  for (sample in other_samples) {
    for (taxon in colnames(Binary_Data)) {
      ref_val <- Ref_Consensus[taxon] * 1
      other_val <- Binary_Data[sample, taxon]
      
      if (ref_val == 1 && other_val == 1) {
        Labelled_Matrix[sample, taxon] <- 10
      } else if (ref_val == 1 && other_val == 0) {
        Labelled_Matrix[sample, taxon] <- 20
      } else if (ref_val == 0 && other_val == 1) {
        Labelled_Matrix[sample, taxon] <- 30
      } else {
        Labelled_Matrix[sample, taxon] <- 0
      }
    }
  }
}

Labelled_Matrix <- Labelled_Matrix[meta$Location_Grouping != Ref_Group, ]
meta_filtered <- meta[rownames(Labelled_Matrix), ]

Label_Colors <- colorRamp2(
  c(0, 10, 20, 30),
  c("white", "#F8766D", "#F6E09B", "#00BFC4")
)

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
  "Cage" = "#D55E00"
)

location_colors <- c("Phnom_Penh" = "grey", "Takeo" = "black")

all_taxa <- colnames(Labelled_Matrix)
experiments <- unique(meta_filtered$Experiment)

heatmap_list <- lapply(experiments, function(exp_id) {
  samples <- rownames(meta_filtered[meta_filtered$Experiment == exp_id, ])
  sub_matrix <- t(Labelled_Matrix[samples, , drop = FALSE]) 
  
  missing_rows <- setdiff(all_taxa, rownames(sub_matrix))
  if (length(missing_rows) > 0) {
    for (v in missing_rows) {
      sub_matrix <- rbind(sub_matrix, setNames(rep(0, ncol(sub_matrix)), colnames(sub_matrix)))
      rownames(sub_matrix)[nrow(sub_matrix)] <- v
    }
  }
  sub_matrix <- sub_matrix[all_taxa, , drop = FALSE]
  
  anno_subset <- meta_filtered[samples, ]
  
  top_anno <- HeatmapAnnotation(
    "Sample Type" = anno_subset$Location_Grouping,
    col = list("Sample Type" = group_palette),
    show_annotation_name = FALSE
  )
  
  bottom_anno <- HeatmapAnnotation(
    Location = anno_subset$Location,
    col = list(Location = location_colors),
    show_annotation_name = FALSE,
    annotation_name_side = "right",
    which = "column"
  )
  
  Heatmap(
    sub_matrix,
    col = Label_Colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontface = "bold", fontsize = 6),
    show_column_names = FALSE,
    rect_gp = gpar(col = "black", lwd = 0.5),
    top_annotation = top_anno,
    bottom_annotation = bottom_anno,
    name = NULL,
    column_title = paste("Experiment", exp_id),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_title_side = "bottom",
    show_heatmap_legend = FALSE
  )
})

status_legend <- Legend(
  title = "Detection",
  at = c(10, 20, 30),
  labels = c("Detected", "Undetected", "ES"),
  legend_gp = gpar(fill = Label_Colors(c(10, 20, 30))),
  border = "black",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

Duck_Cloacal_Final_Heatmap <- draw(
  Reduce(`+`, heatmap_list),
  padding = unit(c(10, 10, 10, 10), "mm"),
  heatmap_legend_list = list(status_legend),
  heatmap_legend_side = "right"
)

pdf("Duck_Cloacal_Heatmap_A4.pdf", width = 10, height = 6)
draw(Duck_Cloacal_Final_Heatmap)
dev.off()

#================================================
#### Chicken Throat Summary - Barplot/Boxplot
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Chicken_Throat_Samples <- rownames(Metadata[Metadata$Location_Grouping == "Chicken_Throat", ])

Groups <- list(
  Air_Slaughter = rownames(Metadata[Metadata$Location_Grouping == "Air_Slaughter", ]),
  Air_Holding = rownames(Metadata[Metadata$Location_Grouping == "Air_Holding", ]),
  Air_Outside = rownames(Metadata[Metadata$Location_Grouping == "Air_Outside", ]),
  Cage = rownames(Metadata[Metadata$Location_Grouping == "Cage", ]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping == "Drinking_Water", ]),
  Wash_Water = rownames(Metadata[Metadata$Location_Grouping == "Wash_Water", ])
)

Chicken_Throat_Species <- Species[Chicken_Throat_Samples, ]
Chicken_Throat_Species <- Chicken_Throat_Species[, colSums(Chicken_Throat_Species != 0) > 0]
Chicken_Throat_Taxa <- colnames(Chicken_Throat_Species)

Chicken_Throat_Data <- Species[, colnames(Species) %in% Chicken_Throat_Taxa]
Chicken_Throat_Data <- Chicken_Throat_Data[!rownames(Chicken_Throat_Data) %in% unlist(Groups), ]

Metadata_Filtered <- Metadata[rownames(Chicken_Throat_Data), ]
Binary_Data <- (Chicken_Throat_Data > 0) * 1

all_experiments <- unique(Metadata$Experiment)
summary_list <- list()

for (exp_id in all_experiments) {
  exp_samples <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  CT_samples <- intersect(exp_samples, Chicken_Throat_Samples)
  CT_block <- Binary_Data[rownames(Binary_Data) %in% CT_samples, , drop = FALSE]
  if (nrow(CT_block) == 0) next
  CT_consensus <- colSums(CT_block) > 0
  
  ct_total <- sum(CT_consensus)
  results <- c(Chicken_Throat = ct_total)
  
  for (group_name in names(Groups)) {
    group_samples <- intersect(exp_samples, Groups[[group_name]])
    group_block <- Species[group_samples, colnames(CT_block), drop = FALSE]
    group_binary <- (group_block > 0) * 1
    group_detected <- colSums(group_binary) > 0
    
    match_count <- sum(CT_consensus & group_detected)
    unique_count <- sum(!CT_consensus & group_detected)
    
    results[paste0(group_name, "_Match")] <- match_count
    results[paste0(group_name, "_Unique")] <- unique_count
  }
  
  summary_list[[as.character(exp_id)]] <- results
}

summary_df <- do.call(cbind, summary_list)
summary_df <- as.data.frame(summary_df)
colnames(summary_df) <- paste0("Experiment_", colnames(summary_df))
summary_df$Group <- rownames(summary_df)
#write.xlsx(summary_df, file = "Chicken_Throat_Species_Detected.xlsx", rowNames = TRUE)

df_match <- summary_df[grep("Chicken_Throat|_Match$", summary_df$Group), ]
df_unique <- summary_df[grep("_Unique$", summary_df$Group), ]

long_df_match <- df_match %>%
  pivot_longer(-Group, names_to = "Experiment", values_to = "Detections") %>%
  mutate(
    Group_Label = gsub("_Match", "", Group),
    Experiment_Number = as.numeric(gsub("Experiment_", "", Experiment))
  )

long_df_unique <- df_unique %>%
  pivot_longer(-Group, names_to = "Experiment", values_to = "Detections") %>%
  mutate(
    Group_Label = gsub("_Unique", "", Group),
    Experiment_Number = as.numeric(gsub("Experiment_", "", Experiment))
  )

collected_groups <- Metadata %>%
  group_by(Experiment, Location_Grouping) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Experiment = paste0("Experiment_", Experiment))

long_df_match_filtered <- long_df_match %>%
  mutate(Experiment = paste0("Experiment_", Experiment_Number)) %>%
  inner_join(collected_groups, by = c("Experiment" = "Experiment", "Group_Label" = "Location_Grouping"))

long_df_unique_filtered <- long_df_unique %>%
  mutate(Experiment = paste0("Experiment_", Experiment_Number)) %>%
  inner_join(collected_groups, by = c("Experiment" = "Experiment", "Group_Label" = "Location_Grouping"))

desired_order <- c(
  "Chicken_Throat",
  "Air_Slaughter",
  "Air_Holding",
  "Air_Outside",
  "Wash_Water",
  "Drinking_Water",
  "Cage"
)

long_df_match_filtered$Group_Label <- factor(long_df_match_filtered$Group_Label, levels = desired_order)
long_df_unique_filtered$Group_Label <- factor(long_df_unique_filtered$Group_Label, levels = desired_order[-1])

group_palette <- c(
  "Chicken_Throat" = "#009E73",
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",
  "Air_Outside" = "#228B22",
  "Cage" = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water" = "#F0E442"
)

Chicken_Throat_Detected <- ggplot(long_df_match_filtered, aes(x = Group_Label, y = Detections, fill = Group_Label)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_match_filtered$Experiment_Number)))), 
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Detected Viruses (Chicken_Throat + Matches)",
    y = "Number of Detections"
  )
Chicken_Throat_Detected

Chicken_Throat_Unique <- ggplot(long_df_unique_filtered, aes(x = Group_Label, y = Detections, fill = Group_Label)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_unique_filtered$Experiment_Number)))), 
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Unique Detections (Chicken_Throat)",
    y = "Number of Detections"
  )
Chicken_Throat_Unique

group_palette <- c(
  "Chicken_Throat" = "#009E73",
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",
  "Air_Outside" = "#228B22",
  "Cage" = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water" = "#F0E442"
)

summary_long <- summary_df %>%
  rownames_to_column("Group_Name") %>%
  pivot_longer(-c(Group, Group_Name), names_to = "Experiment", values_to = "Count") %>%
  mutate(
    Type = case_when(
      grepl("_Match", Group_Name) ~ "Match",
      grepl("_Unique", Group_Name) ~ "Unique",
      Group_Name == "Chicken_Throat" ~ "Chicken_Throat"
    ),
    Group_Label = gsub("_(Match|Unique)", "", Group_Name),
    Combined_Label = paste(Group_Label, Type, sep = "_"),
    Group_Label = factor(Group_Label, levels = names(group_palette)),
    Combined_Label = factor(Combined_Label, levels = c(
      "Chicken_Throat_Chicken_Throat",
      "Air_Slaughter_Match", "Air_Slaughter_Unique",
      "Air_Holding_Match", "Air_Holding_Unique",
      "Air_Outside_Match", "Air_Outside_Unique",
      "Wash_Water_Match", "Wash_Water_Unique",
      "Drinking_Water_Match", "Drinking_Water_Unique",
      "Cage_Match", "Cage_Unique"
    ))
  )

summary_long <- summary_long %>%
  mutate(
    fill_color = group_palette[Group_Label],
    pattern = ifelse(Type == "Unique", "stripe", "none"),
    pattern_type = ifelse(Type == "Unique", "Unique", "Match")  
  )

Kruskall_Chicken_Throat_Unique <- kwAllPairsDunnTest(summary_long$Count, summary_long$Combined_Label, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chicken_Throat_Unique$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Chicken_Throat_Unique.xlsx", rowNames = TRUE)

Chicken_Throat_Boxplot <- ggplot(summary_long, aes(x = Combined_Label, y = Count)) +
  geom_boxplot_pattern(
    aes(fill = Group_Label, pattern = pattern, pattern_type = pattern_type),
    outlier.shape = NA,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.15,
    pattern_spacing = 0.02,
    color = "black",
    width = 0.6
  ) +
  geom_jitter(
    aes(color = Group_Label),
    width = 0.2,
    size = 2,
    alpha = 0.6
  ) +
  scale_fill_manual(
    name = "Group",
    values = group_palette
  ) +
  scale_color_manual(
    values = group_palette,
    guide = "none"
  ) +
  scale_pattern_identity() +
  scale_pattern_type_manual(
    values = c("Match" = "none", "Unique" = "stripe"),
    guide = guide_legend(
      title = "Detection Type",
      override.aes = list(
        fill = c("grey80", "grey80"),
        pattern = c("none", "stripe"),
        pattern_fill = c(NA, "black"),
        pattern_angle = c(0, 45),
        pattern_density = c(0, 0.15),
        pattern_spacing = c(0, 0.02),
        color = "black"
      )
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    title = "Chicken Throat",
    x = NULL,
    y = "Number of Detections"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"), 
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0, color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "right"
  )
Chicken_Throat_Boxplot

summary_df_meta <- summary_df[c(13)]
summary_df <- summary_df[-c(13)]
keep <- c(1, seq(2, nrow(summary_df), by = 2))      
out  <- summary_df[keep, ]                          
num_cols <- sapply(out, is.numeric)                 
for (i in keep[-1]) {                               
  out[rownames(out) == rownames(summary_df)[i], num_cols] <-
    summary_df[i,   num_cols] +
    summary_df[i+1, num_cols]
  rownames(out)[rownames(out) == rownames(summary_df)[i]] <-
    sub("_(Match|Unique)$", "", rownames(summary_df)[i])
}

out_long <- out %>%                               
  rownames_to_column("Group_Label") %>%
  mutate(
    Group_Label = factor(Group_Label,             
                         levels = names(group_palette))
  ) %>%
  pivot_longer(
    -Group_Label,
    names_to  = "Experiment",
    values_to = "Count"
  )

Kruskall_Chicken_Throat_Total <- kwAllPairsDunnTest(out_long$Count, out_long$Group_Label, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chicken_Throat_Total$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Chicken_Throat_Total.xlsx", rowNames = TRUE)

Total_Chicken_Throat_Boxplot <- ggplot(out_long,
                                        aes(x = Group_Label, y = Count, fill = Group_Label)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(aes(color = Group_Label),
              width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(name = "Group", values = group_palette) +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(title = "Chicken_Throat", x = NULL,
       y = "Number of Poultry Viruses") +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color = "black"),
    axis.text.y      = element_text(face = "bold", color = "black"),
    axis.title.y     = element_text(face = "bold", color = "black"),
    plot.title       = element_text(face = "bold", hjust = 0, color = "black"),
    axis.line.x      = element_blank(),
    axis.line.y      = element_blank(),
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title     = element_text(face = "bold", color = "black"),
    legend.text      = element_text(color = "black"),
    legend.position  = "right"
  )
Total_Chicken_Throat_Boxplot

out_long_filtered <- out %>%                               
  rownames_to_column("Group_Label") %>%
  pivot_longer(
    -Group_Label,                           
    names_to  = "Experiment",
    values_to = "Detections"
  ) %>%
  mutate(
    Experiment_Number = as.numeric(sub("Experiment_", "", Experiment))
  ) %>%
  inner_join(                                
    collected_groups,                        
    by = c("Experiment"      = "Experiment",
           "Group_Label"     = "Location_Grouping")
  )

desired_order <- c(
  "Chicken_Throat",
  "Air_Slaughter",
  "Air_Holding",
  "Air_Outside",
  "Wash_Water",
  "Drinking_Water",
  "Cage"
)

out_long_filtered$Group_Label <- factor(out_long_filtered$Group_Label,
                                        levels = desired_order)

Total_Chicken_Throat_Barplot <- ggplot(
  out_long_filtered,
  aes(x = Group_Label, y = Detections, fill = Group_Label)
) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(
    ~ factor(
      Experiment,
      levels = paste0("Experiment_", sort(unique(out_long_filtered$Experiment_Number)))
    ),
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y     = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x     = element_blank(),
    axis.line.y      = element_line(color = "black", linewidth = 0.8),
    axis.line.x      = element_blank(),
    strip.text       = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title       = element_text(face = "bold", size = 13,
                                    hjust = 0.5, color = "black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  ) +
  labs(
    title = "Detected Viruses (Collapsed Groups)",
    y     = "Number of Detections"
  )
Total_Chicken_Throat_Barplot

#================================================
#### Chicken Cloacal Summary - Barplot/Boxplot
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Chicken_Cloacal_Samples <- rownames(Metadata[Metadata$Location_Grouping == "Chicken_Cloacal", ])
Groups <- list(
  Air_Slaughter = rownames(Metadata[Metadata$Location_Grouping == "Air_Slaughter", ]),
  Air_Holding = rownames(Metadata[Metadata$Location_Grouping == "Air_Holding", ]),
  Air_Outside = rownames(Metadata[Metadata$Location_Grouping == "Air_Outside", ]),
  Cage = rownames(Metadata[Metadata$Location_Grouping == "Cage", ]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping == "Drinking_Water", ]),
  Wash_Water = rownames(Metadata[Metadata$Location_Grouping == "Wash_Water", ])
)

Chicken_Cloacal_Species <- Species[Chicken_Cloacal_Samples, ]
Chicken_Cloacal_Species <- Chicken_Cloacal_Species[, colSums(Chicken_Cloacal_Species != 0) > 0]
Chicken_Cloacal_Taxa <- colnames(Chicken_Cloacal_Species)

Chicken_Cloacal_Data <- Species[, colnames(Species) %in% Chicken_Cloacal_Taxa]
Chicken_Cloacal_Data <- Chicken_Cloacal_Data[!rownames(Chicken_Cloacal_Data) %in% unlist(Groups), ]

Metadata_Filtered <- Metadata[rownames(Chicken_Cloacal_Data), ]
Binary_Data <- (Chicken_Cloacal_Data > 0) * 1

all_experiments <- unique(Metadata$Experiment)
summary_list <- list()

for (exp_id in all_experiments) {
  exp_samples <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  CC_samples <- intersect(exp_samples, Chicken_Cloacal_Samples)
  CC_block <- Binary_Data[rownames(Binary_Data) %in% CC_samples, , drop = FALSE]
  if (nrow(CC_block) == 0) next
  CC_consensus <- colSums(CC_block) > 0
  
  cc_total <- sum(CC_consensus)
  results <- c(Chicken_Cloacal = cc_total)
  
  for (group_name in names(Groups)) {
    group_samples <- intersect(exp_samples, Groups[[group_name]])
    group_block <- Species[group_samples, colnames(CC_block), drop = FALSE]
    group_binary <- (group_block > 0) * 1
    group_detected <- colSums(group_binary) > 0
    
    match_count <- sum(CC_consensus & group_detected)
    unique_count <- sum(!CC_consensus & group_detected)
    
    results[paste0(group_name, "_Match")] <- match_count
    results[paste0(group_name, "_Unique")] <- unique_count
  }
  
  summary_list[[as.character(exp_id)]] <- results
}

summary_df <- do.call(cbind, summary_list)
summary_df <- as.data.frame(summary_df)
colnames(summary_df) <- paste0("Experiment_", colnames(summary_df))
summary_df$Group <- rownames(summary_df)
#write.xlsx(summary_df, file = "Chicken_Cloacal_Species_Detected.xlsx", rowNames = TRUE)

df_match <- summary_df[grep("Chicken_Cloacal|_Match$", summary_df$Group), ]
df_unique <- summary_df[grep("_Unique$", summary_df$Group), ]

long_df_match <- df_match %>%
  pivot_longer(-Group, names_to = "Experiment", values_to = "Detections") %>%
  mutate(
    Group_Label = gsub("_Match", "", Group),
    Experiment_Number = as.numeric(gsub("Experiment_", "", Experiment))
  )

long_df_unique <- df_unique %>%
  pivot_longer(-Group, names_to = "Experiment", values_to = "Detections") %>%
  mutate(
    Group_Label = gsub("_Unique", "", Group),
    Experiment_Number = as.numeric(gsub("Experiment_", "", Experiment))
  )

collected_groups <- Metadata %>%
  group_by(Experiment, Location_Grouping) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Experiment = paste0("Experiment_", Experiment))

long_df_match_filtered <- long_df_match %>%
  mutate(Experiment = paste0("Experiment_", Experiment_Number)) %>%
  inner_join(collected_groups, by = c("Experiment" = "Experiment", "Group_Label" = "Location_Grouping"))

long_df_unique_filtered <- long_df_unique %>%
  mutate(Experiment = paste0("Experiment_", Experiment_Number)) %>%
  inner_join(collected_groups, by = c("Experiment" = "Experiment", "Group_Label" = "Location_Grouping"))

desired_order <- c(
  "Chicken_Cloacal",
  "Air_Slaughter",
  "Air_Holding",
  "Air_Outside",
  "Wash_Water",
  "Drinking_Water",
  "Cage"
)

long_df_match_filtered$Group_Label <- factor(long_df_match_filtered$Group_Label, levels = desired_order)
long_df_unique_filtered$Group_Label <- factor(long_df_unique_filtered$Group_Label, levels = desired_order[-1])

group_palette <- c(
  "Chicken_Cloacal" = "#800000",
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",
  "Air_Outside" = "#228B22",
  "Cage" = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water" = "#F0E442"
)

Chicken_Cloacal_Detected <- ggplot(long_df_match_filtered, aes(x = Group_Label, y = Detections, fill = Group_Label)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_match_filtered$Experiment_Number)))), 
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Detected Viruses (Chicken_Cloacal + Matches)",
    y = "Number of Detections"
  )
Chicken_Cloacal_Detected

Chicken_Cloacal_Unique <- ggplot(long_df_unique_filtered, aes(x = Group_Label, y = Detections, fill = Group_Label)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_unique_filtered$Experiment_Number)))), 
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Unique Detections (Chicken_Cloacal)",
    y = "Number of Detections"
  )
Chicken_Cloacal_Unique

group_palette <- c(
  "Chicken_Cloacal" = "#800000",
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",
  "Air_Outside" = "#228B22",
  "Cage" = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water" = "#F0E442"
)

summary_long <- summary_df %>%
  rownames_to_column("Group_Name") %>%
  pivot_longer(-c(Group, Group_Name), names_to = "Experiment", values_to = "Count") %>%
  mutate(
    Type = case_when(
      grepl("_Match", Group_Name) ~ "Match",
      grepl("_Unique", Group_Name) ~ "Unique",
      Group_Name == "Chicken_Cloacal" ~ "Chicken_Cloacal"
    ),
    Group_Label = gsub("_(Match|Unique)", "", Group_Name),
    Combined_Label = paste(Group_Label, Type, sep = "_"),
    Group_Label = factor(Group_Label, levels = names(group_palette)),
    Combined_Label = factor(Combined_Label, levels = c(
      "Chicken_Cloacal_Chicken_Cloacal",
      "Air_Slaughter_Match", "Air_Slaughter_Unique",
      "Air_Holding_Match", "Air_Holding_Unique",
      "Air_Outside_Match", "Air_Outside_Unique",
      "Wash_Water_Match", "Wash_Water_Unique",
      "Drinking_Water_Match", "Drinking_Water_Unique",
      "Cage_Match", "Cage_Unique"
    ))
  )

summary_long <- summary_long %>%
  mutate(
    fill_color = group_palette[Group_Label],
    pattern = ifelse(Type == "Unique", "stripe", "none"),
    pattern_type = ifelse(Type == "Unique", "Unique", "Match")
  )

Kruskall_Chicken_Cloacal_Unique <- kwAllPairsDunnTest(summary_long$Count, summary_long$Combined_Label, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chicken_Cloacal_Unique$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Chicken_Cloacal_Unique.xlsx", rowNames = TRUE)

Chicken_Cloacal_Boxplot <- ggplot(summary_long, aes(x = Combined_Label, y = Count)) +
  geom_boxplot_pattern(
    aes(fill = Group_Label, pattern = pattern, pattern_type = pattern_type),
    outlier.shape = NA,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.15,
    pattern_spacing = 0.02,
    color = "black",
    width = 0.6
  ) +
  geom_jitter(
    aes(color = Group_Label),
    width = 0.2,
    size = 2,
    alpha = 0.6
  ) +
  scale_fill_manual(
    name = "Group",
    values = group_palette
  ) +
  scale_color_manual(
    values = group_palette,
    guide = "none"
  ) +
  scale_pattern_identity() +
  scale_pattern_type_manual(
    values = c("Match" = "none", "Unique" = "stripe"),
    guide = guide_legend(
      title = "Detection Type",
      override.aes = list(
        fill = c("grey80", "grey80"),
        pattern = c("none", "stripe"),
        pattern_fill = c(NA, "black"),
        pattern_angle = c(0, 45),
        pattern_density = c(0, 0.15),
        pattern_spacing = c(0, 0.02),
        color = "black"
      )
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    title = "Chicken Cloacal",
    x = NULL,
    y = "Number of Detections"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0, color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "right"
  )
Chicken_Cloacal_Boxplot

summary_df_meta <- summary_df[c(13)]
summary_df <- summary_df[-c(13)]
keep <- c(1, seq(2, nrow(summary_df), by = 2))      
out  <- summary_df[keep, ]                          
num_cols <- sapply(out, is.numeric)                 
for (i in keep[-1]) {                               
  out[rownames(out) == rownames(summary_df)[i], num_cols] <-
    summary_df[i,   num_cols] +
    summary_df[i+1, num_cols]
  rownames(out)[rownames(out) == rownames(summary_df)[i]] <-
    sub("_(Match|Unique)$", "", rownames(summary_df)[i])
}

out_long <- out %>%                               
  rownames_to_column("Group_Label") %>%
  mutate(
    Group_Label = factor(Group_Label,             
                         levels = names(group_palette))
  ) %>%
  pivot_longer(
    -Group_Label,
    names_to  = "Experiment",
    values_to = "Count"
  )

Kruskall_Chicken_Cloacal_Total <- kwAllPairsDunnTest(out_long$Count, out_long$Group_Label, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chicken_Cloacal_Total$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Chicken_Cloacal_Total.xlsx", rowNames = TRUE)

Total_Chicken_Cloacal_Boxplot <- ggplot(out_long,
                        aes(x = Group_Label, y = Count, fill = Group_Label)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(aes(color = Group_Label),
              width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(name = "Group", values = group_palette) +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(title = "Chicken_Cloacal", x = NULL,
       y = "Number of Poultry Viruses") +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color = "black"),
    axis.text.y      = element_text(face = "bold", color = "black"),
    axis.title.y     = element_text(face = "bold", color = "black"),
    plot.title       = element_text(face = "bold", hjust = 0, color = "black"),
    axis.line.x      = element_blank(),
    axis.line.y      = element_blank(),
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title     = element_text(face = "bold", color = "black"),
    legend.text      = element_text(color = "black"),
    legend.position  = "right"
  )
Total_Chicken_Cloacal_Boxplot

out_long_filtered <- out %>%                               
  rownames_to_column("Group_Label") %>%
  pivot_longer(
    -Group_Label,                           
    names_to  = "Experiment",
    values_to = "Detections"
  ) %>%
  mutate(
    Experiment_Number = as.numeric(sub("Experiment_", "", Experiment))
  ) %>%
  inner_join(                                
    collected_groups,                        
    by = c("Experiment"      = "Experiment",
           "Group_Label"     = "Location_Grouping")
  )

desired_order <- c(
  "Chicken_Cloacal",
  "Air_Slaughter",
  "Air_Holding",
  "Air_Outside",
  "Wash_Water",
  "Drinking_Water",
  "Cage"
)

out_long_filtered$Group_Label <- factor(out_long_filtered$Group_Label,
                                        levels = desired_order)

Total_Chicken_Cloacal_Barplot <- ggplot(
  out_long_filtered,
  aes(x = Group_Label, y = Detections, fill = Group_Label)
) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(
    ~ factor(
      Experiment,
      levels = paste0("Experiment_", sort(unique(out_long_filtered$Experiment_Number)))
    ),
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y     = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x     = element_blank(),
    axis.line.y      = element_line(color = "black", linewidth = 0.8),
    axis.line.x      = element_blank(),
    strip.text       = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title       = element_text(face = "bold", size = 13,
                                    hjust = 0.5, color = "black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  ) +
  labs(
    title = "Detected Viruses (Collapsed Groups)",
    y     = "Number of Detections"
  )
Total_Chicken_Cloacal_Barplot

#================================================
#### Duck Throat Summary - Barplot/Boxplot
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Duck_Throat_Samples <- rownames(Metadata[Metadata$Location_Grouping == "Duck_Throat", ])

Groups <- list(
  Air_Slaughter = rownames(Metadata[Metadata$Location_Grouping == "Air_Slaughter", ]),
  Air_Holding = rownames(Metadata[Metadata$Location_Grouping == "Air_Holding", ]),
  Air_Outside = rownames(Metadata[Metadata$Location_Grouping == "Air_Outside", ]),
  Cage = rownames(Metadata[Metadata$Location_Grouping == "Cage", ]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping == "Drinking_Water", ]),
  Wash_Water = rownames(Metadata[Metadata$Location_Grouping == "Wash_Water", ])
)

Duck_Throat_Species <- Species[Duck_Throat_Samples, ]
Duck_Throat_Species <- Duck_Throat_Species[, colSums(Duck_Throat_Species != 0) > 0]
Duck_Throat_Taxa <- colnames(Duck_Throat_Species)

Duck_Throat_Data <- Species[, colnames(Species) %in% Duck_Throat_Taxa]
Duck_Throat_Data <- Duck_Throat_Data[!rownames(Duck_Throat_Data) %in% unlist(Groups), ]

Metadata_Filtered <- Metadata[rownames(Duck_Throat_Data), ]
Binary_Data <- (Duck_Throat_Data > 0) * 1

all_experiments <- unique(Metadata$Experiment)
summary_list <- list()

for (exp_id in all_experiments) {
  exp_samples <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  DT_samples <- intersect(exp_samples, Duck_Throat_Samples)
  DT_block <- Binary_Data[rownames(Binary_Data) %in% DT_samples, , drop = FALSE]
  if (nrow(DT_block) == 0) next
  DT_consensus <- colSums(DT_block) > 0
  
  dt_total <- sum(DT_consensus)
  results <- c(Duck_Throat = dt_total)
  
  for (group_name in names(Groups)) {
    group_samples <- intersect(exp_samples, Groups[[group_name]])
    group_block <- Species[group_samples, colnames(DT_block), drop = FALSE]
    group_binary <- (group_block > 0) * 1
    group_detected <- colSums(group_binary) > 0
    
    match_count <- sum(DT_consensus & group_detected)
    unique_count <- sum(!DT_consensus & group_detected)
    
    results[paste0(group_name, "_Match")] <- match_count
    results[paste0(group_name, "_Unique")] <- unique_count
  }
  
  summary_list[[as.character(exp_id)]] <- results
}

summary_df <- do.call(cbind, summary_list)
summary_df <- as.data.frame(summary_df)
colnames(summary_df) <- paste0("Experiment_", colnames(summary_df))
summary_df$Group <- rownames(summary_df)
#write.xlsx(summary_df, file = "Duck_Throat_Species_Detected.xlsx", rowNames = TRUE)

df_match <- summary_df[grep("Duck_Throat|_Match$", summary_df$Group), ]
df_unique <- summary_df[grep("_Unique$", summary_df$Group), ]

long_df_match <- df_match %>%
  pivot_longer(-Group, names_to = "Experiment", values_to = "Detections") %>%
  mutate(
    Group_Label = gsub("_Match", "", Group),
    Experiment_Number = as.numeric(gsub("Experiment_", "", Experiment))
  )

long_df_unique <- df_unique %>%
  pivot_longer(-Group, names_to = "Experiment", values_to = "Detections") %>%
  mutate(
    Group_Label = gsub("_Unique", "", Group),
    Experiment_Number = as.numeric(gsub("Experiment_", "", Experiment))
  )

collected_groups <- Metadata %>%
  group_by(Experiment, Location_Grouping) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Experiment = paste0("Experiment_", Experiment))

long_df_match_filtered <- long_df_match %>%
  mutate(Experiment = paste0("Experiment_", Experiment_Number)) %>%
  inner_join(collected_groups, by = c("Experiment" = "Experiment", "Group_Label" = "Location_Grouping"))

long_df_unique_filtered <- long_df_unique %>%
  mutate(Experiment = paste0("Experiment_", Experiment_Number)) %>%
  inner_join(collected_groups, by = c("Experiment" = "Experiment", "Group_Label" = "Location_Grouping"))

desired_order <- c(
  "Duck_Throat",
  "Air_Slaughter",
  "Air_Holding",
  "Air_Outside",
  "Wash_Water",
  "Drinking_Water",
  "Cage"
)

long_df_match_filtered$Group_Label <- factor(long_df_match_filtered$Group_Label, levels = desired_order)
long_df_unique_filtered$Group_Label <- factor(long_df_unique_filtered$Group_Label, levels = desired_order[-1])

group_palette <- c(
  "Duck_Throat" = "#0072B2",
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",
  "Air_Outside" = "#228B22",
  "Cage" = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water" = "#F0E442"
)

Duck_Throat_Detected <- ggplot(long_df_match_filtered, aes(x = Group_Label, y = Detections, fill = Group_Label)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_match_filtered$Experiment_Number)))), 
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Detected Viruses (Duck_Throat + Matches)",
    y = "Number of Detections"
  )
Duck_Throat_Detected

Duck_Throat_Unique <- ggplot(long_df_unique_filtered, aes(x = Group_Label, y = Detections, fill = Group_Label)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_unique_filtered$Experiment_Number)))), 
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Unique Detections (Duck_Throat)",
    y = "Number of Detections"
  )
Duck_Throat_Unique

summary_long <- summary_df %>%
  rownames_to_column("Group_Name") %>%
  pivot_longer(-c(Group, Group_Name), names_to = "Experiment", values_to = "Count") %>%
  mutate(
    Type = case_when(
      grepl("_Match", Group_Name) ~ "Match",
      grepl("_Unique", Group_Name) ~ "Unique",
      Group_Name == "Duck_Throat" ~ "Duck_Throat"
    ),
    Group_Label = gsub("_(Match|Unique)", "", Group_Name),
    Combined_Label = paste(Group_Label, Type, sep = "_"),
    Group_Label = factor(Group_Label, levels = names(group_palette)),
    Combined_Label = factor(Combined_Label, levels = c(
      "Duck_Throat_Duck_Throat",
      "Air_Slaughter_Match", "Air_Slaughter_Unique",
      "Air_Holding_Match", "Air_Holding_Unique",
      "Air_Outside_Match", "Air_Outside_Unique",
      "Wash_Water_Match", "Wash_Water_Unique",
      "Drinking_Water_Match", "Drinking_Water_Unique",
      "Cage_Match", "Cage_Unique"
    ))
  ) %>%
  mutate(
    fill_color = group_palette[Group_Label],
    pattern = ifelse(Type == "Unique", "stripe", "none"),
    pattern_type = ifelse(Type == "Unique", "Unique", "Match")
  )

Kruskall_Duck_Throat_Unique <- kwAllPairsDunnTest(summary_long$Count, summary_long$Combined_Label, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Duck_Throat_Unique$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Duck_Throat_Unique.xlsx", rowNames = TRUE)

Duck_Throat_Boxplot <- ggplot(summary_long, aes(x = Combined_Label, y = Count)) +
  geom_boxplot_pattern(
    aes(fill = Group_Label, pattern = pattern, pattern_type = pattern_type),
    outlier.shape = NA,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.15,
    pattern_spacing = 0.02,
    color = "black",
    width = 0.6
  ) +
  geom_jitter(
    aes(color = Group_Label),
    width = 0.2,
    size = 2,
    alpha = 0.6
  ) +
  scale_fill_manual(
    name = "Group",
    values = group_palette
  ) +
  scale_color_manual(
    values = group_palette,
    guide = "none"
  ) +
  scale_pattern_identity() +
  scale_pattern_type_manual(
    values = c("Match" = "none", "Unique" = "stripe"),
    guide = guide_legend(
      title = "Detection Type",
      override.aes = list(
        fill = c("grey80", "grey80"),
        pattern = c("none", "stripe"),
        pattern_fill = c(NA, "black"),
        pattern_angle = c(0, 45),
        pattern_density = c(0, 0.15),
        pattern_spacing = c(0, 0.02),
        color = "black"
      )
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Duck Throat",
    x = NULL,
    y = "Number of Detections"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"), 
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0, color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "right"
  )

Duck_Throat_Detected
Duck_Throat_Unique
Duck_Throat_Boxplot

summary_df_meta <- summary_df[c(13)]
summary_df <- summary_df[-c(13)]
keep <- c(1, seq(2, nrow(summary_df), by = 2))      
out  <- summary_df[keep, ]                          
num_cols <- sapply(out, is.numeric)                 
for (i in keep[-1]) {                               
  out[rownames(out) == rownames(summary_df)[i], num_cols] <-
    summary_df[i,   num_cols] +
    summary_df[i+1, num_cols]
  rownames(out)[rownames(out) == rownames(summary_df)[i]] <-
    sub("_(Match|Unique)$", "", rownames(summary_df)[i])
}

out_long <- out %>%                               
  rownames_to_column("Group_Label") %>%
  mutate(
    Group_Label = factor(Group_Label,             
                         levels = names(group_palette))
  ) %>%
  pivot_longer(
    -Group_Label,
    names_to  = "Experiment",
    values_to = "Count"
  )

Kruskall_Duck_Throat_Total <- kwAllPairsDunnTest(out_long$Count, out_long$Group_Label, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Duck_Throat_Total$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Duck_Throat_Total.xlsx", rowNames = TRUE)

Total_Duck_Throat_Boxplot <- ggplot(out_long,
                                       aes(x = Group_Label, y = Count, fill = Group_Label)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(aes(color = Group_Label),
              width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(name = "Group", values = group_palette) +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(title = "Duck_Throat", x = NULL,
       y = "Number of Poultry Viruses") +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color = "black"),
    axis.text.y      = element_text(face = "bold", color = "black"),
    axis.title.y     = element_text(face = "bold", color = "black"),
    plot.title       = element_text(face = "bold", hjust = 0, color = "black"),
    axis.line.x      = element_blank(),
    axis.line.y      = element_blank(),
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title     = element_text(face = "bold", color = "black"),
    legend.text      = element_text(color = "black"),
    legend.position  = "right"
  )
Total_Duck_Throat_Boxplot

out_long_filtered <- out %>%                               
  rownames_to_column("Group_Label") %>%
  pivot_longer(
    -Group_Label,                           
    names_to  = "Experiment",
    values_to = "Detections"
  ) %>%
  mutate(
    Experiment_Number = as.numeric(sub("Experiment_", "", Experiment))
  ) %>%
  inner_join(                                
    collected_groups,                        
    by = c("Experiment"      = "Experiment",
           "Group_Label"     = "Location_Grouping")
  )

desired_order <- c(
  "Duck_Throat",
  "Air_Slaughter",
  "Air_Holding",
  "Air_Outside",
  "Wash_Water",
  "Drinking_Water",
  "Cage"
)

out_long_filtered$Group_Label <- factor(out_long_filtered$Group_Label,
                                        levels = desired_order)

Total_Duck_Throat_Barplot <- ggplot(
  out_long_filtered,
  aes(x = Group_Label, y = Detections, fill = Group_Label)
) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(
    ~ factor(
      Experiment,
      levels = paste0("Experiment_", sort(unique(out_long_filtered$Experiment_Number)))
    ),
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y     = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x     = element_blank(),
    axis.line.y      = element_line(color = "black", linewidth = 0.8),
    axis.line.x      = element_blank(),
    strip.text       = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title       = element_text(face = "bold", size = 13,
                                    hjust = 0.5, color = "black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  ) +
  labs(
    title = "Detected Viruses (Collapsed Groups)",
    y     = "Number of Detections"
  )
Total_Duck_Throat_Barplot

#================================================
#### Duck Cloacal Summary - Barplot/Boxplot
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Duck_Cloacal_Samples <- rownames(Metadata[Metadata$Location_Grouping == "Duck_Cloacal", ])

Groups <- list(
  Air_Slaughter = rownames(Metadata[Metadata$Location_Grouping == "Air_Slaughter", ]),
  Air_Holding = rownames(Metadata[Metadata$Location_Grouping == "Air_Holding", ]),
  Air_Outside = rownames(Metadata[Metadata$Location_Grouping == "Air_Outside", ]),
  Cage = rownames(Metadata[Metadata$Location_Grouping == "Cage", ]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping == "Drinking_Water", ]),
  Wash_Water = rownames(Metadata[Metadata$Location_Grouping == "Wash_Water", ])
)

Duck_Cloacal_Species <- Species[Duck_Cloacal_Samples, ]
Duck_Cloacal_Species <- Duck_Cloacal_Species[, colSums(Duck_Cloacal_Species != 0) > 0]
Duck_Cloacal_Taxa <- colnames(Duck_Cloacal_Species)

Duck_Cloacal_Data <- Species[, colnames(Species) %in% Duck_Cloacal_Taxa]
Duck_Cloacal_Data <- Duck_Cloacal_Data[!rownames(Duck_Cloacal_Data) %in% unlist(Groups), ]

Metadata_Filtered <- Metadata[rownames(Duck_Cloacal_Data), ]
Binary_Data <- (Duck_Cloacal_Data > 0) * 1

all_experiments <- unique(Metadata$Experiment)
summary_list <- list()

for (exp_id in all_experiments) {
  exp_samples <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  DC_samples <- intersect(exp_samples, Duck_Cloacal_Samples)
  DC_block <- Binary_Data[rownames(Binary_Data) %in% DC_samples, , drop = FALSE]
  if (nrow(DC_block) == 0) next
  DC_consensus <- colSums(DC_block) > 0
  
  dc_total <- sum(DC_consensus)
  results <- c(Duck_Cloacal = dc_total)
  
  for (group_name in names(Groups)) {
    group_samples <- intersect(exp_samples, Groups[[group_name]])
    group_block <- Species[group_samples, colnames(DC_block), drop = FALSE]
    group_binary <- (group_block > 0) * 1
    group_detected <- colSums(group_binary) > 0
    
    match_count <- sum(DC_consensus & group_detected)
    unique_count <- sum(!DC_consensus & group_detected)
    
    results[paste0(group_name, "_Match")] <- match_count
    results[paste0(group_name, "_Unique")] <- unique_count
  }
  
  summary_list[[as.character(exp_id)]] <- results
}

summary_df <- do.call(cbind, summary_list)
summary_df <- as.data.frame(summary_df)
colnames(summary_df) <- paste0("Experiment_", colnames(summary_df))
summary_df$Group <- rownames(summary_df)
#write.xlsx(summary_df, file = "Duck_Cloacal_Species_Detected.xlsx", rowNames = TRUE)

df_match <- summary_df[grep("Duck_Cloacal|_Match$", summary_df$Group), ]
df_unique <- summary_df[grep("_Unique$", summary_df$Group), ]

long_df_match <- df_match %>%
  pivot_longer(-Group, names_to = "Experiment", values_to = "Detections") %>%
  mutate(
    Group_Label = gsub("_Match", "", Group),
    Experiment_Number = as.numeric(gsub("Experiment_", "", Experiment))
  )

long_df_unique <- df_unique %>%
  pivot_longer(-Group, names_to = "Experiment", values_to = "Detections") %>%
  mutate(
    Group_Label = gsub("_Unique", "", Group),
    Experiment_Number = as.numeric(gsub("Experiment_", "", Experiment))
  )

collected_groups <- Metadata %>%
  group_by(Experiment, Location_Grouping) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Experiment = paste0("Experiment_", Experiment))

long_df_match_filtered <- long_df_match %>%
  mutate(Experiment = paste0("Experiment_", Experiment_Number)) %>%
  inner_join(collected_groups, by = c("Experiment" = "Experiment", "Group_Label" = "Location_Grouping"))

long_df_unique_filtered <- long_df_unique %>%
  mutate(Experiment = paste0("Experiment_", Experiment_Number)) %>%
  inner_join(collected_groups, by = c("Experiment" = "Experiment", "Group_Label" = "Location_Grouping"))

desired_order <- c(
  "Duck_Cloacal",
  "Air_Slaughter",
  "Air_Holding",
  "Air_Outside",
  "Wash_Water",
  "Drinking_Water",
  "Cage"
)

group_palette <- c(
  "Duck_Cloacal" = "#964B00",
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",
  "Air_Outside" = "#228B22",
  "Cage" = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water" = "#F0E442"
)

long_df_match_filtered$Group_Label <- factor(long_df_match_filtered$Group_Label, levels = desired_order)
long_df_unique_filtered$Group_Label <- factor(long_df_unique_filtered$Group_Label, levels = desired_order[-1])

Duck_Cloacal_Detected <- ggplot(long_df_match_filtered, aes(x = Group_Label, y = Detections, fill = Group_Label)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_match_filtered$Experiment_Number)))), 
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Detected Viruses (Duck_Cloacal + Matches)",
    y = "Number of Detections"
  )

Duck_Cloacal_Unique <- ggplot(long_df_unique_filtered, aes(x = Group_Label, y = Detections, fill = Group_Label)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_unique_filtered$Experiment_Number)))), 
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Unique Detections (Duck_Cloacal)",
    y = "Number of Detections"
  )

summary_long <- summary_df %>%
  rownames_to_column("Group_Name") %>%
  pivot_longer(-c(Group, Group_Name), names_to = "Experiment", values_to = "Count") %>%
  mutate(
    Type = case_when(
      grepl("_Match", Group_Name) ~ "Match",
      grepl("_Unique", Group_Name) ~ "Unique",
      Group_Name == "Duck_Cloacal" ~ "Duck_Cloacal"
    ),
    Group_Label = gsub("_(Match|Unique)", "", Group_Name),
    Combined_Label = paste(Group_Label, Type, sep = "_"),
    Group_Label = factor(Group_Label, levels = names(group_palette)),
    Combined_Label = factor(Combined_Label, levels = c(
      "Duck_Cloacal_Duck_Cloacal",
      "Air_Slaughter_Match", "Air_Slaughter_Unique",
      "Air_Holding_Match", "Air_Holding_Unique",
      "Air_Outside_Match", "Air_Outside_Unique",
      "Wash_Water_Match", "Wash_Water_Unique",
      "Drinking_Water_Match", "Drinking_Water_Unique",
      "Cage_Match", "Cage_Unique"
    ))
  ) %>%
  mutate(
    fill_color = group_palette[Group_Label],
    pattern = ifelse(Type == "Unique", "stripe", "none"),
    pattern_type = ifelse(Type == "Unique", "Unique", "Match")
  )

Kruskall_Duck_Cloacal_Unique <- kwAllPairsDunnTest(summary_long$Count, summary_long$Combined_Label, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Duck_Cloacal_Unique$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Duck_Cloacal_Unique.xlsx", rowNames = TRUE)

Duck_Cloacal_Boxplot <- ggplot(summary_long, aes(x = Combined_Label, y = Count)) +
  geom_boxplot_pattern(
    aes(fill = Group_Label, pattern = pattern, pattern_type = pattern_type),
    outlier.shape = NA,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.15,
    pattern_spacing = 0.02,
    color = "black",
    width = 0.6
  ) +
  geom_jitter(
    aes(color = Group_Label),
    width = 0.2,
    size = 2,
    alpha = 0.6
  ) +
  scale_fill_manual(name = "Group", values = group_palette) +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_pattern_identity() +
  scale_pattern_type_manual(
    values = c("Match" = "none", "Unique" = "stripe"),
    guide = guide_legend(
      title = "Detection Type",
      override.aes = list(
        fill = c("grey80", "grey80"),
        pattern = c("none", "stripe"),
        pattern_fill = c(NA, "black"),
        pattern_angle = c(0, 45),
        pattern_density = c(0, 0.15),
        pattern_spacing = c(0, 0.02),
        color = "black"
      )
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    title = "Duck Cloacal",
    x = NULL,
    y = "Number of Detections"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0, color = "black"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "right"
  )

Duck_Cloacal_Detected
Duck_Cloacal_Unique
Duck_Cloacal_Boxplot

summary_df_meta <- summary_df[c(13)]
summary_df <- summary_df[-c(13)]
keep <- c(1, seq(2, nrow(summary_df), by = 2))      
out  <- summary_df[keep, ]                          
num_cols <- sapply(out, is.numeric)                 
for (i in keep[-1]) {                               
  out[rownames(out) == rownames(summary_df)[i], num_cols] <-
    summary_df[i,   num_cols] +
    summary_df[i+1, num_cols]
  rownames(out)[rownames(out) == rownames(summary_df)[i]] <-
    sub("_(Match|Unique)$", "", rownames(summary_df)[i])
}

out_long <- out %>%                               
  rownames_to_column("Group_Label") %>%
  mutate(
    Group_Label = factor(Group_Label,             
                         levels = names(group_palette))
  ) %>%
  pivot_longer(
    -Group_Label,
    names_to  = "Experiment",
    values_to = "Count"
  )

Kruskall_Duck_Cloacal_Total <- kwAllPairsDunnTest(out_long$Count, out_long$Group_Label, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Duck_Cloacal_Total$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Duck_Cloacal_Total.xlsx", rowNames = TRUE)

Total_Duck_Cloacal_Boxplot <- ggplot(out_long,
                                    aes(x = Group_Label, y = Count, fill = Group_Label)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(aes(color = Group_Label),
              width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(name = "Group", values = group_palette) +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(title = "Duck_Cloacal", x = NULL,
       y = "Number of Poultry Viruses") +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color = "black"),
    axis.text.y      = element_text(face = "bold", color = "black"),
    axis.title.y     = element_text(face = "bold", color = "black"),
    plot.title       = element_text(face = "bold", hjust = 0, color = "black"),
    axis.line.x      = element_blank(),
    axis.line.y      = element_blank(),
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title     = element_text(face = "bold", color = "black"),
    legend.text      = element_text(color = "black"),
    legend.position  = "right"
  )
Total_Duck_Cloacal_Boxplot

out_long_filtered <- out %>%                               
  rownames_to_column("Group_Label") %>%
  pivot_longer(
    -Group_Label,                           
    names_to  = "Experiment",
    values_to = "Detections"
  ) %>%
  mutate(
    Experiment_Number = as.numeric(sub("Experiment_", "", Experiment))
  ) %>%
  inner_join(                                
    collected_groups,                        
    by = c("Experiment"      = "Experiment",
           "Group_Label"     = "Location_Grouping")
  )

desired_order <- c(
  "Duck_Cloacal",
  "Air_Slaughter",
  "Air_Holding",
  "Air_Outside",
  "Wash_Water",
  "Drinking_Water",
  "Cage"
)

out_long_filtered$Group_Label <- factor(out_long_filtered$Group_Label,
                                        levels = desired_order)

Total_Duck_Cloacal_Barplot <- ggplot(
  out_long_filtered,
  aes(x = Group_Label, y = Detections, fill = Group_Label)
) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(
    ~ factor(
      Experiment,
      levels = paste0("Experiment_", sort(unique(out_long_filtered$Experiment_Number)))
    ),
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y     = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x     = element_blank(),
    axis.line.y      = element_line(color = "black", linewidth = 0.8),
    axis.line.x      = element_blank(),
    strip.text       = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title       = element_text(face = "bold", size = 13,
                                    hjust = 0.5, color = "black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  ) +
  labs(
    title = "Detected Viruses (Collapsed Groups)",
    y     = "Number of Detections"
  )
Total_Duck_Cloacal_Barplot

Chicken_Boxplots = Chicken_Throat_Boxplot / Chicken_Cloacal_Boxplot
ggsave(filename="Chicken_Unique_Boxplots.pdf", plot = Chicken_Boxplots, width=200, height=280, units="mm")

Duck_Boxplots = Duck_Throat_Boxplot / Duck_Cloacal_Boxplot
ggsave(filename="Duck_Unique_Boxplots.pdf", plot = Duck_Boxplots, width=200, height=280, units="mm")

Total_Chicken_Boxplots = Total_Chicken_Throat_Boxplot | Total_Chicken_Cloacal_Boxplot
Total_Duck_Boxplots = Total_Duck_Throat_Boxplot | Total_Duck_Cloacal_Boxplot
Total_Boxplots = Total_Chicken_Boxplots / Total_Duck_Boxplots
ggsave(filename="Species_Total_Boxplots.pdf", plot = Boxplots, width=400, height=200, units="mm")

Vertical_Boxplots = Total_Chicken_Throat_Boxplot / Total_Chicken_Cloacal_Boxplot / Total_Duck_Throat_Boxplot / Total_Duck_Cloacal_Boxplot
ggsave(filename="Species_Total_Boxplots_Vertical.pdf", plot = Vertical_Boxplots, width=200, height=400, units="mm")

Barplot_Detected = Chicken_Throat_Detected / Chicken_Cloacal_Detected / Duck_Throat_Detected / Duck_Cloacal_Detected
ggsave(filename="Species_Recaptured_Barplots.pdf", plot = Barplot_Detected, width=200, height=280, units="mm")

Barplot_Unique = Chicken_Throat_Unique / Chicken_Cloacal_Unique / Duck_Throat_Unique / Duck_Cloacal_Unique
ggsave(filename="Species_Unique_Environment_Barplots.pdf", plot = Barplot_Unique, width=200, height=280, units="mm")

Barplot_Total = Total_Chicken_Throat_Barplot / Total_Chicken_Cloacal_Barplot / Total_Duck_Throat_Barplot / Total_Duck_Cloacal_Barplot
ggsave(filename="Total_Virus_Barplots.pdf", plot = Barplot_Total, width=200, height=280, units="mm")

#================================================
#### Different Environment Samples Capture - Stacked Barplot - Chicken Throat
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Chicken_Throat_Samples <- rownames(Metadata[Metadata$Location_Grouping == "Chicken_Throat", ])

Groups <- list(
  Air_Slaughter = rownames(Metadata[Metadata$Location_Grouping == "Air_Slaughter", ]),
  Air_Holding = rownames(Metadata[Metadata$Location_Grouping == "Air_Holding", ]),
  Air_Outside = rownames(Metadata[Metadata$Location_Grouping == "Air_Outside", ]),
  Cage = rownames(Metadata[Metadata$Location_Grouping == "Cage", ]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping == "Drinking_Water", ]),
  Wash_Water = rownames(Metadata[Metadata$Location_Grouping == "Wash_Water", ])
)

Chicken_Throat_Species <- Species[Chicken_Throat_Samples, ]
Chicken_Throat_Species <- Chicken_Throat_Species[, colSums(Chicken_Throat_Species != 0) > 0]
Chicken_Throat_Taxa <- colnames(Chicken_Throat_Species)

Chicken_Throat_Data <- Species[, colnames(Species) %in% Chicken_Throat_Taxa]
Chicken_Throat_Data <- Chicken_Throat_Data[!rownames(Chicken_Throat_Data) %in% unlist(Groups), ]

Metadata_Filtered <- Metadata[rownames(Chicken_Throat_Data), ]
Binary_Data <- (Chicken_Throat_Data > 0) * 1
all_experiments <- unique(Metadata$Experiment)
match_only_summary <- list()

for (exp_id in all_experiments) {
  exp_samples <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  CT_samples <- intersect(exp_samples, Chicken_Throat_Samples)
  CT_block <- Binary_Data[rownames(Binary_Data) %in% CT_samples, , drop = FALSE]
  if (nrow(CT_block) == 0) next
  CT_consensus <- colSums(CT_block) > 0
  ct_total <- sum(CT_consensus)
  
  results <- c(Chicken_Throat = ct_total)
  
  virus_match_matrix <- matrix(0, nrow = sum(CT_consensus), ncol = length(Groups))
  colnames(virus_match_matrix) <- names(Groups)
  rownames(virus_match_matrix) <- names(CT_consensus[CT_consensus])
  
  i <- 1
  for (group_name in names(Groups)) {
    group_samples <- intersect(exp_samples, Groups[[group_name]])
    group_block <- Species[group_samples, colnames(CT_block), drop = FALSE]
    group_binary <- (group_block > 0) * 1
    group_detected <- colSums(group_binary) > 0
    
    matched_taxa <- names(group_detected[group_detected & CT_consensus])
    virus_match_matrix[matched_taxa, group_name] <- 1
  }
  
  group_counts <- colSums(virus_match_matrix[rowSums(virus_match_matrix) == 1, , drop = FALSE])
  multiple_count <- sum(rowSums(virus_match_matrix) > 1)
  
  for (group_name in names(Groups)) {
    results[group_name] <- group_counts[group_name]
  }
  
  results["Multiple"] <- multiple_count
  match_only_summary[[as.character(exp_id)]] <- results
}

match_summary_df <- do.call(cbind, match_only_summary)
match_summary_df <- as.data.frame(match_summary_df)
colnames(match_summary_df) <- paste0("Experiment_", colnames(match_summary_df))
match_summary_df$Group <- rownames(match_summary_df)

desired_order <- c("Chicken_Throat", names(Groups), "Multiple")
match_summary_df <- match_summary_df[match_summary_df$Group %in% desired_order, ]
match_summary_df$Group <- factor(match_summary_df$Group, levels = desired_order)
long_df_match <- pivot_longer(match_summary_df, cols = starts_with("Experiment_"),
                              names_to = "Experiment", values_to = "Detections")
long_df_match$Experiment_Number <- as.numeric(gsub("Experiment_", "", long_df_match$Experiment))
long_df_match$Group_Label <- ifelse(long_df_match$Group == "Chicken_Throat", 
                                    "Chicken_Throat", 
                                    "Environmental")
long_df_match$Group <- factor(long_df_match$Group)
long_df_stack <- long_df_match[long_df_match$Group != "Chicken_Throat", ]

group_palette <- c(
  "Chicken_Throat" = "#009E73",
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",
  "Air_Outside" = "#228B22",
  "Cage" = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water" = "#F0E442",
  "Multiple" = "grey"
)

Chicken_Throat_Stacked <- ggplot() +
  geom_bar(
    data = long_df_match[long_df_match$Group == "Chicken_Throat", ],
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  geom_bar(
    data = long_df_stack,
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_match$Experiment_Number)))),
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "Detected Viruses (Chicken_Throat + Matches)",
    y = "Number of Detections"
  )

Chicken_Throat_Stacked

#================================================
#### Different Environment Samples Capture - Stacked Barplot - Chicken Cloacal
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Chicken_Cloacal_Samples <- rownames(Metadata[Metadata$Location_Grouping == "Chicken_Cloacal", ])

Groups <- list(
  Air_Slaughter = rownames(Metadata[Metadata$Location_Grouping == "Air_Slaughter", ]),
  Air_Holding = rownames(Metadata[Metadata$Location_Grouping == "Air_Holding", ]),
  Air_Outside = rownames(Metadata[Metadata$Location_Grouping == "Air_Outside", ]),
  Cage = rownames(Metadata[Metadata$Location_Grouping == "Cage", ]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping == "Drinking_Water", ]),
  Wash_Water = rownames(Metadata[Metadata$Location_Grouping == "Wash_Water", ])
)

Chicken_Cloacal_Species <- Species[Chicken_Cloacal_Samples, ]
Chicken_Cloacal_Species <- Chicken_Cloacal_Species[, colSums(Chicken_Cloacal_Species != 0) > 0]
Chicken_Cloacal_Taxa <- colnames(Chicken_Cloacal_Species)

Chicken_Cloacal_Data <- Species[, colnames(Species) %in% Chicken_Cloacal_Taxa]
Chicken_Cloacal_Data <- Chicken_Cloacal_Data[!rownames(Chicken_Cloacal_Data) %in% unlist(Groups), ]

Metadata_Filtered <- Metadata[rownames(Chicken_Cloacal_Data), ]
Binary_Data <- (Chicken_Cloacal_Data > 0) * 1
all_experiments <- unique(Metadata$Experiment)
match_only_summary <- list()

for (exp_id in all_experiments) {
  exp_samples <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  CC_samples <- intersect(exp_samples, Chicken_Cloacal_Samples)
  CC_block <- Binary_Data[rownames(Binary_Data) %in% CC_samples, , drop = FALSE]
  if (nrow(CC_block) == 0) next
  CC_consensus <- colSums(CC_block) > 0
  cc_total <- sum(CC_consensus)
  
  results <- c(Chicken_Cloacal = cc_total)
  
  virus_match_matrix <- matrix(0, nrow = sum(CC_consensus), ncol = length(Groups))
  colnames(virus_match_matrix) <- names(Groups)
  rownames(virus_match_matrix) <- names(CC_consensus[CC_consensus])
  
  i <- 1
  for (group_name in names(Groups)) {
    group_samples <- intersect(exp_samples, Groups[[group_name]])
    group_block <- Species[group_samples, colnames(CC_block), drop = FALSE]
    group_binary <- (group_block > 0) * 1
    group_detected <- colSums(group_binary) > 0
    
    matched_taxa <- names(group_detected[group_detected & CC_consensus])
    virus_match_matrix[matched_taxa, group_name] <- 1
  }
  
  group_counts <- colSums(virus_match_matrix[rowSums(virus_match_matrix) == 1, , drop = FALSE])
  multiple_count <- sum(rowSums(virus_match_matrix) > 1)
  
  for (group_name in names(Groups)) {
    results[group_name] <- group_counts[group_name]
  }
  
  results["Multiple"] <- multiple_count
  match_only_summary[[as.character(exp_id)]] <- results
}

match_summary_df <- do.call(cbind, match_only_summary)
match_summary_df <- as.data.frame(match_summary_df)
colnames(match_summary_df) <- paste0("Experiment_", colnames(match_summary_df))
match_summary_df$Group <- rownames(match_summary_df)

desired_order <- c("Chicken_Cloacal", names(Groups), "Multiple")
match_summary_df <- match_summary_df[match_summary_df$Group %in% desired_order, ]
match_summary_df$Group <- factor(match_summary_df$Group, levels = desired_order)
long_df_match <- pivot_longer(match_summary_df, cols = starts_with("Experiment_"),
                              names_to = "Experiment", values_to = "Detections")
long_df_match$Experiment_Number <- as.numeric(gsub("Experiment_", "", long_df_match$Experiment))
long_df_match$Group_Label <- ifelse(long_df_match$Group == "Chicken_Cloacal", 
                                    "Chicken_Cloacal", 
                                    "Environmental")
long_df_match$Group <- factor(long_df_match$Group)
long_df_stack <- long_df_match[long_df_match$Group != "Chicken_Cloacal", ]

group_palette <- c(
  "Chicken_Cloacal" = "#800000",
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",
  "Air_Outside" = "#228B22",
  "Cage" = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water" = "#F0E442",
  "Multiple" = "grey"
)

Chicken_Cloacal_Stacked <- ggplot() +
  geom_bar(
    data = long_df_match[long_df_match$Group == "Chicken_Cloacal", ],
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  geom_bar(
    data = long_df_stack,
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_match$Experiment_Number)))),
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "Detected Viruses (Chicken_Cloacal + Matches)",
    y = "Number of Detections"
  )

Chicken_Cloacal_Stacked

#================================================
#### Different Environment Samples Capture - Stacked Barplot - Duck Throat
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Duck_Throat_Samples <- rownames(Metadata[Metadata$Location_Grouping == "Duck_Throat", ])

Groups <- list(
  Air_Slaughter = rownames(Metadata[Metadata$Location_Grouping == "Air_Slaughter", ]),
  Air_Holding = rownames(Metadata[Metadata$Location_Grouping == "Air_Holding", ]),
  Air_Outside = rownames(Metadata[Metadata$Location_Grouping == "Air_Outside", ]),
  Cage = rownames(Metadata[Metadata$Location_Grouping == "Cage", ]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping == "Drinking_Water", ]),
  Wash_Water = rownames(Metadata[Metadata$Location_Grouping == "Wash_Water", ])
)

Duck_Throat_Species <- Species[Duck_Throat_Samples, ]
Duck_Throat_Species <- Duck_Throat_Species[, colSums(Duck_Throat_Species != 0) > 0]
Duck_Throat_Taxa <- colnames(Duck_Throat_Species)

Duck_Throat_Data <- Species[, colnames(Species) %in% Duck_Throat_Taxa]
Duck_Throat_Data <- Duck_Throat_Data[!rownames(Duck_Throat_Data) %in% unlist(Groups), ]

Metadata_Filtered <- Metadata[rownames(Duck_Throat_Data), ]
Binary_Data <- (Duck_Throat_Data > 0) * 1
all_experiments <- unique(Metadata$Experiment)
match_only_summary <- list()

for (exp_id in all_experiments) {
  exp_samples <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  DT_samples <- intersect(exp_samples, Duck_Throat_Samples)
  DT_block <- Binary_Data[rownames(Binary_Data) %in% DT_samples, , drop = FALSE]
  if (nrow(DT_block) == 0) next
  DT_consensus <- colSums(DT_block) > 0
  dt_total <- sum(DT_consensus)
  
  results <- c(Duck_Throat = dt_total)
  
  virus_match_matrix <- matrix(0, nrow = sum(DT_consensus), ncol = length(Groups))
  colnames(virus_match_matrix) <- names(Groups)
  rownames(virus_match_matrix) <- names(DT_consensus[DT_consensus])
  
  i <- 1
  for (group_name in names(Groups)) {
    group_samples <- intersect(exp_samples, Groups[[group_name]])
    group_block <- Species[group_samples, colnames(DT_block), drop = FALSE]
    group_binary <- (group_block > 0) * 1
    group_detected <- colSums(group_binary) > 0
    
    matched_taxa <- names(group_detected[group_detected & DT_consensus])
    virus_match_matrix[matched_taxa, group_name] <- 1
  }
  
  group_counts <- colSums(virus_match_matrix[rowSums(virus_match_matrix) == 1, , drop = FALSE])
  multiple_count <- sum(rowSums(virus_match_matrix) > 1)
  
  for (group_name in names(Groups)) {
    results[group_name] <- group_counts[group_name]
  }
  
  results["Multiple"] <- multiple_count
  match_only_summary[[as.character(exp_id)]] <- results
}

match_summary_df <- do.call(cbind, match_only_summary)
match_summary_df <- as.data.frame(match_summary_df)
colnames(match_summary_df) <- paste0("Experiment_", colnames(match_summary_df))
match_summary_df$Group <- rownames(match_summary_df)

desired_order <- c("Duck_Throat", names(Groups), "Multiple")
match_summary_df <- match_summary_df[match_summary_df$Group %in% desired_order, ]
match_summary_df$Group <- factor(match_summary_df$Group, levels = desired_order)
long_df_match <- pivot_longer(match_summary_df, cols = starts_with("Experiment_"),
                              names_to = "Experiment", values_to = "Detections")
long_df_match$Experiment_Number <- as.numeric(gsub("Experiment_", "", long_df_match$Experiment))
long_df_match$Group_Label <- ifelse(long_df_match$Group == "Duck_Throat", 
                                    "Duck_Throat", 
                                    "Environmental")
long_df_match$Group <- factor(long_df_match$Group)
long_df_stack <- long_df_match[long_df_match$Group != "Duck_Throat", ]

group_palette <- c(
  "Duck_Throat" = "#0072B2",
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",
  "Air_Outside" = "#228B22",
  "Cage" = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water" = "#F0E442",
  "Multiple" = "grey"
)

Duck_Throat_Stacked <- ggplot() +
  geom_bar(
    data = long_df_match[long_df_match$Group == "Duck_Throat", ],
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  geom_bar(
    data = long_df_stack,
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_match$Experiment_Number)))),
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "Detected Viruses (Duck_Throat + Matches)",
    y = "Number of Detections"
  )

Duck_Throat_Stacked

#================================================
#### Different Environment Samples Capture - Stacked Barplot - Duck Cloacal
#================================================
Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Duck_Cloacal_Samples <- rownames(Metadata[Metadata$Location_Grouping == "Duck_Cloacal", ])

Groups <- list(
  Air_Slaughter = rownames(Metadata[Metadata$Location_Grouping == "Air_Slaughter", ]),
  Air_Holding = rownames(Metadata[Metadata$Location_Grouping == "Air_Holding", ]),
  Air_Outside = rownames(Metadata[Metadata$Location_Grouping == "Air_Outside", ]),
  Cage = rownames(Metadata[Metadata$Location_Grouping == "Cage", ]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping == "Drinking_Water", ]),
  Wash_Water = rownames(Metadata[Metadata$Location_Grouping == "Wash_Water", ])
)

Duck_Cloacal_Species <- Species[Duck_Cloacal_Samples, ]
Duck_Cloacal_Species <- Duck_Cloacal_Species[, colSums(Duck_Cloacal_Species != 0) > 0]
Duck_Cloacal_Taxa <- colnames(Duck_Cloacal_Species)

Duck_Cloacal_Data <- Species[, colnames(Species) %in% Duck_Cloacal_Taxa]
Duck_Cloacal_Data <- Duck_Cloacal_Data[!rownames(Duck_Cloacal_Data) %in% unlist(Groups), ]

Metadata_Filtered <- Metadata[rownames(Duck_Cloacal_Data), ]
Binary_Data <- (Duck_Cloacal_Data > 0) * 1
all_experiments <- unique(Metadata$Experiment)
match_only_summary <- list()

for (exp_id in all_experiments) {
  exp_samples <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  DC_samples <- intersect(exp_samples, Duck_Cloacal_Samples)
  DC_block <- Binary_Data[rownames(Binary_Data) %in% DC_samples, , drop = FALSE]
  if (nrow(DC_block) == 0) next
  DC_consensus <- colSums(DC_block) > 0
  dc_total <- sum(DC_consensus)
  
  results <- c(Duck_Cloacal = dc_total)
  
  virus_match_matrix <- matrix(0, nrow = sum(DC_consensus), ncol = length(Groups))
  colnames(virus_match_matrix) <- names(Groups)
  rownames(virus_match_matrix) <- names(DC_consensus[DC_consensus])
  
  i <- 1
  for (group_name in names(Groups)) {
    group_samples <- intersect(exp_samples, Groups[[group_name]])
    group_block <- Species[group_samples, colnames(DC_block), drop = FALSE]
    group_binary <- (group_block > 0) * 1
    group_detected <- colSums(group_binary) > 0
    
    matched_taxa <- names(group_detected[group_detected & DC_consensus])
    virus_match_matrix[matched_taxa, group_name] <- 1
  }
  
  group_counts <- colSums(virus_match_matrix[rowSums(virus_match_matrix) == 1, , drop = FALSE])
  multiple_count <- sum(rowSums(virus_match_matrix) > 1)
  
  for (group_name in names(Groups)) {
    results[group_name] <- group_counts[group_name]
  }
  
  results["Multiple"] <- multiple_count
  match_only_summary[[as.character(exp_id)]] <- results
}

match_summary_df <- do.call(cbind, match_only_summary)
match_summary_df <- as.data.frame(match_summary_df)
colnames(match_summary_df) <- paste0("Experiment_", colnames(match_summary_df))
match_summary_df$Group <- rownames(match_summary_df)

desired_order <- c("Duck_Cloacal", names(Groups), "Multiple")
match_summary_df <- match_summary_df[match_summary_df$Group %in% desired_order, ]
match_summary_df$Group <- factor(match_summary_df$Group, levels = desired_order)
long_df_match <- pivot_longer(match_summary_df, cols = starts_with("Experiment_"),
                              names_to = "Experiment", values_to = "Detections")
long_df_match$Experiment_Number <- as.numeric(gsub("Experiment_", "", long_df_match$Experiment))
long_df_match$Group_Label <- ifelse(long_df_match$Group == "Duck_Cloacal", 
                                    "Duck_Cloacal", 
                                    "Environmental")
long_df_match$Group <- factor(long_df_match$Group)
long_df_stack <- long_df_match[long_df_match$Group != "Duck_Cloacal", ]

group_palette <- c(
  "Duck_Cloacal" = "#964B00",
  "Air_Slaughter" = "#CC79A7",
  "Air_Holding" = "#6A5ACD",
  "Air_Outside" = "#228B22",
  "Cage" = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water" = "#F0E442",
  "Multiple" = "grey"
)

Duck_Cloacal_Stacked <- ggplot() +
  geom_bar(
    data = long_df_match[long_df_match$Group == "Duck_Cloacal", ],
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  geom_bar(
    data = long_df_stack,
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  facet_grid(~factor(Experiment, levels = paste0("Experiment_", sort(unique(long_df_match$Experiment_Number)))),
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y = element_text(face = "bold", size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.line.x = element_blank(),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "Detected Viruses (Duck_Cloacal + Matches)",
    y = "Number of Detections"
  )

Duck_Cloacal_Stacked

Barplot_Stacked = Chicken_Throat_Stacked / Chicken_Cloacal_Stacked / Duck_Throat_Stacked / Duck_Cloacal_Stacked
ggsave(filename="Species_Recaptured_Stacked_Barplots.pdf", plot = Barplot_Stacked, width=200, height=280, units="mm")

