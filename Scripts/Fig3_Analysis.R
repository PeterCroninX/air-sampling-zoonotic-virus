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
#### Contig Mapping Percentages - Heatmap
#================================================
Contigs <- read.xlsx("Avian_Contigs_Mapping_Percentages.xlsx")

group_palette <- c(
  "Chicken_Throat"   = "#009E73",
  "Chicken_Cloacal"  = "#800000",
  "Duck_Throat"      = "#0072B2",
  "Duck_Cloacal"     = "#964B00",
  "Air_Slaughter"    = "#CC79A7",
  "Air_Holding"      = "#6A5ACD",
  "Air_Outside"      = "#228B22",
  "Wash_Water"       = "#F0E442",
  "Drinking_Water"   = "#56B4E9",
  "Cage"             = "#D55E00"
)

Contigs$Chicken_Cloacal <- as.numeric(Contigs$Chicken_Cloacal)
Contigs_Mat <- as.matrix(Contigs[-c(1, 2, 3, 4)])      
rownames(Contigs_Mat) <- Contigs$Sample_ID             

Split <- factor(Contigs$Location_Grouping,
                levels = names(group_palette))        

col_fun1 <- colorRamp2(c(0, 100), c("white", "red"))   

row_anno <- rowAnnotation(
  Group = Split,
  col   = list(Group = group_palette),                 
  show_annotation_name = FALSE,
  border = TRUE
)

ht <- Heatmap(
  Contigs_Mat,
  name              = "% Mapping",
  col               = col_fun1,
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  row_split         = Split,           
  right_annotation  = row_anno,        
  border            = TRUE,
  row_title_rot     = 0,
  show_row_names    = FALSE,          
  column_names_rot  = 45,
  column_names_gp   = gpar(fontsize = 8, fontface = "bold")
)

draw(ht, heatmap_legend_side = "right")

pdf("Contigs_Rectangular_Heatmap_A4.pdf", width = 5, height = 11.7)
draw(ht, heatmap_legend_side = "right")
dev.off()

#================================================
#### Contig Mapping Percentages - Boxplots 
#================================================
Contigs <- read.xlsx("Avian_Contigs_Mapping_Percentages.xlsx")

group_palette <- c(
  "Chicken_Throat"   = "#009E73",
  "Chicken_Cloacal"  = "#800000",
  "Duck_Throat"      = "#0072B2",
  "Duck_Cloacal"     = "#964B00",
  "Air_Slaughter"    = "#CC79A7",
  "Air_Holding"      = "#6A5ACD",
  "Air_Outside"      = "#228B22",
  "Wash_Water"       = "#F0E442",
  "Drinking_Water"   = "#56B4E9",
  "Cage"             = "#D55E00"
)

poultry_groups <- c("Chicken_Throat", "Chicken_Cloacal",
                    "Duck_Throat",    "Duck_Cloacal")
env_groups <- setdiff(names(group_palette), poultry_groups)   

for (pg in poultry_groups) {
  df_long <- Contigs %>% 
    filter(Location_Grouping %in% env_groups) %>%              
    select(Group_Env  = Location_Grouping,                     
           PctMapped  = all_of(pg)) %>%                        
    drop_na(PctMapped) %>%                                     
    mutate(Group_Env  = factor(Group_Env, levels = env_groups))
  
  plt <- ggplot(df_long,
                aes(x = Group_Env, y = PctMapped,
                    fill = Group_Env, colour = Group_Env)) +
    geom_boxplot(outlier.shape = NA, colour = "black", width = 0.6) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    scale_fill_manual(values = group_palette) +
    scale_colour_manual(values = group_palette, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(title = pg, x = NULL,
         y = "Percentage Mapped (%)") +
    theme_minimal() +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.ticks.y     = element_line(colour = "black"),
      axis.text.y      = element_text(face = "bold", colour = "black"),
      axis.title.y     = element_text(face = "bold", colour = "black"),
      plot.title       = element_text(face = "bold", hjust = 0, colour = "black"),
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = "grey50", fill = NA, size = 1),
      legend.position  = "right",
      legend.title     = element_blank()
    )
  assign(paste0(pg, "_Pct_Boxplot"), plt, envir = .GlobalEnv)
}

Chicken <- Chicken_Throat_Pct_Boxplot / Chicken_Cloacal_Pct_Boxplot
ggsave(filename="Chicken_Percentage_Contigs_Mapped.pdf", plot = Chicken, width=200, height=280, units="mm")

Duck <- Duck_Throat_Pct_Boxplot / Duck_Cloacal_Pct_Boxplot
ggsave(filename="hDuck_Percentage_Contigs_Mapped.pdf", plot = Duck, width=200, height=280, units="mm")

Contigs$Location_Grouping <- as.factor(Contigs$Location_Grouping)
Kruskall_Chicken_Throat_Pct_Contigs <- kwAllPairsDunnTest(Contigs$Chicken_Throat, Contigs$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chicken_Throat_Pct_Contigs$p, file ="Supplementary_Files/Kruskall_Chicken_Throat_Pct_Contigs.xlsx", rowNames = TRUE)

Kruskall_Chicken_Cloacal_Pct_Contigs <- kwAllPairsDunnTest(Contigs$Chicken_Cloacal, Contigs$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chicken_Cloacal_Pct_Contigs$p, file ="Supplementary_Files/Kruskall_Chicken_Cloacal_Pct_Contigs.xlsx", rowNames = TRUE)

Kruskall_Duck_Throat_Pct_Contigs <- kwAllPairsDunnTest(Contigs$Duck_Throat, Contigs$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Duck_Throat_Pct_Contigs$p, file ="Supplementary_Files/Kruskall_Duck_Throat_Pct_Contigs.xlsx", rowNames = TRUE)

Kruskall_Duck_Cloacal_Pct_Contigs <- kwAllPairsDunnTest(Contigs$Duck_Cloacal, Contigs$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Duck_Cloacal_Pct_Contigs$p, file ="Supplementary_Files/Kruskall_Duck_Cloacal_Pct_Contigs.xlsx", rowNames = TRUE)

#================================================
#### DIAMOND - Gene - Coverage Statistics
#================================================
Species <- read.xlsx("DIAMOND/DIAMOND_PID_Coverage_Annotated.xlsx")
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Sample_Metadata <- Species[c(1)]
Sample_Metadata <- Sample_Metadata %>%
  separate(Query_Seq_ID, into = paste0("V", 1:13), sep = "\\|", remove = TRUE)
Contig_Taxonomy <- Sample_Metadata[c(11,12)]
Contig_Taxonomy$V12 <- gsub("PREDICT", "sp.", Contig_Taxonomy$V12)
Contig_Taxonomy$Contig_Species <- paste(Contig_Taxonomy$V11, Contig_Taxonomy$V12, sep = "_")
Contig_Species <- Contig_Taxonomy[c(3)]
Important_Sample_Metadata <- Sample_Metadata[c(1,2,3)]
DIAMOND <- Species[c(1,2,3,15,16,18,27,28,29)]
Final_Working_Dataset <- cbind(Important_Sample_Metadata, DIAMOND)
Final_Working_Dataset <- Final_Working_Dataset %>%
  rename_with(~ c("Sample_ID", "Group", "Experiment"), .cols = c("V1", "V2", "V3"))
#write.xlsx(Final_Working_Dataset, file="DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)

#================================================
#### DIAMOND Hits - Same Genes Identified?
#================================================
DIAMOND  <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)  
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)      
Gene_Meta <- read.xlsx("DIAMOND/Genes_Metadata_Final_Version.xlsx", rowNames = TRUE)   

Gene_Count <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>%
  distinct() %>%                   
  mutate(hit = 1) %>%
  pivot_wider(names_from  = Subject_Seq_ID,
              values_from = hit,
              values_fill = 0) %>%
  column_to_rownames("Sample_ID")

Binary_Genes <- (Gene_Count > 0) * 1  

avian_groups <- c("Chicken_Throat", "Chicken_Cloacal",
                  "Duck_Throat",    "Duck_Cloacal")

all_experiments <- sort(unique(Metadata$Experiment))
all_genes       <- colnames(Binary_Genes)

summary_matrix <- matrix(0,
                         nrow = length(all_genes),
                         ncol = length(all_experiments),
                         dimnames = list(all_genes,
                                         paste0("Experiment_", all_experiments)))

for (exp_id in all_experiments) {
  
  exp_samples_all <- rownames(Metadata[Metadata$Experiment == exp_id, ])
  exp_samples     <- intersect(exp_samples_all, rownames(Binary_Genes))
  if (length(exp_samples) == 0) next          
  
  avian_samples_all <- rownames(Metadata[Metadata$Experiment == exp_id &
                                           Metadata$Location_Grouping %in% avian_groups, ])
  avian_samples <- intersect(avian_samples_all, rownames(Binary_Genes))
  env_samples   <- setdiff(exp_samples, avian_samples)
  if (length(avian_samples) == 0 || length(env_samples) == 0) next
  
  avian_data <- Binary_Genes[avian_samples, , drop = FALSE]
  env_data   <- Binary_Genes[env_samples,   , drop = FALSE]
  
  
  avian_detected <- colSums(avian_data) > 0
  env_detected   <- colSums(env_data)   > 0
  
  col_label <- paste0("Experiment_", exp_id)
  
  summary_matrix[avian_detected &  env_detected, col_label] <- 1   
  summary_matrix[avian_detected & !env_detected, col_label] <- 2   
  summary_matrix[!avian_detected & env_detected, col_label] <- 3   
}

Label_Colors <- colorRamp2(
  c(0, 1, 2, 3),
  c("white", "#F8766D", "#F6E09B", "#00BFC4")
)

status_legend <- Legend(
  title   = "Detection",
  at      = c(1, 2, 3),
  labels  = c("Both Poultry & Env", "Poultry Only", "Environment Only"),
  legend_gp  = gpar(fill = Label_Colors(c(1, 2, 3))),
  border     = "black",
  title_gp   = gpar(fontsize = 10, fontface = "bold"),
  labels_gp  = gpar(fontsize = 9)
)

summary_matrix <- summary_matrix[, paste0("Experiment_", all_experiments)]

top_n <- 10                                   
anno_df <- Gene_Meta[rownames(summary_matrix), , drop = FALSE]
genus_freq  <- sort(table(anno_df$Genus), decreasing = TRUE)
top_genera  <- names(genus_freq)[seq_len(min(top_n, length(genus_freq)))]

anno_df$Genus_collapsed <- ifelse(anno_df$Genus %in% top_genera,
                                  anno_df$Genus, "Other")

n_cols   <- length(unique(anno_df$Genus_collapsed))   
base_pal <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(n_cols)
genus_pal <- setNames(base_pal, sort(unique(anno_df$Genus_collapsed)))

gp_full <- Gene_Meta[rownames(summary_matrix), "Gene_Product"]

rl <- rle(gp_full)                                        
starts <- cumsum(c(1, head(rl$lengths, -1)))              

gp_labels <- rep("", length(gp_full))                     
idx_keep  <- starts[ rl$lengths > 10 ]                    
gp_labels[idx_keep] <- gp_full[idx_keep]                  


row_annot <- rowAnnotation(
  Genus        = anno_df$Genus_collapsed,                 
  Gene_Product = anno_text(gp_labels, gp = gpar(fontsize = 6, just = "left")),
  col   = list(Genus = genus_pal),
  show_annotation_name = TRUE,
  annotation_width     = unit(c(4, 25), "mm")
)

heat_gene_exp <- Heatmap(
  summary_matrix,
  name              = "Detection",
  col               = Label_Colors,
  cluster_rows      = TRUE,
  cluster_columns   = FALSE,
  show_row_names    = FALSE,
  show_row_dend     = FALSE,
  show_column_names = TRUE,
  column_names_rot  = 45,
  column_names_gp   = gpar(fontface = "bold", fontsize = 8),
  right_annotation  = row_annot,
  use_raster = FALSE)

Gene_Exp_Heatmap <- draw(
  heat_gene_exp,
  padding             = unit(c(10, 20, 10, 10), "mm"),
  heatmap_legend_list = list(status_legend),
  heatmap_legend_side = "right"
)

pdf("Gene_Experiment_Summary_Heatmap_A4.pdf", width = 6, height = 12)
draw(Gene_Exp_Heatmap)
dev.off()

#================================================
#### DIAMOND Hits - Same Genes Identified - Chicken Throat
#================================================
DIAMOND   <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)
Metadata  <- read.xlsx("Metadata.xlsx",              rowNames = TRUE)
Gene_Meta <- read.xlsx("DIAMOND/Genes_Metadata_Final_Version.xlsx",         rowNames = TRUE)  

Gene_Count <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>%
  distinct() %>%
  mutate(hit = 1) %>%
  pivot_wider(names_from  = Subject_Seq_ID,
              values_from = hit,
              values_fill = 0) %>%
  column_to_rownames("Sample_ID")

Binary_Genes <- (Gene_Count > 0) * 1

Ref_Group      <- "Chicken_Throat"
poultry_groups <- c("Chicken_Throat","Chicken_Cloacal",
                    "Duck_Throat","Duck_Cloacal")

Ref_Samples <- rownames(Metadata[Metadata$Location_Grouping == Ref_Group, ])
Env_Samples <- rownames(Metadata)[!(Metadata$Location_Grouping %in% poultry_groups)]

Ref_Taxa    <- names(which(colSums(Binary_Genes[Ref_Samples, , drop = FALSE]) > 0))

Binary_Data <- Binary_Genes[c(Ref_Samples, Env_Samples), Ref_Taxa, drop = FALSE]
meta        <- Metadata[rownames(Binary_Data), ]

Labelled_Matrix <- matrix(NA,
                          nrow = nrow(Binary_Data),
                          ncol = ncol(Binary_Data),
                          dimnames = dimnames(Binary_Data))

for (exp_id in unique(meta$Experiment)) {
  exp_samples     <- rownames(meta[meta$Experiment == exp_id, ])
  Ref_Exp_Samples <- intersect(exp_samples, Ref_Samples)
  other_samples   <- setdiff(exp_samples, Ref_Exp_Samples)
  if (length(Ref_Exp_Samples) == 0 || length(other_samples) == 0) next
  
  Ref_Block     <- Binary_Data[Ref_Exp_Samples, , drop = FALSE]
  Ref_Consensus <- (colSums(Ref_Block) > 0) * 1
  
  for (sample in other_samples) {
    ref_vec   <- Ref_Consensus
    other_vec <- Binary_Data[sample, ]
    
    Labelled_Matrix[sample, ] <-
      ifelse(ref_vec == 1 & other_vec == 1, 10,
             ifelse(ref_vec == 1 & other_vec == 0, 20,
                    ifelse(ref_vec == 0 & other_vec == 1, 30, 0)))
  }
}

Labelled_Matrix <- Labelled_Matrix[meta$Location_Grouping != Ref_Group, ]
Labelled_Matrix <- Labelled_Matrix[, colSums(Labelled_Matrix != 0) > 0]
meta_filtered   <- meta[rownames(Labelled_Matrix), ]

Label_Colors <- colorRamp2(c(0, 10, 20, 30),
                           c("white", "#F8766D", "#F6E09B", "#00BFC4"))

group_palette <- c(
  "Chicken_Cloacal" = "#009E73",
  "Chicken_Throat"  = "#800000",
  "Duck_Throat"     = "#0072B2",
  "Duck_Cloacal"    = "#964B00",
  "Air_Slaughter"   = "#CC79A7",
  "Air_Holding"     = "#6A5ACD",
  "Air_Outside"     = "#228B22",
  "Wash_Water"      = "#F0E442",
  "Drinking_Water"  = "#56B4E9",
  "Cage"            = "#D55E00"
)

location_colors <- c("Phnom_Penh" = "grey", "Takeo" = "black")

all_taxa <- colnames(Labelled_Matrix)

anno_df  <- Gene_Meta[all_taxa, , drop = FALSE]

genus_freq <- sort(table(anno_df$Genus), decreasing = TRUE)
top_n      <- 10
top_genera <- names(genus_freq)[seq_len(min(top_n, length(genus_freq)))]

anno_df$Genus_collapsed <- ifelse(anno_df$Genus %in% top_genera,
                                  anno_df$Genus, "Other")

n_cols    <- length(unique(anno_df$Genus_collapsed))
genus_pal <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(n_cols),
  sort(unique(anno_df$Genus_collapsed))
)

row_annot <- rowAnnotation(
  Genus = anno_df$Genus_collapsed,
  col   = list(Genus = genus_pal),
  show_annotation_name = TRUE,
  annotation_width     = unit(4, "mm")
)

experiments  <- unique(meta_filtered$Experiment)

heatmap_list <- lapply(experiments, function(exp_id) {
  samples    <- rownames(meta_filtered[meta_filtered$Experiment == exp_id, ])
  sub_matrix <- t(Labelled_Matrix[samples, , drop = FALSE])
  
  missing <- setdiff(all_taxa, rownames(sub_matrix))
  if (length(missing))
    sub_matrix <- rbind(sub_matrix,
                        matrix(0, nrow = length(missing),
                               ncol = ncol(sub_matrix),
                               dimnames = list(missing, colnames(sub_matrix))))
  sub_matrix <- sub_matrix[all_taxa, , drop = FALSE]
  
  anno_subset <- meta_filtered[samples, ]
  
  Heatmap(
    sub_matrix,
    col               = Label_Colors,
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_row_dend     = FALSE,
    show_column_names = FALSE,
    top_annotation    = HeatmapAnnotation(
      "Sample Type" = anno_subset$Location_Grouping,
      col           = list("Sample Type" = group_palette),
      show_annotation_name = FALSE),
    bottom_annotation = HeatmapAnnotation(
      Location      = anno_subset$Location,
      col           = list(Location = location_colors),
      which         = "column",
      show_annotation_name = FALSE),
    column_title      = paste("Experiment", exp_id),
    column_title_side = "bottom",
    column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
    name              = NULL,
    show_heatmap_legend = FALSE
  )
})

combined_ht <- Reduce(`+`, c(heatmap_list, list(row_annot)))

status_legend <- Legend(
  title      = "Detection",
  at         = c(10, 20, 30),
  labels     = c("Detected", "Undetected", "ES"),
  legend_gp  = gpar(fill = Label_Colors(c(10, 20, 30))),
  border     = "black",
  title_gp   = gpar(fontsize = 10, fontface = "bold"),
  labels_gp  = gpar(fontsize = 9)
)

final_plot <- draw(
  combined_ht,
  padding             = unit(c(10, 10, 10, 10), "mm"),
  heatmap_legend_list = list(status_legend),
  heatmap_legend_side = "right"
)

pdf("Genes_Chicken_Throat_Heatmap_A4.pdf", width = 8, height = 11.5)
draw(final_plot)
dev.off()

#================================================
#### DIAMOND Hits - Same Genes Identified - Chicken Cloacal
#================================================
DIAMOND   <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)
Metadata  <- read.xlsx("Metadata.xlsx",              rowNames = TRUE)
Gene_Meta <- read.xlsx("DIAMOND/Genes_Metadata_Final_Version.xlsx",         rowNames = TRUE)  

Gene_Count <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>%
  distinct() %>%
  mutate(hit = 1) %>%
  pivot_wider(names_from  = Subject_Seq_ID,
              values_from = hit,
              values_fill = 0) %>%
  column_to_rownames("Sample_ID")

Binary_Genes <- (Gene_Count > 0) * 1

Ref_Group      <- "Chicken_Cloacal"
poultry_groups <- c("Chicken_Throat","Chicken_Cloacal",
                    "Duck_Throat","Duck_Cloacal")

Ref_Samples <- rownames(Metadata[Metadata$Location_Grouping == Ref_Group, ])
Env_Samples <- rownames(Metadata)[!(Metadata$Location_Grouping %in% poultry_groups)]

Ref_Taxa    <- names(which(colSums(Binary_Genes[Ref_Samples, , drop = FALSE]) > 0))

Binary_Data <- Binary_Genes[c(Ref_Samples, Env_Samples), Ref_Taxa, drop = FALSE]
meta        <- Metadata[rownames(Binary_Data), ]

Labelled_Matrix <- matrix(NA,
                          nrow = nrow(Binary_Data),
                          ncol = ncol(Binary_Data),
                          dimnames = dimnames(Binary_Data))

for (exp_id in unique(meta$Experiment)) {
  exp_samples     <- rownames(meta[meta$Experiment == exp_id, ])
  Ref_Exp_Samples <- intersect(exp_samples, Ref_Samples)
  other_samples   <- setdiff(exp_samples, Ref_Exp_Samples)
  if (length(Ref_Exp_Samples) == 0 || length(other_samples) == 0) next
  
  Ref_Block     <- Binary_Data[Ref_Exp_Samples, , drop = FALSE]
  Ref_Consensus <- (colSums(Ref_Block) > 0) * 1
  
  for (sample in other_samples) {
    ref_vec   <- Ref_Consensus
    other_vec <- Binary_Data[sample, ]
    
    Labelled_Matrix[sample, ] <-
      ifelse(ref_vec == 1 & other_vec == 1, 10,
             ifelse(ref_vec == 1 & other_vec == 0, 20,
                    ifelse(ref_vec == 0 & other_vec == 1, 30, 0)))
  }
}

Labelled_Matrix <- Labelled_Matrix[meta$Location_Grouping != Ref_Group, ]
Labelled_Matrix <- Labelled_Matrix[, colSums(Labelled_Matrix != 0) > 0]
meta_filtered   <- meta[rownames(Labelled_Matrix), ]

Label_Colors <- colorRamp2(c(0, 10, 20, 30),
                           c("white", "#F8766D", "#F6E09B", "#00BFC4"))

group_palette <- c(
  "Chicken_Cloacal" = "#009E73",
  "Chicken_Throat"  = "#800000",
  "Duck_Throat"     = "#0072B2",
  "Duck_Cloacal"    = "#964B00",
  "Air_Slaughter"   = "#CC79A7",
  "Air_Holding"     = "#6A5ACD",
  "Air_Outside"     = "#228B22",
  "Wash_Water"      = "#F0E442",
  "Drinking_Water"  = "#56B4E9",
  "Cage"            = "#D55E00"
)

location_colors <- c("Phnom_Penh" = "grey", "Takeo" = "black")

all_taxa <- colnames(Labelled_Matrix)

anno_df  <- Gene_Meta[all_taxa, , drop = FALSE]

genus_freq <- sort(table(anno_df$Genus), decreasing = TRUE)
top_n      <- 10
top_genera <- names(genus_freq)[seq_len(min(top_n, length(genus_freq)))]

anno_df$Genus_collapsed <- ifelse(anno_df$Genus %in% top_genera,
                                  anno_df$Genus, "Other")

n_cols    <- length(unique(anno_df$Genus_collapsed))
genus_pal <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(n_cols),
  sort(unique(anno_df$Genus_collapsed))
)

row_annot <- rowAnnotation(
  Genus = anno_df$Genus_collapsed,
  col   = list(Genus = genus_pal),
  show_annotation_name = TRUE,
  annotation_width     = unit(4, "mm")
)

experiments  <- unique(meta_filtered$Experiment)

heatmap_list <- lapply(experiments, function(exp_id) {
  samples    <- rownames(meta_filtered[meta_filtered$Experiment == exp_id, ])
  sub_matrix <- t(Labelled_Matrix[samples, , drop = FALSE])
  
  missing <- setdiff(all_taxa, rownames(sub_matrix))
  if (length(missing))
    sub_matrix <- rbind(sub_matrix,
                        matrix(0, nrow = length(missing),
                               ncol = ncol(sub_matrix),
                               dimnames = list(missing, colnames(sub_matrix))))
  sub_matrix <- sub_matrix[all_taxa, , drop = FALSE]
  
  anno_subset <- meta_filtered[samples, ]
  
  Heatmap(
    sub_matrix,
    col               = Label_Colors,
    cluster_rows      = TRUE,
    show_row_dend     = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = FALSE,
    top_annotation    = HeatmapAnnotation(
      "Sample Type" = anno_subset$Location_Grouping,
      col           = list("Sample Type" = group_palette),
      show_annotation_name = FALSE),
    bottom_annotation = HeatmapAnnotation(
      Location      = anno_subset$Location,
      col           = list(Location = location_colors),
      which         = "column",
      show_annotation_name = FALSE),
    column_title      = paste("Experiment", exp_id),
    column_title_side = "bottom",
    column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
    name              = NULL,
    show_heatmap_legend = FALSE
  )
})

combined_ht <- Reduce(`+`, c(heatmap_list, list(row_annot)))

status_legend <- Legend(
  title      = "Detection",
  at         = c(10, 20, 30),
  labels     = c("Detected", "Undetected", "ES"),
  legend_gp  = gpar(fill = Label_Colors(c(10, 20, 30))),
  border     = "black",
  title_gp   = gpar(fontsize = 10, fontface = "bold"),
  labels_gp  = gpar(fontsize = 9)
)

final_plot <- draw(
  combined_ht,
  padding             = unit(c(10, 10, 10, 10), "mm"),
  heatmap_legend_list = list(status_legend),
  heatmap_legend_side = "right"
)

pdf("Genes_Chicken_Cloacal_Heatmap_A4.pdf", width = 8, height = 11.5)
draw(final_plot)
dev.off()

#================================================
#### DIAMOND Hits - Same Genes Identified - Duck Throat
#================================================
DIAMOND   <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)
Metadata  <- read.xlsx("Metadata.xlsx",              rowNames = TRUE)
Gene_Meta <- read.xlsx("DIAMOND/Genes_Metadata_Final_Version.xlsx",         rowNames = TRUE)  

Gene_Count <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>%
  distinct() %>%
  mutate(hit = 1) %>%
  pivot_wider(names_from  = Subject_Seq_ID,
              values_from = hit,
              values_fill = 0) %>%
  column_to_rownames("Sample_ID")

Binary_Genes <- (Gene_Count > 0) * 1

Ref_Group      <- "Duck_Throat"
poultry_groups <- c("Chicken_Throat","Chicken_Cloacal",
                    "Duck_Throat","Duck_Cloacal")

Ref_Samples <- rownames(Metadata[Metadata$Location_Grouping == Ref_Group, ])
Env_Samples <- rownames(Metadata)[!(Metadata$Location_Grouping %in% poultry_groups)]

Ref_Taxa    <- names(which(colSums(Binary_Genes[Ref_Samples, , drop = FALSE]) > 0))

Binary_Data <- Binary_Genes[c(Ref_Samples, Env_Samples), Ref_Taxa, drop = FALSE]
meta        <- Metadata[rownames(Binary_Data), ]

Labelled_Matrix <- matrix(NA,
                          nrow = nrow(Binary_Data),
                          ncol = ncol(Binary_Data),
                          dimnames = dimnames(Binary_Data))

for (exp_id in unique(meta$Experiment)) {
  exp_samples     <- rownames(meta[meta$Experiment == exp_id, ])
  Ref_Exp_Samples <- intersect(exp_samples, Ref_Samples)
  other_samples   <- setdiff(exp_samples, Ref_Exp_Samples)
  if (length(Ref_Exp_Samples) == 0 || length(other_samples) == 0) next
  
  Ref_Block     <- Binary_Data[Ref_Exp_Samples, , drop = FALSE]
  Ref_Consensus <- (colSums(Ref_Block) > 0) * 1
  
  for (sample in other_samples) {
    ref_vec   <- Ref_Consensus
    other_vec <- Binary_Data[sample, ]
    
    Labelled_Matrix[sample, ] <-
      ifelse(ref_vec == 1 & other_vec == 1, 10,
             ifelse(ref_vec == 1 & other_vec == 0, 20,
                    ifelse(ref_vec == 0 & other_vec == 1, 30, 0)))
  }
}

Labelled_Matrix <- Labelled_Matrix[meta$Location_Grouping != Ref_Group, ]
Labelled_Matrix <- Labelled_Matrix[, colSums(Labelled_Matrix != 0) > 0]
meta_filtered   <- meta[rownames(Labelled_Matrix), ]

Label_Colors <- colorRamp2(c(0, 10, 20, 30),
                           c("white", "#F8766D", "#F6E09B", "#00BFC4"))

group_palette <- c(
  "Chicken_Cloacal" = "#009E73",
  "Chicken_Throat"  = "#800000",
  "Duck_Throat"     = "#0072B2",
  "Duck_Cloacal"    = "#964B00",
  "Air_Slaughter"   = "#CC79A7",
  "Air_Holding"     = "#6A5ACD",
  "Air_Outside"     = "#228B22",
  "Wash_Water"      = "#F0E442",
  "Drinking_Water"  = "#56B4E9",
  "Cage"            = "#D55E00"
)

location_colors <- c("Phnom_Penh" = "grey", "Takeo" = "black")

all_taxa <- colnames(Labelled_Matrix)

anno_df  <- Gene_Meta[all_taxa, , drop = FALSE]

genus_freq <- sort(table(anno_df$Genus), decreasing = TRUE)
top_n      <- 10
top_genera <- names(genus_freq)[seq_len(min(top_n, length(genus_freq)))]

anno_df$Genus_collapsed <- ifelse(anno_df$Genus %in% top_genera,
                                  anno_df$Genus, "Other")

n_cols    <- length(unique(anno_df$Genus_collapsed))
genus_pal <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(n_cols),
  sort(unique(anno_df$Genus_collapsed))
)

row_annot <- rowAnnotation(
  Genus = anno_df$Genus_collapsed,
  col   = list(Genus = genus_pal),
  show_annotation_name = TRUE,
  annotation_width     = unit(4, "mm")
)

experiments  <- unique(meta_filtered$Experiment)

heatmap_list <- lapply(experiments, function(exp_id) {
  samples    <- rownames(meta_filtered[meta_filtered$Experiment == exp_id, ])
  sub_matrix <- t(Labelled_Matrix[samples, , drop = FALSE])
  
  missing <- setdiff(all_taxa, rownames(sub_matrix))
  if (length(missing))
    sub_matrix <- rbind(sub_matrix,
                        matrix(0, nrow = length(missing),
                               ncol = ncol(sub_matrix),
                               dimnames = list(missing, colnames(sub_matrix))))
  sub_matrix <- sub_matrix[all_taxa, , drop = FALSE]
  
  anno_subset <- meta_filtered[samples, ]
  
  Heatmap(
    sub_matrix,
    col               = Label_Colors,
    cluster_rows      = TRUE,
    show_row_dend     = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = FALSE,
    top_annotation    = HeatmapAnnotation(
      "Sample Type" = anno_subset$Location_Grouping,
      col           = list("Sample Type" = group_palette),
      show_annotation_name = FALSE),
    bottom_annotation = HeatmapAnnotation(
      Location      = anno_subset$Location,
      col           = list(Location = location_colors),
      which         = "column",
      show_annotation_name = FALSE),
    column_title      = paste("Experiment", exp_id),
    column_title_side = "bottom",
    column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
    name              = NULL,
    show_heatmap_legend = FALSE
  )
})

combined_ht <- Reduce(`+`, c(heatmap_list, list(row_annot)))

status_legend <- Legend(
  title      = "Detection",
  at         = c(10, 20, 30),
  labels     = c("Detected", "Undetected", "ES"),
  legend_gp  = gpar(fill = Label_Colors(c(10, 20, 30))),
  border     = "black",
  title_gp   = gpar(fontsize = 10, fontface = "bold"),
  labels_gp  = gpar(fontsize = 9)
)

final_plot <- draw(
  combined_ht,
  padding             = unit(c(10, 10, 10, 10), "mm"),
  heatmap_legend_list = list(status_legend),
  heatmap_legend_side = "right"
)

pdf("Genes_Duck_Throat_Heatmap_A4.pdf", width = 8, height = 11.5)
draw(final_plot)
dev.off()

#================================================
#### DIAMOND Hits - Same Genes Identified - Duck Cloacal
#================================================
DIAMOND   <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)
Metadata  <- read.xlsx("Metadata.xlsx",              rowNames = TRUE)
Gene_Meta <- read.xlsx("DIAMOND/Genes_Metadata_Final_Version.xlsx",         rowNames = TRUE)  

Gene_Count <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>%
  distinct() %>%
  mutate(hit = 1) %>%
  pivot_wider(names_from  = Subject_Seq_ID,
              values_from = hit,
              values_fill = 0) %>%
  column_to_rownames("Sample_ID")

Binary_Genes <- (Gene_Count > 0) * 1

Ref_Group      <- "Duck_Cloacal"
poultry_groups <- c("Chicken_Throat","Chicken_Cloacal",
                    "Duck_Throat","Duck_Cloacal")

Ref_Samples <- rownames(Metadata[Metadata$Location_Grouping == Ref_Group, ])
Env_Samples <- rownames(Metadata)[!(Metadata$Location_Grouping %in% poultry_groups)]

Ref_Taxa    <- names(which(colSums(Binary_Genes[Ref_Samples, , drop = FALSE]) > 0))

Binary_Data <- Binary_Genes[c(Ref_Samples, Env_Samples), Ref_Taxa, drop = FALSE]
meta        <- Metadata[rownames(Binary_Data), ]

Labelled_Matrix <- matrix(NA,
                          nrow = nrow(Binary_Data),
                          ncol = ncol(Binary_Data),
                          dimnames = dimnames(Binary_Data))

for (exp_id in unique(meta$Experiment)) {
  exp_samples     <- rownames(meta[meta$Experiment == exp_id, ])
  Ref_Exp_Samples <- intersect(exp_samples, Ref_Samples)
  other_samples   <- setdiff(exp_samples, Ref_Exp_Samples)
  if (length(Ref_Exp_Samples) == 0 || length(other_samples) == 0) next
  
  Ref_Block     <- Binary_Data[Ref_Exp_Samples, , drop = FALSE]
  Ref_Consensus <- (colSums(Ref_Block) > 0) * 1
  
  for (sample in other_samples) {
    ref_vec   <- Ref_Consensus
    other_vec <- Binary_Data[sample, ]
    
    Labelled_Matrix[sample, ] <-
      ifelse(ref_vec == 1 & other_vec == 1, 10,
             ifelse(ref_vec == 1 & other_vec == 0, 20,
                    ifelse(ref_vec == 0 & other_vec == 1, 30, 0)))
  }
}

Labelled_Matrix <- Labelled_Matrix[meta$Location_Grouping != Ref_Group, ]
Labelled_Matrix <- Labelled_Matrix[, colSums(Labelled_Matrix != 0) > 0]
meta_filtered   <- meta[rownames(Labelled_Matrix), ]

Label_Colors <- colorRamp2(c(0, 10, 20, 30),
                           c("white", "#F8766D", "#F6E09B", "#00BFC4"))

group_palette <- c(
  "Chicken_Cloacal" = "#009E73",
  "Chicken_Throat"  = "#800000",
  "Duck_Throat"     = "#0072B2",
  "Duck_Cloacal"    = "#964B00",
  "Air_Slaughter"   = "#CC79A7",
  "Air_Holding"     = "#6A5ACD",
  "Air_Outside"     = "#228B22",
  "Wash_Water"      = "#F0E442",
  "Drinking_Water"  = "#56B4E9",
  "Cage"            = "#D55E00"
)

location_colors <- c("Phnom_Penh" = "grey", "Takeo" = "black")

all_taxa <- colnames(Labelled_Matrix)

anno_df  <- Gene_Meta[all_taxa, , drop = FALSE]

genus_freq <- sort(table(anno_df$Genus), decreasing = TRUE)
top_n      <- 10
top_genera <- names(genus_freq)[seq_len(min(top_n, length(genus_freq)))]

anno_df$Genus_collapsed <- ifelse(anno_df$Genus %in% top_genera,
                                  anno_df$Genus, "Other")

n_cols    <- length(unique(anno_df$Genus_collapsed))
genus_pal <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(n_cols),
  sort(unique(anno_df$Genus_collapsed))
)

row_annot <- rowAnnotation(
  Genus = anno_df$Genus_collapsed,
  col   = list(Genus = genus_pal),
  show_annotation_name = TRUE,
  annotation_width     = unit(4, "mm")
)

experiments  <- unique(meta_filtered$Experiment)

heatmap_list <- lapply(experiments, function(exp_id) {
  samples    <- rownames(meta_filtered[meta_filtered$Experiment == exp_id, ])
  sub_matrix <- t(Labelled_Matrix[samples, , drop = FALSE])
  
  missing <- setdiff(all_taxa, rownames(sub_matrix))
  if (length(missing))
    sub_matrix <- rbind(sub_matrix,
                        matrix(0, nrow = length(missing),
                               ncol = ncol(sub_matrix),
                               dimnames = list(missing, colnames(sub_matrix))))
  sub_matrix <- sub_matrix[all_taxa, , drop = FALSE]
  
  anno_subset <- meta_filtered[samples, ]
  
  Heatmap(
    sub_matrix,
    col               = Label_Colors,
    cluster_rows      = TRUE,
    show_row_dend     = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = TRUE,
    show_column_names = FALSE,
    top_annotation    = HeatmapAnnotation(
      "Sample Type" = anno_subset$Location_Grouping,
      col           = list("Sample Type" = group_palette),
      show_annotation_name = FALSE),
    bottom_annotation = HeatmapAnnotation(
      Location      = anno_subset$Location,
      col           = list(Location = location_colors),
      which         = "column",
      show_annotation_name = FALSE),
    column_title      = paste("Experiment", exp_id),
    column_title_side = "bottom",
    column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
    name              = NULL,
    show_heatmap_legend = FALSE
  )
})

combined_ht <- Reduce(`+`, c(heatmap_list, list(row_annot)))

status_legend <- Legend(
  title      = "Detection",
  at         = c(10, 20, 30),
  labels     = c("Detected", "Undetected", "ES"),
  legend_gp  = gpar(fill = Label_Colors(c(10, 20, 30))),
  border     = "black",
  title_gp   = gpar(fontsize = 10, fontface = "bold"),
  labels_gp  = gpar(fontsize = 9)
)

final_plot <- draw(
  combined_ht,
  padding             = unit(c(10, 10, 10, 10), "mm"),
  heatmap_legend_list = list(status_legend),
  heatmap_legend_side = "right"
)

pdf("Genes_Duck_Cloacal_Heatmap_A4.pdf", width = 8, height = 11.5)
draw(final_plot)
dev.off()

#================================================
#### DIAMOND Gene Hits - Barplot/Boxplot Summary - Chicken Throat
#================================================
DIAMOND  <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames=TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames=TRUE)

Gene_Count <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>% distinct() %>%
  mutate(hit=1) %>%
  pivot_wider(names_from=Subject_Seq_ID, values_from=hit, values_fill=0) %>%
  column_to_rownames("Sample_ID")
Binary_Genes <- (Gene_Count>0)*1

CT_samples <- Metadata %>% 
  filter(Location_Grouping=="Chicken_Throat") %>% rownames()

Groups <- list(
  Air_Slaughter   = Metadata %>% filter(Location_Grouping=="Air_Slaughter")   %>% rownames(),
  Air_Holding     = Metadata %>% filter(Location_Grouping=="Air_Holding")     %>% rownames(),
  Air_Outside     = Metadata %>% filter(Location_Grouping=="Air_Outside")     %>% rownames(),
  Cage            = Metadata %>% filter(Location_Grouping=="Cage")            %>% rownames(),
  Drinking_Water  = Metadata %>% filter(Location_Grouping=="Drinking_Water")  %>% rownames(),
  Wash_Water      = Metadata %>% filter(Location_Grouping=="Wash_Water")      %>% rownames()
)

CT_mat  <- Binary_Genes[intersect(CT_samples, rownames(Binary_Genes)), , drop=FALSE]
CT_mat  <- CT_mat[, colSums(CT_mat)>0, drop=FALSE]
CT_Taxa <- colnames(CT_mat)

summary_list <- list()
for(exp_id in sort(unique(Metadata$Experiment))){
  exp_meta  <- Metadata %>% filter(Experiment==exp_id)
  CT_exp    <- intersect(rownames(exp_meta), CT_samples)
  CT_block  <- Binary_Genes[CT_exp, CT_Taxa, drop=FALSE]
  if(nrow(CT_block)==0) next
  CT_cons   <- colSums(CT_block)>0
  
  res <- tibble(
    Experiment = exp_id,
    Group      = "Chicken_Throat",
    Type       = "Chicken_Throat",
    Detections = sum(CT_cons)
  )
  
  for(g in names(Groups)){
    g_exp   <- intersect(rownames(exp_meta), Groups[[g]])
    g_block <- Binary_Genes[g_exp, CT_Taxa, drop=FALSE]
    if(nrow(g_block)==0) next
    g_cons  <- colSums(g_block)>0
    
    res <- bind_rows(res,
                     tibble(Experiment=exp_id, Group=g,    Type="Match",  Detections=sum(CT_cons & g_cons)),
                     tibble(Experiment=exp_id, Group=g,    Type="Unique", Detections=sum(!CT_cons & g_cons))
    )
  }
  
  summary_list[[as.character(exp_id)]] <- res
}
summary_df <- bind_rows(summary_list)

collected_groups <- Metadata %>%
  group_by(Experiment, Location_Grouping) %>%
  summarise(n=n(), .groups="drop")

long_match <- summary_df %>%
  filter(Type=="Match" | Type=="Chicken_Throat") %>%
  inner_join(collected_groups, by=c("Experiment","Group"="Location_Grouping")) %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=c("Chicken_Throat", names(Groups)))
  )

long_unique <- summary_df %>%
  filter(Type=="Unique") %>%
  inner_join(collected_groups, by=c("Experiment","Group"="Location_Grouping")) %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=names(Groups))
  )

group_palette <- c(
  "Chicken_Throat"  = "#009E73",
  "Air_Slaughter"   = "#CC79A7",
  "Air_Holding"     = "#6A5ACD",
  "Air_Outside"     = "#228B22",
  "Wash_Water"      = "#F0E442",
  "Drinking_Water"  = "#56B4E9",
  "Cage"            = "#D55E00"
)

Chicken_Throat_Gene_p_match <- ggplot(long_match, aes(Group, Detections, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Chicken Throat Matches", y="Detections", x=NULL) +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(face="bold", size=12, color="black"),
    axis.title.y = element_text(face="bold", size=13, color="black"),
    axis.line.y = element_line(color="black", linewidth=0.8),
    strip.text = element_text(size=13, face="bold", color="black"),
    plot.title = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  )

Chicken_Throat_Gene_p_unique <- ggplot(long_unique, aes(Group, Detections, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Chicken Throat Unique Detections", y="Detections", x=NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face="bold", size=12, color="black"),
    axis.title.y = element_text(face="bold", size=13, color="black"),
    axis.line.y = element_line(color="black", linewidth=0.8),
    strip.text = element_text(size=13, face="bold", color="black"),
    plot.title = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  )

box_df <- summary_df %>%
  mutate(
    Combined = paste(Group, Type, sep = "_"),
    Combined = factor(Combined, levels = c(
      "Chicken_Throat_Chicken_Throat",
      "Air_Slaughter_Match",   "Air_Slaughter_Unique",
      "Air_Holding_Match",     "Air_Holding_Unique",
      "Air_Outside_Match",     "Air_Outside_Unique",
      "Wash_Water_Match",      "Wash_Water_Unique",
      "Drinking_Water_Match",  "Drinking_Water_Unique",
      "Cage_Match",            "Cage_Unique"
    ))
  )

box_df$Type[ box_df$Type == "Chicken_Throat" ] <- "Match"
box_plot_df <- box_df %>% 
  filter(Type %in% c("Match", "Unique"))

box_plot_df$Group <- factor(
  box_plot_df$Group,
  levels = c(
    "Chicken_Throat",
    "Air_Slaughter",
    "Air_Holding",
    "Air_Outside",
    "Wash_Water",
    "Drinking_Water",
    "Cage"
  ))

Kruskall_Genes_Chicken_Throat_Recapture_Unique <- kwAllPairsDunnTest(box_plot_df$Detections, box_plot_df$Combined, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Genes_Chicken_Throat_Recapture_Unique$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Genes_Chicken_Throat_Recapture_Unique.xlsx", rowNames = TRUE)

Chicken_Throat_Gene_p_box <- ggplot(box_plot_df, aes(x = Combined, y = Detections, fill = Group)) +
  geom_boxplot_pattern(
    aes(pattern = Type),
    outlier.shape   = NA,
    pattern_fill    = "black",
    pattern_angle   = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.02,
    color           = "black",
    width           = 0.6,
    show.legend     = TRUE
  ) +
  geom_jitter(
    aes(color = Group),
    width = 0.2,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_palette) +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_pattern_manual(
    name   = "Detection Type",
    values = c(
      Chicken_Throat = "none",  
      Match          = "none",
      Unique         = "stripe"
    ),
    breaks = c("Match", "Unique"),  
    guide  = guide_legend(
      override.aes = list(
        fill            = c("grey80", "grey80"),
        pattern         = c("none",    "stripe"),
        pattern_fill    = c(NA,        "black"),
        pattern_angle   = c(0,         45),
        pattern_density = c(0,         0.1),
        pattern_spacing = c(0,         0.02),
        color           = "black"
      )
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    title = "Match vs Unique Boxplots",
    x     = NULL,
    y     = "Detections"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color = "black"),
    axis.text.y      = element_text(face = "bold", color = "black", size = 12),
    axis.title.y     = element_text(face = "bold", color = "black", size = 13),
    plot.title       = element_text(face = "bold", hjust = 0, color = "black", size = 13),
    axis.line.x      = element_blank(),
    axis.line.y      = element_blank(),
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(color = "grey50", fill = NA, size = 1),
    legend.title     = element_text(face = "bold", color = "black"),
    legend.text      = element_text(color = "black"),
    legend.position  = "right"
  )

total_df <- summary_df %>%
  group_by(Experiment, Group) %>%
  summarise(
    Total = sum(Detections),
    .groups = "drop"
  ) %>%
  mutate(
    Experiment = factor(Experiment, levels = sort(unique(Experiment))),
    Group      = factor(Group, levels = c("Chicken_Throat", names(Groups)))
  )

Kruskall_Genes_Chicken_Throat_Total <- kwAllPairsDunnTest(total_df$Total, total_df$Group, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Genes_Chicken_Throat_Total$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Genes_Chicken_Throat_Total.xlsx", rowNames = TRUE)

Chicken_Throat_Gene_p_total_box <- ggplot(total_df, aes(x=Group, y=Total, fill=Group)) +
  geom_boxplot(outlier.shape=NA, color="black", width=0.6) +
  geom_jitter(aes(color=Group), width=0.2, alpha=0.6) +
  scale_fill_manual(values=group_palette) +
  scale_color_manual(values=group_palette, guide="none") +
  labs(title="Chicken Throat Total Detections", y="Detections", x=NULL) +
  theme_minimal() +
  theme(
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.ticks.y    = element_line(color = "black"), 
    axis.text.y     = element_text(face = "bold", color = "black"),
    axis.title.y    = element_text(face = "bold", color = "black"),
    plot.title      = element_text(face = "bold", hjust = 0, color = "black"),
    axis.line.x     = element_blank(),
    axis.line.y     = element_blank(),
    panel.grid      = element_blank(),
    panel.background= element_blank(),
    panel.border    = element_rect(color = "grey50", fill = NA, size = 1),
    legend.position = "right"
  )

Chicken_Throat_Gene_p_total_bar <- ggplot(total_df, aes(x=Group, y=Total, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Detected Viruses (Collapsed Groups)", y="Detections", x=NULL) +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y    = element_line(color = "black"),
    axis.text.y = element_text(face="bold", size=12, color="black"),
    axis.title.y = element_text(face="bold", size=13, color="black"),
    axis.line.y = element_line(color="black", linewidth=0.8),
    strip.text = element_text(size=13, face="bold", color="black"),
    plot.title = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  )

print(Chicken_Throat_Gene_p_match)
print(Chicken_Throat_Gene_p_unique)
print(Chicken_Throat_Gene_p_box)
print(Chicken_Throat_Gene_p_total_box)
print(Chicken_Throat_Gene_p_total_bar)

#================================================
#### DIAMOND Gene Hits - Barplot/Boxplot Summary - Chicken Cloacal
#================================================
DIAMOND  <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames=TRUE)
Metadata <- read.xlsx("Metadata.xlsx",                rowNames=TRUE)

Gene_Count <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>% distinct() %>%
  mutate(hit=1) %>%
  pivot_wider(names_from=Subject_Seq_ID, values_from=hit, values_fill=0) %>%
  column_to_rownames("Sample_ID")
Binary_Genes <- (Gene_Count>0)*1

CC_samples <- Metadata %>% 
  filter(Location_Grouping=="Chicken_Cloacal") %>% rownames()

Groups <- list(
  Air_Slaughter   = Metadata %>% filter(Location_Grouping=="Air_Slaughter")   %>% rownames(),
  Air_Holding     = Metadata %>% filter(Location_Grouping=="Air_Holding")     %>% rownames(),
  Air_Outside     = Metadata %>% filter(Location_Grouping=="Air_Outside")     %>% rownames(),
  Cage            = Metadata %>% filter(Location_Grouping=="Cage")            %>% rownames(),
  Drinking_Water  = Metadata %>% filter(Location_Grouping=="Drinking_Water")  %>% rownames(),
  Wash_Water      = Metadata %>% filter(Location_Grouping=="Wash_Water")      %>% rownames()
)

CC_mat  <- Binary_Genes[intersect(CC_samples, rownames(Binary_Genes)), , drop=FALSE]
CC_mat  <- CC_mat[, colSums(CC_mat)>0, drop=FALSE]
CC_Taxa <- colnames(CC_mat)

summary_list <- list()
for(exp_id in sort(unique(Metadata$Experiment))){
  exp_meta  <- Metadata %>% filter(Experiment==exp_id)
  CC_exp    <- intersect(rownames(exp_meta), CC_samples)
  CC_block  <- Binary_Genes[CC_exp, CC_Taxa, drop=FALSE]
  if(nrow(CC_block)==0) next
  CC_cons   <- colSums(CC_block)>0
  
  res <- tibble(
    Experiment = exp_id,
    Group      = "Chicken_Cloacal",
    Type       = "Chicken_Cloacal",
    Detections = sum(CC_cons)
  )
  
  for(g in names(Groups)){
    g_exp   <- intersect(rownames(exp_meta), Groups[[g]])
    g_block <- Binary_Genes[g_exp, CC_Taxa, drop=FALSE]
    if(nrow(g_block)==0) next
    g_cons  <- colSums(g_block)>0
    
    res <- bind_rows(res,
                     tibble(Experiment=exp_id, Group=g,    Type="Match",  Detections=sum(CC_cons & g_cons)),
                     tibble(Experiment=exp_id, Group=g,    Type="Unique", Detections=sum(!CC_cons & g_cons))
    )
  }
  
  summary_list[[as.character(exp_id)]] <- res
}
summary_df <- bind_rows(summary_list)

collected_groups <- Metadata %>%
  group_by(Experiment, Location_Grouping) %>%
  summarise(n=n(), .groups="drop")

long_match <- summary_df %>%
  filter(Type=="Match" | Type=="Chicken_Cloacal") %>%
  inner_join(collected_groups, by=c("Experiment","Group"="Location_Grouping")) %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=c("Chicken_Cloacal", names(Groups)))
  )

long_unique <- summary_df %>%
  filter(Type=="Unique") %>%
  inner_join(collected_groups, by=c("Experiment","Group"="Location_Grouping")) %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=names(Groups))
  )

group_palette <- c(
  "Chicken_Cloacal" = "#800000",
  "Air_Slaughter"   = "#CC79A7",
  "Air_Holding"     = "#6A5ACD",
  "Air_Outside"     = "#228B22",
  "Wash_Water"      = "#F0E442",
  "Drinking_Water"  = "#56B4E9",
  "Cage"            = "#D55E00"
)

Chicken_Cloacal_Gene_p_match <- ggplot(long_match, aes(Group, Detections, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Chicken Cloacal Matches", y="Detections", x=NULL) +
  theme_minimal() +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color="black"),
    axis.text.y      = element_text(face="bold", size=12, color="black"),
    axis.title.y     = element_text(face="bold", size=13, color="black"),
    axis.line.y      = element_line(color="black", linewidth=0.8),
    strip.text       = element_text(size=13, face="bold", color="black"),
    plot.title       = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )

Chicken_Cloacal_Gene_p_unique <- ggplot(long_unique, aes(Group, Detections, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Chicken Cloacal Unique Detections", y="Detections", x=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face="bold", size=12, color="black"),
    axis.title.y     = element_text(face="bold", size=13, color="black"),
    axis.line.y      = element_line(color="black", linewidth=0.8),
    strip.text       = element_text(size=13, face="bold", color="black"),
    plot.title       = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )

box_df <- summary_df %>%
  mutate(
    Combined = paste(Group, Type, sep="_"),
    Combined = factor(Combined, levels=c(
      "Chicken_Cloacal_Chicken_Cloacal",
      "Air_Slaughter_Match",   "Air_Slaughter_Unique",
      "Air_Holding_Match",     "Air_Holding_Unique",
      "Air_Outside_Match",     "Air_Outside_Unique",
      "Wash_Water_Match",      "Wash_Water_Unique",
      "Drinking_Water_Match",  "Drinking_Water_Unique",
      "Cage_Match",            "Cage_Unique"
    ))
  )

box_df$Type[box_df$Type=="Chicken_Cloacal"] <- "Match"
box_plot_df <- box_df %>% filter(Type %in% c("Match","Unique"))

box_plot_df$Group <- factor(
  box_plot_df$Group,
  levels = c("Chicken_Cloacal", names(Groups))
)

Kruskall_Genes_Chicken_Cloacal_Recapture_Unique <- kwAllPairsDunnTest(box_plot_df$Detections, box_plot_df$Combined, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Genes_Chicken_Cloacal_Recapture_Unique$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Genes_Chicken_Cloacal_Recapture_Unique.xlsx", rowNames = TRUE)

Chicken_Cloacal_Gene_p_box <- ggplot(box_plot_df, aes(x=Combined, y=Detections, fill=Group)) +
  geom_boxplot_pattern(
    aes(pattern=Type),
    outlier.shape   = NA,
    pattern_fill    = "black",
    pattern_angle   = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.02,
    color           = "black",
    width           = 0.6,
    show.legend     = TRUE
  ) +
  geom_jitter(aes(color=Group), width=0.2, alpha=0.6, show.legend=FALSE) +
  scale_fill_manual(values=group_palette) +
  scale_color_manual(values=group_palette, guide="none") +
  scale_pattern_manual(
    name   = "Detection Type",
    values = c(Match="none", Unique="stripe"),
    breaks = c("Match","Unique"),
    guide  = guide_legend(
      override.aes = list(
        fill            = c("grey80","grey80"),
        pattern         = c("none","stripe"),
        pattern_fill    = c(NA,"black"),
        pattern_angle   = c(0,45),
        pattern_density = c(0,0.1),
        pattern_spacing = c(0,0.02),
        color           = "black"
      )
    )
  ) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.05))) +
  labs(title="Match vs Unique Boxplots", x=NULL, y="Detections") +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(face="bold", color="black", size=12),
    axis.ticks.x     = element_line(color="black"),
    axis.ticks.y     = element_line(color="black"),
    axis.text.y      = element_text(face="bold", color="black", size=12),
    axis.title.y     = element_text(face="bold", color="black", size=13),
    plot.title       = element_text(face="bold", hjust=0, color="black", size=13),
    axis.line.x      = element_blank(),
    axis.line.y      = element_blank(),
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(color="grey50", fill=NA, size=1),
    legend.title     = element_text(face="bold", color="black"),
    legend.text      = element_text(color="black"),
    legend.position  = "right"
  )

total_df <- summary_df %>%
  group_by(Experiment, Group) %>%
  summarise(Total=sum(Detections), .groups="drop") %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=c("Chicken_Cloacal", names(Groups)))
  )

Kruskall_Genes_Chicken_Cloacal_Total <- kwAllPairsDunnTest(total_df$Total, total_df$Group, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Genes_Chicken_Cloacal_Total$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Genes_Chicken_Cloacal_Total.xlsx", rowNames = TRUE)

Chicken_Cloacal_Gene_p_total_box <- ggplot(total_df, aes(x=Group, y=Total, fill=Group)) +
  geom_boxplot(outlier.shape=NA, color="black", width=0.6) +
  geom_jitter(aes(color=Group), width=0.2, alpha=0.6) +
  scale_fill_manual(values=group_palette) +
  scale_color_manual(values=group_palette, guide="none") +
  labs(title="Chicken Cloacal Total Detections", y="Detections", x=NULL) +
  theme_minimal() +
  theme(
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.ticks.y    = element_line(color="black"),
    axis.text.y     = element_text(face="bold", color="black"),
    axis.title.y    = element_text(face="bold", color="black"),
    plot.title      = element_text(face="bold", hjust=0, color="black"),
    axis.line.x     = element_blank(),
    axis.line.y     = element_blank(),
    panel.grid      = element_blank(),
    panel.background= element_blank(),
    panel.border    = element_rect(color="grey50", fill=NA, size=1),
    legend.position = "right"
  )

Chicken_Cloacal_Gene_p_total_bar <- ggplot(total_df, aes(x=Group, y=Total, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Detected Viruses (Collapsed Groups)", y="Detections", x=NULL) +
  theme_minimal() +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color="black"),
    axis.text.y      = element_text(face="bold", size=12, color="black"),
    axis.title.y     = element_text(face="bold", size=13, color="black"),
    axis.line.y      = element_line(color="black", linewidth=0.8),
    strip.text       = element_text(size=13, face="bold", color="black"),
    plot.title       = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing    = unit(1,"lines"),
    legend.position  = "none"
  )

print(Chicken_Cloacal_Gene_p_match)
print(Chicken_Cloacal_Gene_p_unique)
print(Chicken_Cloacal_Gene_p_box)
print(Chicken_Cloacal_Gene_p_total_box)
print(Chicken_Cloacal_Gene_p_total_bar)

#================================================
#### DIAMOND Gene Hits - Barplot/Boxplot Summary - Duck Throat
#================================================
DIAMOND  <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames=TRUE)
Metadata <- read.xlsx("Metadata.xlsx",                rowNames=TRUE)

Gene_Count   <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>% distinct() %>%
  mutate(hit=1) %>%
  pivot_wider(names_from=Subject_Seq_ID, values_from=hit, values_fill=0) %>%
  column_to_rownames("Sample_ID")
Binary_Genes <- (Gene_Count>0)*1

DT_samples <- Metadata %>% 
  filter(Location_Grouping=="Duck_Throat") %>% rownames()

Groups <- list(
  Air_Slaughter   = Metadata %>% filter(Location_Grouping=="Air_Slaughter")   %>% rownames(),
  Air_Holding     = Metadata %>% filter(Location_Grouping=="Air_Holding")     %>% rownames(),
  Air_Outside     = Metadata %>% filter(Location_Grouping=="Air_Outside")     %>% rownames(),
  Cage            = Metadata %>% filter(Location_Grouping=="Cage")            %>% rownames(),
  Drinking_Water  = Metadata %>% filter(Location_Grouping=="Drinking_Water")  %>% rownames(),
  Wash_Water      = Metadata %>% filter(Location_Grouping=="Wash_Water")      %>% rownames()
)

DT_mat  <- Binary_Genes[intersect(DT_samples, rownames(Binary_Genes)), , drop=FALSE]
DT_mat  <- DT_mat[, colSums(DT_mat)>0, drop=FALSE]
DT_Taxa <- colnames(DT_mat)

summary_list <- list()
for(exp_id in sort(unique(Metadata$Experiment))){
  exp_meta <- Metadata %>% filter(Experiment==exp_id)
  DT_exp   <- intersect(rownames(exp_meta), DT_samples)
  DT_block <- Binary_Genes[DT_exp, DT_Taxa, drop=FALSE]
  if(nrow(DT_block)==0) next
  DT_cons  <- colSums(DT_block)>0
  
  res <- tibble(
    Experiment = exp_id,
    Group      = "Duck_Throat",
    Type       = "Duck_Throat",
    Detections = sum(DT_cons)
  )
  
  for(g in names(Groups)){
    g_exp   <- intersect(rownames(exp_meta), Groups[[g]])
    g_block <- Binary_Genes[g_exp, DT_Taxa, drop=FALSE]
    if(nrow(g_block)==0) next
    g_cons  <- colSums(g_block)>0
    
    res <- bind_rows(res,
                     tibble(Experiment=exp_id, Group=g,    Type="Match",  Detections=sum(DT_cons & g_cons)),
                     tibble(Experiment=exp_id, Group=g,    Type="Unique", Detections=sum(!DT_cons & g_cons))
    )
  }
  
  summary_list[[as.character(exp_id)]] <- res
}
summary_df <- bind_rows(summary_list)

collected_groups <- Metadata %>%
  group_by(Experiment, Location_Grouping) %>%
  summarise(n=n(), .groups="drop")

long_match <- summary_df %>%
  filter(Type=="Match" | Type=="Duck_Throat") %>%
  inner_join(collected_groups, by=c("Experiment","Group"="Location_Grouping")) %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=c("Duck_Throat", names(Groups)))
  )

long_unique <- summary_df %>%
  filter(Type=="Unique") %>%
  inner_join(collected_groups, by=c("Experiment","Group"="Location_Grouping")) %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=names(Groups))
  )

group_palette <- c(
  "Duck_Throat"     = "#0072B2",
  "Air_Slaughter"   = "#CC79A7",
  "Air_Holding"     = "#6A5ACD",
  "Air_Outside"     = "#228B22",
  "Wash_Water"      = "#F0E442",
  "Drinking_Water"  = "#56B4E9",
  "Cage"            = "#D55E00"
)

Duck_Throat_Gene_p_match <- ggplot(long_match, aes(Group, Detections, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Duck Throat Matches", y="Detections", x=NULL) +
  theme_minimal() +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color="black"),
    axis.text.y      = element_text(face="bold", size=12, color="black"),
    axis.title.y     = element_text(face="bold", size=13, color="black"),
    axis.line.y      = element_line(color="black", linewidth=0.8),
    strip.text       = element_text(size=13, face="bold", color="black"),
    plot.title       = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )

Duck_Throat_Gene_p_unique <- ggplot(long_unique, aes(Group, Detections, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Duck Throat Unique Detections", y="Detections", x=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face="bold", size=12, color="black"),
    axis.title.y     = element_text(face="bold", size=13, color="black"),
    axis.line.y      = element_line(color="black", linewidth=0.8),
    strip.text       = element_text(size=13, face="bold", color="black"),
    plot.title       = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )

box_df <- summary_df %>%
  mutate(
    Combined = paste(Group, Type, sep="_"),
    Combined = factor(Combined, levels=c(
      "Duck_Throat_Duck_Throat",
      "Air_Slaughter_Match",   "Air_Slaughter_Unique",
      "Air_Holding_Match",     "Air_Holding_Unique",
      "Air_Outside_Match",     "Air_Outside_Unique",
      "Wash_Water_Match",      "Wash_Water_Unique",
      "Drinking_Water_Match",  "Drinking_Water_Unique",
      "Cage_Match",            "Cage_Unique"
    ))
  )

box_df$Type[box_df$Type=="Duck_Throat"] <- "Match"
box_plot_df <- box_df %>% filter(Type %in% c("Match","Unique"))

box_plot_df$Group <- factor(
  box_plot_df$Group,
  levels = c("Duck_Throat", names(Groups))
)

Kruskall_Genes_Duck_Throat_Recapture_Unique <- kwAllPairsDunnTest(box_plot_df$Detections, box_plot_df$Combined, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Genes_Duck_Throat_Recapture_Unique$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Genes_Duck_Throat_Recapture_Unique.xlsx", rowNames = TRUE)

Duck_Throat_Gene_p_box <- ggplot(box_plot_df, aes(x=Combined, y=Detections, fill=Group)) +
  geom_boxplot_pattern(
    aes(pattern=Type),
    outlier.shape   = NA,
    pattern_fill    = "black",
    pattern_angle   = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.02,
    color           = "black",
    width           = 0.6,
    show.legend     = TRUE
  ) +
  geom_jitter(aes(color=Group), width=0.2, alpha=0.6, show.legend=FALSE) +
  scale_fill_manual(values=group_palette) +
  scale_color_manual(values=group_palette, guide="none") +
  scale_pattern_manual(
    name   = "Detection Type",
    values = c(Match="none", Unique="stripe"),
    breaks = c("Match","Unique"),
    guide  = guide_legend(
      override.aes = list(
        fill            = c("grey80","grey80"),
        pattern         = c("none","stripe"),
        pattern_fill    = c(NA,"black"),
        pattern_angle   = c(0,45),
        pattern_density = c(0,0.1),
        pattern_spacing = c(0,0.02),
        color           = "black"
      )
    )
  ) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.05))) +
  labs(title="Match vs Unique Boxplots", x=NULL, y="Detections") +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color="black"),
    axis.text.y      = element_text(face="bold", color="black", size=12),
    axis.title.y     = element_text(face="bold", color="black", size=13),
    plot.title       = element_text(face="bold", hjust=0, color="black", size=13),
    axis.line.x      = element_blank(),
    axis.line.y      = element_blank(),
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(color="grey50", fill=NA, size=1),
    legend.title     = element_text(face="bold", color="black"),
    legend.text      = element_text(color="black"),
    legend.position  = "right"
  )

total_df <- summary_df %>%
  group_by(Experiment, Group) %>%
  summarise(Total=sum(Detections), .groups="drop") %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=c("Duck_Throat", names(Groups)))
  )

Kruskall_Genes_Duck_Throat_Total <- kwAllPairsDunnTest(total_df$Total, total_df$Group, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Genes_Duck_Throat_Total$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Genes_Duck_Throat_Total.xlsx", rowNames = TRUE)

Duck_Throat_Gene_p_total_box <- ggplot(total_df, aes(x=Group, y=Total, fill=Group)) +
  geom_boxplot(outlier.shape=NA, color="black", width=0.6) +
  geom_jitter(aes(color=Group), width=0.2, alpha=0.6) +
  scale_fill_manual(values=group_palette) +
  scale_color_manual(values=group_palette, guide="none") +
  labs(title="Duck Throat Total Detections", y="Detections", x=NULL) +
  theme_minimal() +
  theme(
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.ticks.y    = element_line(color="black"),
    axis.text.y     = element_text(face="bold", color="black"),
    axis.title.y    = element_text(face="bold", color="black"),
    plot.title      = element_text(face="bold", hjust=0, color="black"),
    axis.line.x     = element_blank(),
    axis.line.y     = element_blank(),
    panel.grid      = element_blank(),
    panel.background= element_blank(),
    panel.border    = element_rect(color="grey50", fill=NA, size=1),
    legend.position = "right"
  )

Duck_Throat_Gene_p_total_bar <- ggplot(total_df, aes(x=Group, y=Total, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Detected Viruses (Collapsed Groups)", y="Detections", x=NULL) +
  theme_minimal() +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color="black"),
    axis.text.y      = element_text(face="bold", size=12, color="black"),
    axis.title.y     = element_text(face="bold", size=13, color="black"),
    axis.line.y      = element_line(color="black", linewidth=0.8),
    strip.text       = element_text(size=13, face="bold", color="black"),
    plot.title       = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing    = unit(1,"lines"),
    legend.position  = "none"
  )

print(Duck_Throat_Gene_p_match)
print(Duck_Throat_Gene_p_unique)
print(Duck_Throat_Gene_p_box)
print(Duck_Throat_Gene_p_total_box)
print(Duck_Throat_Gene_p_total_bar)

#================================================
#### DIAMOND Gene Hits - Barplot/Boxplot Summary - Duck Cloacal
#================================================
DIAMOND  <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames=TRUE)
Metadata <- read.xlsx("Metadata.xlsx",                rowNames=TRUE)

Gene_Count   <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>% distinct() %>%
  mutate(hit=1) %>%
  pivot_wider(names_from=Subject_Seq_ID, values_from=hit, values_fill=0) %>%
  column_to_rownames("Sample_ID")
Binary_Genes <- (Gene_Count>0)*1

DC_samples <- Metadata %>% 
  filter(Location_Grouping=="Duck_Cloacal") %>% rownames()

Groups <- list(
  Air_Slaughter   = Metadata %>% filter(Location_Grouping=="Air_Slaughter")   %>% rownames(),
  Air_Holding     = Metadata %>% filter(Location_Grouping=="Air_Holding")     %>% rownames(),
  Air_Outside     = Metadata %>% filter(Location_Grouping=="Air_Outside")     %>% rownames(),
  Cage            = Metadata %>% filter(Location_Grouping=="Cage")            %>% rownames(),
  Drinking_Water  = Metadata %>% filter(Location_Grouping=="Drinking_Water")  %>% rownames(),
  Wash_Water      = Metadata %>% filter(Location_Grouping=="Wash_Water")      %>% rownames()
)

DC_mat  <- Binary_Genes[intersect(DC_samples, rownames(Binary_Genes)), , drop=FALSE]
DC_mat  <- DC_mat[, colSums(DC_mat)>0, drop=FALSE]
DC_Taxa <- colnames(DC_mat)

summary_list <- list()
for(exp_id in sort(unique(Metadata$Experiment))){
  exp_meta <- Metadata %>% filter(Experiment==exp_id)
  DC_exp   <- intersect(rownames(exp_meta), DC_samples)
  DC_block <- Binary_Genes[DC_exp, DC_Taxa, drop=FALSE]
  if(nrow(DC_block)==0) next
  DC_cons  <- colSums(DC_block)>0
  
  res <- tibble(
    Experiment = exp_id,
    Group      = "Duck_Cloacal",
    Type       = "Duck_Cloacal",
    Detections = sum(DC_cons)
  )
  
  for(g in names(Groups)){
    g_exp   <- intersect(rownames(exp_meta), Groups[[g]])
    g_block <- Binary_Genes[g_exp, DC_Taxa, drop=FALSE]
    if(nrow(g_block)==0) next
    g_cons  <- colSums(g_block)>0
    
    res <- bind_rows(res,
                     tibble(Experiment=exp_id, Group=g,    Type="Match",  Detections=sum(DC_cons & g_cons)),
                     tibble(Experiment=exp_id, Group=g,    Type="Unique", Detections=sum(!DC_cons & g_cons))
    )
  }
  
  summary_list[[as.character(exp_id)]] <- res
}
summary_df <- bind_rows(summary_list)

collected_groups <- Metadata %>%
  group_by(Experiment, Location_Grouping) %>%
  summarise(n=n(), .groups="drop")

long_match <- summary_df %>%
  filter(Type=="Match" | Type=="Duck_Cloacal") %>%
  inner_join(collected_groups, by=c("Experiment","Group"="Location_Grouping")) %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=c("Duck_Cloacal", names(Groups)))
  )

long_unique <- summary_df %>%
  filter(Type=="Unique") %>%
  inner_join(collected_groups, by=c("Experiment","Group"="Location_Grouping")) %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=names(Groups))
  )

group_palette <- c(
  "Duck_Cloacal"    = "#964B00",
  "Air_Slaughter"   = "#CC79A7",
  "Air_Holding"     = "#6A5ACD",
  "Air_Outside"     = "#228B22",
  "Wash_Water"      = "#F0E442",
  "Drinking_Water"  = "#56B4E9",
  "Cage"            = "#D55E00"
)

Duck_Cloacal_Gene_p_match <- ggplot(long_match, aes(Group, Detections, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Duck Cloacal Matches", y="Detections", x=NULL) +
  theme_minimal() +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color="black"),
    axis.text.y      = element_text(face="bold", size=12, color="black"),
    axis.title.y     = element_text(face="bold", size=13, color="black"),
    axis.line.y      = element_line(color="black", linewidth=0.8),
    strip.text       = element_text(size=13, face="bold", color="black"),
    plot.title       = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )

Duck_Cloacal_Gene_p_unique <- ggplot(long_unique, aes(Group, Detections, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Duck Cloacal Unique Detections", y="Detections", x=NULL) +
  theme_minimal() +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face="bold", size=12, color="black"),
    axis.title.y     = element_text(face="bold", size=13, color="black"),
    axis.line.y      = element_line(color="black", linewidth=0.8),
    strip.text       = element_text(size=13, face="bold", color="black"),
    plot.title       = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )

box_df <- summary_df %>%
  mutate(
    Combined = paste(Group, Type, sep="_"),
    Combined = factor(Combined, levels=c(
      "Duck_Cloacal_Duck_Cloacal",
      "Air_Slaughter_Match",   "Air_Slaughter_Unique",
      "Air_Holding_Match",     "Air_Holding_Unique",
      "Air_Outside_Match",     "Air_Outside_Unique",
      "Wash_Water_Match",      "Wash_Water_Unique",
      "Drinking_Water_Match",  "Drinking_Water_Unique",
      "Cage_Match",            "Cage_Unique"
    ))
  )

box_df$Type[box_df$Type=="Duck_Cloacal"] <- "Match"
box_plot_df <- box_df %>% filter(Type %in% c("Match","Unique"))

box_plot_df$Group <- factor(
  box_plot_df$Group,
  levels = c("Duck_Cloacal", names(Groups))
)

Kruskall_Genes_Duck_Cloacal_Recapture_Unique <- kwAllPairsDunnTest(box_plot_df$Detections, box_plot_df$Combined, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Genes_Duck_Cloacal_Recapture_Unique$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Genes_Duck_Cloacal_Recapture_Unique.xlsx", rowNames = TRUE)

Duck_Cloacal_Gene_p_box <- ggplot(box_plot_df, aes(x=Combined, y=Detections, fill=Group)) +
  geom_boxplot_pattern(
    aes(pattern=Type),
    outlier.shape   = NA,
    pattern_fill    = "black",
    pattern_angle   = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.02,
    color           = "black",
    width           = 0.6,
    show.legend     = TRUE
  ) +
  geom_jitter(aes(color=Group), width=0.2, alpha=0.6, show.legend=FALSE) +
  scale_fill_manual(values=group_palette) +
  scale_color_manual(values=group_palette, guide="none") +
  scale_pattern_manual(
    name   = "Detection Type",
    values = c(Match="none", Unique="stripe"),
    breaks = c("Match","Unique"),
    guide  = guide_legend(
      override.aes = list(
        fill            = c("grey80","grey80"),
        pattern         = c("none","stripe"),
        pattern_fill    = c(NA,"black"),
        pattern_angle   = c(0,45),
        pattern_density = c(0,0.1),
        pattern_spacing = c(0,0.02),
        color           = "black"
      )
    )
  ) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.05))) +
  labs(title="Match vs Unique Boxplots", x=NULL, y="Detections") +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color="black"),
    axis.text.y      = element_text(face="bold", color="black", size=12),
    axis.title.y     = element_text(face="bold", color="black", size=13),
    plot.title       = element_text(face="bold", hjust=0, color="black", size=13),
    axis.line.x      = element_blank(),
    axis.line.y      = element_blank(),
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(color="grey50", fill=NA, size=1),
    legend.title     = element_text(face="bold", color="black"),
    legend.text      = element_text(color="black"),
    legend.position  = "right"
  )

total_df <- summary_df %>%
  group_by(Experiment, Group) %>%
  summarise(Total=sum(Detections), .groups="drop") %>%
  mutate(
    Experiment = factor(Experiment, levels=sort(unique(Experiment))),
    Group      = factor(Group, levels=c("Duck_Cloacal", names(Groups)))
  )

Kruskall_Genes_Duck_Cloacal_Total <- kwAllPairsDunnTest(total_df$Total, total_df$Group, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Genes_Duck_Cloacal_Total$p.value, file ="Supplementary_Files/Supplementary_Kruskall_Genes_Duck_Cloacal_Total.xlsx", rowNames = TRUE)

Duck_Cloacal_Gene_p_total_box <- ggplot(total_df, aes(x=Group, y=Total, fill=Group)) +
  geom_boxplot(outlier.shape=NA, color="black", width=0.6) +
  geom_jitter(aes(color=Group), width=0.2, alpha=0.6) +
  scale_fill_manual(values=group_palette) +
  scale_color_manual(values=group_palette, guide="none") +
  labs(title="Duck Cloacal Total Detections", y="Detections", x=NULL) +
  theme_minimal() +
  theme(
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.ticks.y    = element_line(color="black"),
    axis.text.y     = element_text(face="bold", color="black"),
    axis.title.y    = element_text(face="bold", color="black"),
    plot.title      = element_text(face="bold", hjust=0, color="black"),
    axis.line.x     = element_blank(),
    axis.line.y     = element_blank(),
    panel.grid      = element_blank(),
    panel.background= element_blank(),
    panel.border    = element_rect(color="grey50", fill=NA, size=1),
    legend.position = "right"
  )

Duck_Cloacal_Gene_p_total_bar <- ggplot(total_df, aes(x=Group, y=Total, fill=Group)) +
  geom_bar(stat="identity", color="black") +
  facet_grid(~Experiment, scales="free_x", space="free_x") +
  scale_fill_manual(values=group_palette) +
  labs(title="Detected Viruses (Collapsed Groups)", y="Detections", x=NULL) +
  theme_minimal() +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color="black"),
    axis.text.y      = element_text(face="bold", size=12, color="black"),
    axis.title.y     = element_text(face="bold", size=13, color="black"),
    axis.line.y      = element_line(color="black", linewidth=0.8),
    strip.text       = element_text(size=13, face="bold", color="black"),
    plot.title       = element_text(face="bold", size=13, hjust=0.5, color="black"),
    panel.spacing    = unit(1,"lines"),
    legend.position  = "none"
  )

print(Duck_Cloacal_Gene_p_match)
print(Duck_Cloacal_Gene_p_unique)
print(Duck_Cloacal_Gene_p_box)
print(Duck_Cloacal_Gene_p_total_box)
print(Duck_Cloacal_Gene_p_total_bar)

Recapture_Barplots = Chicken_Throat_Gene_p_match / Chicken_Cloacal_Gene_p_match / Duck_Throat_Gene_p_match / Duck_Cloacal_Gene_p_match
ggsave(filename="Gene_Recapture_Barplots.pdf", plot = Recapture_Barplots, width=200, height=280, units="mm")

Unique_Barplots = Chicken_Throat_Gene_p_unique / Chicken_Cloacal_Gene_p_unique / Duck_Throat_Gene_p_unique / Duck_Cloacal_Gene_p_unique
ggsave(filename="Gene_Unique_Barplots.pdf", plot = Unique_Barplots, width=200, height=280, units="mm")

Total_Barplots = Chicken_Throat_Gene_p_total_bar / Chicken_Cloacal_Gene_p_total_bar / Duck_Throat_Gene_p_total_bar / Duck_Cloacal_Gene_p_total_bar
ggsave(filename="Gene_Total_Barplots.pdf", plot = Total_Barplots, width=200, height=280, units="mm")

Chicken_Recapture_Unique_Boxplots = Chicken_Throat_Gene_p_box /Chicken_Cloacal_Gene_p_box
ggsave(filename="Chicken_Recapture_Unique_Boxplots.pdf", plot = Chicken_Recapture_Unique_Boxplots, width=200, height=280, units="mm")

Duck_Recapture_Unique_Boxplots = Duck_Throat_Gene_p_box /Duck_Cloacal_Gene_p_box
ggsave(filename="Duck_Recapture_Unique_Boxplots.pdf", plot = Duck_Recapture_Unique_Boxplots, width=200, height=280, units="mm")

Chicken_Total_Boxplots = Chicken_Throat_Gene_p_total_box /Chicken_Cloacal_Gene_p_total_box
ggsave(filename="Gene_Chicken_Total_Boxplots.pdf", plot = Chicken_Total_Boxplots, width=200, height=280, units="mm")

Duck_Total_Boxplots = Duck_Throat_Gene_p_total_box /Duck_Cloacal_Gene_p_total_box
ggsave(filename="Gene_Duck_Total_Boxplots.pdf", plot = Duck_Total_Boxplots, width=200, height=280, units="mm")

#================================================
#### DIAMOND Hits - Same Genes - Stacked Barplot - Chicken Throat
#================================================
DIAMOND  <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx",                            rowNames = TRUE)

Gene_Count <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>% distinct() %>%
  mutate(hit = 1) %>%
  pivot_wider(
    names_from  = Subject_Seq_ID,
    values_from = hit,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample_ID")

Binary_Genes <- (Gene_Count > 0) * 1

CT_samples <- rownames(Metadata[Metadata$Location_Grouping == "Chicken_Throat", ])

Groups <- list(
  Air_Slaughter  = rownames(Metadata[Metadata$Location_Grouping=="Air_Slaughter",]),
  Air_Holding    = rownames(Metadata[Metadata$Location_Grouping=="Air_Holding",]),
  Air_Outside    = rownames(Metadata[Metadata$Location_Grouping=="Air_Outside",]),
  Cage           = rownames(Metadata[Metadata$Location_Grouping=="Cage",]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping=="Drinking_Water",]),
  Wash_Water     = rownames(Metadata[Metadata$Location_Grouping=="Wash_Water",])
)

CT_mat  <- Binary_Genes[intersect(CT_samples, rownames(Binary_Genes)), , drop = FALSE]
CT_mat  <- CT_mat[, colSums(CT_mat) > 0, drop = FALSE]
CT_Taxa <- colnames(CT_mat)

match_only_summary <- list()

for (exp_id in sort(unique(Metadata$Experiment))) {
  exp_meta     <- Metadata %>% filter(Experiment == exp_id)
  this_CT      <- intersect(rownames(exp_meta), CT_samples)
  CT_block     <- Binary_Genes[this_CT, CT_Taxa, drop = FALSE]
  if (nrow(CT_block) == 0) next
  
  CT_consensus <- colSums(CT_block) > 0
  ct_taxa      <- names(which(CT_consensus))
  ct_total     <- length(ct_taxa)
  
  results <- c(Chicken_Throat = ct_total)
  
  virus_match_matrix <- matrix(
    0,
    nrow = length(ct_taxa),
    ncol = length(Groups),
    dimnames = list(ct_taxa, names(Groups))
  )
  for (group_name in names(Groups)) {
    grp_samps    <- intersect(rownames(exp_meta), Groups[[group_name]])
    grp_block    <- Binary_Genes[grp_samps, CT_Taxa, drop = FALSE]
    grp_detected <- colSums(grp_block) > 0
    matched_taxa <- ct_taxa[grp_detected[ct_taxa]]
    virus_match_matrix[matched_taxa, group_name] <- 1
  }
  
  one_group_counts <- colSums(
    virus_match_matrix[rowSums(virus_match_matrix) == 1, , drop = FALSE]
  )
  multiple_count <- sum(rowSums(virus_match_matrix) > 1)
  for (g in names(Groups))  results[g] <- one_group_counts[g]
  results["Multiple"] <- multiple_count
  
  match_only_summary[[as.character(exp_id)]] <- results
}

match_summary_df <- do.call(cbind, match_only_summary) %>%
  as.data.frame()
colnames(match_summary_df) <- paste0("Experiment_", colnames(match_summary_df))
match_summary_df$Group <- rownames(match_summary_df)

desired_order <- c("Chicken_Throat", names(Groups), "Multiple")
match_summary_df <- match_summary_df %>%
  filter(Group %in% desired_order) %>%
  mutate(Group = factor(Group, levels = desired_order))

long_df_match <- match_summary_df %>%
  pivot_longer(
    cols      = starts_with("Experiment_"),
    names_to  = "Experiment",
    values_to = "Detections"
  ) %>%
  mutate(
    Experiment_Number = as.numeric(sub("Experiment_", "", Experiment)),
    Group_Label      = if_else(Group == "Chicken_Throat", "Chicken_Throat", "Environmental"),
    Experiment = paste0("Experiment_", Experiment_Number)
  )

long_df_stack <- long_df_match %>%
  filter(Group != "Chicken_Throat")

group_palette <- c(
  "Chicken_Throat" = "#009E73",
  "Air_Slaughter"  = "#CC79A7",
  "Air_Holding"    = "#6A5ACD",
  "Air_Outside"    = "#228B22",
  "Cage"           = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water"     = "#F0E442",
  "Multiple"       = "grey"
)

Chicken_Throat_Stacked <- ggplot() +
  geom_bar(
    data = long_df_match %>% filter(Group == "Chicken_Throat"),
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  geom_bar(
    data = long_df_stack,
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  facet_grid(
    ~ factor(Experiment, levels = unique(long_df_match$Experiment)),
    scales = "free_x", space = "free_x"
  ) +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Detected Genes (Chicken_Throat + Matches)",
    y     = "Number of Detections",
    x     = NULL
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y     = element_text(face = "bold", size = 13, color = "black"),
    axis.line.y      = element_line(color = "black", linewidth = 0.8),
    strip.text       = element_text(face = "bold", size = 13, color = "black"),
    strip.background = element_blank(),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )

print(Chicken_Throat_Stacked)

#================================================
#### DIAMOND Hits - Same Genes - Stacked Barplot - Chicken Cloacal
#================================================
DIAMOND  <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Gene_Count <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>% distinct() %>%
  mutate(hit = 1) %>%
  pivot_wider(
    names_from  = Subject_Seq_ID,
    values_from = hit,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample_ID")

Binary_Genes <- (Gene_Count > 0) * 1

CC_samples <- rownames(Metadata[Metadata$Location_Grouping == "Chicken_Cloacal", ])

Groups <- list(
  Air_Slaughter  = rownames(Metadata[Metadata$Location_Grouping=="Air_Slaughter",]),
  Air_Holding    = rownames(Metadata[Metadata$Location_Grouping=="Air_Holding",]),
  Air_Outside    = rownames(Metadata[Metadata$Location_Grouping=="Air_Outside",]),
  Cage           = rownames(Metadata[Metadata$Location_Grouping=="Cage",]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping=="Drinking_Water",]),
  Wash_Water     = rownames(Metadata[Metadata$Location_Grouping=="Wash_Water",])
)

CC_mat  <- Binary_Genes[intersect(CC_samples, rownames(Binary_Genes)), , drop = FALSE]
CC_mat  <- CC_mat[, colSums(CC_mat) > 0, drop = FALSE]
CC_Taxa <- colnames(CC_mat)

match_only_summary <- list()

for (exp_id in sort(unique(Metadata$Experiment))) {
  exp_meta     <- Metadata %>% filter(Experiment == exp_id)
  this_CC      <- intersect(rownames(exp_meta), CC_samples)
  CC_block     <- Binary_Genes[this_CC, CC_Taxa, drop = FALSE]
  if (nrow(CC_block) == 0) next
  
  CC_consensus <- colSums(CC_block) > 0
  cc_taxa      <- names(which(CC_consensus))
  cc_total     <- length(cc_taxa)
  
  results <- c(Chicken_Cloacal = cc_total)
  
  virus_match_matrix <- matrix(
    0,
    nrow = length(cc_taxa),
    ncol = length(Groups),
    dimnames = list(cc_taxa, names(Groups))
  )
  
  for (group_name in names(Groups)) {
    grp_samps    <- intersect(rownames(exp_meta), Groups[[group_name]])
    grp_block    <- Binary_Genes[grp_samps, CC_Taxa, drop = FALSE]
    grp_detected <- colSums(grp_block) > 0
    
    matched_taxa <- cc_taxa[grp_detected[cc_taxa]]
    virus_match_matrix[matched_taxa, group_name] <- 1
  }
  
  one_group_counts <- colSums(
    virus_match_matrix[rowSums(virus_match_matrix) == 1, , drop = FALSE]
  )
  multiple_count <- sum(rowSums(virus_match_matrix) > 1)
  
  for (g in names(Groups))  results[g] <- one_group_counts[g]
  results["Multiple"] <- multiple_count
  
  match_only_summary[[as.character(exp_id)]] <- results
}

match_summary_df <- do.call(cbind, match_only_summary) %>% as.data.frame()
colnames(match_summary_df) <- paste0("Experiment_", colnames(match_summary_df))
match_summary_df$Group <- rownames(match_summary_df)

desired_order <- c("Chicken_Cloacal", names(Groups), "Multiple")
match_summary_df <- match_summary_df %>%
  filter(Group %in% desired_order) %>%
  mutate(Group = factor(Group, levels = desired_order))

long_df_match <- match_summary_df %>%
  pivot_longer(
    cols      = starts_with("Experiment_"),
    names_to  = "Experiment",
    values_to = "Detections"
  ) %>%
  mutate(
    Experiment_Number = as.numeric(sub("Experiment_", "", Experiment)),
    Group_Label      = if_else(Group == "Chicken_Cloacal", "Chicken_Cloacal", "Environmental"),
    Experiment       = paste0("Experiment_", Experiment_Number)
  )

long_df_stack <- long_df_match %>%
  filter(Group != "Chicken_Cloacal")

group_palette <- c(
  "Chicken_Cloacal" = "#800000",
  "Air_Slaughter"   = "#CC79A7",
  "Air_Holding"     = "#6A5ACD",
  "Air_Outside"     = "#228B22",
  "Cage"            = "#D55E00",
  "Drinking_Water"  = "#56B4E9",
  "Wash_Water"      = "#F0E442",
  "Multiple"        = "grey"
)

Chicken_Cloacal_Stacked <- ggplot() +
  geom_bar(
    data = long_df_match %>% filter(Group == "Chicken_Cloacal"),
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  geom_bar(
    data = long_df_stack,
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  facet_grid(
    ~ factor(Experiment, levels = unique(long_df_match$Experiment)),
    scales = "free_x", space = "free_x"
  ) +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Detected Genes (Chicken_Cloacal + Matches)",
    y     = "Number of Detections",
    x     = NULL
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y     = element_text(face = "bold", size = 13, color = "black"),
    axis.line.y      = element_line(color = "black", linewidth = 0.8),
    strip.text       = element_text(face = "bold", size = 13, color = "black"),
    strip.background = element_blank(),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )

print(Chicken_Cloacal_Stacked)

#================================================
#### DIAMOND Hits - Same Genes - Stacked Barplot - Duck Throat
#================================================
DIAMOND  <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Gene_Count   <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>% distinct() %>%
  mutate(hit = 1) %>%
  pivot_wider(
    names_from  = Subject_Seq_ID,
    values_from = hit,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample_ID")

Binary_Genes <- (Gene_Count > 0) * 1

DT_samples <- rownames(Metadata[Metadata$Location_Grouping == "Duck_Throat", ])

Groups <- list(
  Air_Slaughter  = rownames(Metadata[Metadata$Location_Grouping=="Air_Slaughter",]),
  Air_Holding    = rownames(Metadata[Metadata$Location_Grouping=="Air_Holding",]),
  Air_Outside    = rownames(Metadata[Metadata$Location_Grouping=="Air_Outside",]),
  Cage           = rownames(Metadata[Metadata$Location_Grouping=="Cage",]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping=="Drinking_Water",]),
  Wash_Water     = rownames(Metadata[Metadata$Location_Grouping=="Wash_Water",])
)

DT_mat  <- Binary_Genes[intersect(DT_samples, rownames(Binary_Genes)), , drop = FALSE]
DT_mat  <- DT_mat[, colSums(DT_mat) > 0, drop = FALSE]
DT_Taxa <- colnames(DT_mat)

match_only_summary <- list()

for (exp_id in sort(unique(Metadata$Experiment))) {
  exp_meta     <- Metadata %>% filter(Experiment == exp_id)
  this_DT      <- intersect(rownames(exp_meta), DT_samples)
  DT_block     <- Binary_Genes[this_DT, DT_Taxa, drop = FALSE]
  if (nrow(DT_block) == 0) next
  
  DT_consensus <- colSums(DT_block) > 0
  dt_taxa      <- names(which(DT_consensus))
  dt_total     <- length(dt_taxa)
  
  results <- c(Duck_Throat = dt_total)

  virus_match_matrix <- matrix(
    0,
    nrow = length(dt_taxa),
    ncol = length(Groups),
    dimnames = list(dt_taxa, names(Groups))
  )
  
  for (group_name in names(Groups)) {
    grp_samps    <- intersect(rownames(exp_meta), Groups[[group_name]])
    grp_block    <- Binary_Genes[grp_samps, DT_Taxa, drop = FALSE]
    grp_detected <- colSums(grp_block) > 0
    
    matched_taxa <- dt_taxa[grp_detected[dt_taxa]]
    virus_match_matrix[matched_taxa, group_name] <- 1
  }
  
  one_group_counts <- colSums(
    virus_match_matrix[rowSums(virus_match_matrix) == 1, , drop = FALSE]
  )
  multiple_count <- sum(rowSums(virus_match_matrix) > 1)
  
  for (g in names(Groups))  results[g] <- one_group_counts[g]
  results["Multiple"] <- multiple_count
  
  match_only_summary[[as.character(exp_id)]] <- results
}

match_summary_df <- do.call(cbind, match_only_summary) %>% as.data.frame()
colnames(match_summary_df) <- paste0("Experiment_", colnames(match_summary_df))
match_summary_df$Group <- rownames(match_summary_df)

desired_order <- c("Duck_Throat", names(Groups), "Multiple")
match_summary_df <- match_summary_df %>%
  filter(Group %in% desired_order) %>%
  mutate(Group = factor(Group, levels = desired_order))

long_df_match <- match_summary_df %>%
  pivot_longer(
    cols      = starts_with("Experiment_"),
    names_to  = "Experiment",
    values_to = "Detections"
  ) %>%
  mutate(
    Experiment_Number = as.numeric(sub("Experiment_", "", Experiment)),
    Group_Label      = if_else(Group == "Duck_Throat", "Duck_Throat", "Environmental"),
    Experiment       = paste0("Experiment_", Experiment_Number)
  )

long_df_stack <- long_df_match %>%
  filter(Group != "Duck_Throat")

group_palette <- c(
  "Duck_Throat"    = "#0072B2",
  "Air_Slaughter"  = "#CC79A7",
  "Air_Holding"    = "#6A5ACD",
  "Air_Outside"    = "#228B22",
  "Cage"           = "#D55E00",
  "Drinking_Water" = "#56B4E9",
  "Wash_Water"     = "#F0E442",
  "Multiple"       = "grey"
)

Duck_Throat_Stacked <- ggplot() +
  geom_bar(
    data = long_df_match %>% filter(Group == "Duck_Throat"),
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  geom_bar(
    data = long_df_stack,
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  facet_grid(
    ~ factor(Experiment, levels = unique(long_df_match$Experiment)),
    scales = "free_x", space = "free_x"
  ) +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Detected Genes (Duck_Throat + Matches)",
    y     = "Number of Detections",
    x     = NULL
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y     = element_text(face = "bold", size = 13, color = "black"),
    axis.line.y      = element_line(color = "black", linewidth = 0.8),
    strip.text       = element_text(face = "bold", size = 13, color = "black"),
    strip.background = element_blank(),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )

print(Duck_Throat_Stacked)

#================================================
#### DIAMOND Hits - Same Genes - Stacked Barplot - Duck Cloacal
#================================================
DIAMOND  <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx",                            rowNames = TRUE)

Gene_Count   <- DIAMOND %>%
  select(Sample_ID, Subject_Seq_ID) %>% distinct() %>%
  mutate(hit = 1) %>%
  pivot_wider(
    names_from  = Subject_Seq_ID,
    values_from = hit,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample_ID")

Binary_Genes <- (Gene_Count > 0) * 1

DC_samples <- rownames(Metadata[Metadata$Location_Grouping == "Duck_Cloacal", ])

Groups <- list(
  Air_Slaughter  = rownames(Metadata[Metadata$Location_Grouping=="Air_Slaughter",]),
  Air_Holding    = rownames(Metadata[Metadata$Location_Grouping=="Air_Holding",]),
  Air_Outside    = rownames(Metadata[Metadata$Location_Grouping=="Air_Outside",]),
  Cage           = rownames(Metadata[Metadata$Location_Grouping=="Cage",]),
  Drinking_Water = rownames(Metadata[Metadata$Location_Grouping=="Drinking_Water",]),
  Wash_Water     = rownames(Metadata[Metadata$Location_Grouping=="Wash_Water",])
)

DC_mat  <- Binary_Genes[intersect(DC_samples, rownames(Binary_Genes)), , drop = FALSE]
DC_mat  <- DC_mat[, colSums(DC_mat) > 0, drop = FALSE]
DC_Taxa <- colnames(DC_mat)

match_only_summary <- list()

for (exp_id in sort(unique(Metadata$Experiment))) {
  exp_meta     <- Metadata %>% filter(Experiment == exp_id)
  this_DC      <- intersect(rownames(exp_meta), DC_samples)
  DC_block     <- Binary_Genes[this_DC, DC_Taxa, drop = FALSE]
  if (nrow(DC_block) == 0) next
  
  DC_consensus <- colSums(DC_block) > 0
  dc_taxa      <- names(which(DC_consensus))
  
  results <- c(Duck_Cloacal = length(dc_taxa))
  
  virus_match_matrix <- matrix(
    0,
    nrow = length(dc_taxa),
    ncol = length(Groups),
    dimnames = list(dc_taxa, names(Groups))
  )
  
  for (group_name in names(Groups)) {
    grp_samps    <- intersect(rownames(exp_meta), Groups[[group_name]])
    grp_block    <- Binary_Genes[grp_samps, DC_Taxa, drop = FALSE]
    grp_detected <- colSums(grp_block) > 0
    
    matched_taxa <- dc_taxa[grp_detected[dc_taxa]]
    virus_match_matrix[matched_taxa, group_name] <- 1
  }
  
  one_group_counts <- colSums(
    virus_match_matrix[rowSums(virus_match_matrix) == 1, , drop = FALSE]
  )
  multiple_count <- sum(rowSums(virus_match_matrix) > 1)
  
  for (g in names(Groups))  results[g] <- one_group_counts[g]
  results["Multiple"] <- multiple_count
  
  match_only_summary[[as.character(exp_id)]] <- results
}

match_summary_df <- do.call(cbind, match_only_summary) %>% as.data.frame()
colnames(match_summary_df) <- paste0("Experiment_", colnames(match_summary_df))
match_summary_df$Group <- rownames(match_summary_df)

desired_order <- c("Duck_Cloacal", names(Groups), "Multiple")
match_summary_df <- match_summary_df %>%
  filter(Group %in% desired_order) %>%
  mutate(Group = factor(Group, levels = desired_order))

long_df_match <- match_summary_df %>%
  pivot_longer(
    cols      = starts_with("Experiment_"),
    names_to  = "Experiment",
    values_to = "Detections"
  ) %>%
  mutate(
    Experiment_Number = as.numeric(sub("Experiment_", "", Experiment)),
    Group_Label      = if_else(Group == "Duck_Cloacal", "Duck_Cloacal", "Environmental"),
    Experiment       = paste0("Experiment_", Experiment_Number)
  )

long_df_stack <- long_df_match %>%
  filter(Group != "Duck_Cloacal")

group_palette <- c(
  "Duck_Cloacal"    = "#964B00",
  "Air_Slaughter"   = "#CC79A7",
  "Air_Holding"     = "#6A5ACD",
  "Air_Outside"     = "#228B22",
  "Cage"            = "#D55E00",
  "Drinking_Water"  = "#56B4E9",
  "Wash_Water"      = "#F0E442",
  "Multiple"        = "grey"
)

Duck_Cloacal_Stacked <- ggplot() +
  geom_bar(
    data = long_df_match %>% filter(Group == "Duck_Cloacal"),
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  geom_bar(
    data = long_df_stack,
    aes(x = Group_Label, y = Detections, fill = Group),
    stat = "identity", color = "black", width = 1
  ) +
  facet_grid(
    ~ factor(Experiment, levels = unique(long_df_match$Experiment)),
    scales = "free_x", space = "free_x"
  ) +
  scale_fill_manual(values = group_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Detected Genes (Duck_Cloacal + Matches)",
    y     = "Number of Detections",
    x     = NULL
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face = "bold", size = 12, color = "black"),
    axis.title.y     = element_text(face = "bold", size = 13, color = "black"),
    axis.line.y      = element_line(color = "black", linewidth = 0.8),
    strip.text       = element_text(face = "bold", size = 13, color = "black"),
    strip.background = element_blank(),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5, color = "black"),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "none"
  )

print(Duck_Cloacal_Stacked)

Recapture_Stacked_Barplots = Chicken_Throat_Stacked / Chicken_Cloacal_Stacked / Duck_Throat_Stacked / Duck_Cloacal_Stacked
ggsave(filename="Gene_Recapture_Stacked_Barplots.pdf", plot = Recapture_Stacked_Barplots, width=200, height=280, units="mm")

#================================================
#### DIAMOND Hits - Percent Identity Similar? - Prepare Dataframes
#================================================
DIAMOND <- read.xlsx("DIAMOND/DIAMOND_Final_Data.xlsx", rowNames = TRUE)

DIAMOND_Meta <- DIAMOND[c(1,2,3,4)]
DIAMOND_Meta <- DIAMOND_Meta %>% 
  distinct(Query_Seq_ID, .keep_all = TRUE)
#write.xlsx(DIAMOND_Meta, file ="DIAMOND/Contig_Metadata.xlsx", rowNames = TRUE)

Gene_Meta <- DIAMOND[c(5,10,11,12)]
Gene_Meta <- Gene_Meta %>% 
  distinct(Subject_Seq_ID, .keep_all = TRUE)
#write.xlsx(Gene_Meta, file ="DIAMOND/GenBank_ID_Metadata.xlsx", rowNames = TRUE)

poultry_groups <- c("Chicken_Throat","Chicken_Cloacal","Duck_Throat","Duck_Cloacal")
poultry_ids <- DIAMOND %>%
  filter(Group %in% poultry_groups) %>%
  pull(Subject_Seq_ID) %>%
  unique()

DIAMOND_Meta <- DIAMOND %>%
  select(Query_Seq_ID, Sample_ID, Group, Experiment) %>%
  distinct(Query_Seq_ID, .keep_all = TRUE)

extra_info <- DIAMOND %>%
  select(Query_Seq_ID, Subject_Seq_ID,
         Query_Length_AA,
         Subject_Coverage_Pct,
         Query_Mismatch_Pct) %>%
  distinct()

PID_Matrix <- DIAMOND %>%
  select(Query_Seq_ID, Subject_Seq_ID, Percent_Identity) %>%
  distinct() %>%
  pivot_wider(
    names_from  = Subject_Seq_ID,
    values_from = Percent_Identity,
    values_fill = 0
  ) %>%
  column_to_rownames("Query_Seq_ID") %>%
  select(all_of(poultry_ids))

Coverage_Matrix <- DIAMOND %>%
  select(Query_Seq_ID, Subject_Seq_ID, Subject_Coverage_Pct) %>%
  distinct() %>%
  pivot_wider(
    names_from  = Subject_Seq_ID,
    values_from = Subject_Coverage_Pct,
    values_fill = 0
  ) %>%
  column_to_rownames("Query_Seq_ID") %>%
  select(all_of(poultry_ids))

pair_one_experiment <- function(exp_id) {
  contig_exp <- DIAMOND_Meta %>% filter(Experiment == exp_id)
  
  PID_meta_exp <- PID_Matrix %>%
    rownames_to_column("Query_Seq_ID") %>%
    semi_join(contig_exp, by = "Query_Seq_ID") %>%
    pivot_longer(
      -Query_Seq_ID,
      names_to  = "Subject_Seq_ID",
      values_to = "PID"
    ) %>%
    left_join(contig_exp, by = "Query_Seq_ID") %>%
    left_join(extra_info, by = c("Query_Seq_ID","Subject_Seq_ID"))
  
  poultry_hits <- PID_meta_exp %>%
    filter(Group %in% poultry_groups, PID > 0)
  env_hits <- PID_meta_exp %>%
    filter(!Group %in% poultry_groups, PID > 0)
  
  poultry_tbl <- poultry_hits %>%
    select(
      Subject_Seq_ID,
      Experiment,
      Sample_ID_Poultry    = Sample_ID,
      Group_Poultry        = Group,
      Query_Seq_ID_Poultry = Query_Seq_ID,
      PID_Poultry          = PID,
      Query_Length_AA_Poultry     = Query_Length_AA,
      Subject_Coverage_Pct_Poultry = Subject_Coverage_Pct,
      Query_Mismatch_Pct_Poultry   = Query_Mismatch_Pct
    )
  
  env_tbl <- env_hits %>%
    select(
      Subject_Seq_ID,
      Experiment,
      Sample_ID_Env        = Sample_ID,
      Group_Env            = Group,
      Query_Seq_ID_Env     = Query_Seq_ID,
      PID_Env              = PID,
      Query_Length_AA_Env     = Query_Length_AA,
      Subject_Coverage_Pct_Env = Subject_Coverage_Pct,
      Query_Mismatch_Pct_Env   = Query_Mismatch_Pct
    )
  
  matched <- inner_join(
    poultry_tbl, env_tbl,
    by = c("Subject_Seq_ID", "Experiment")
  )
  
  env_samples <- env_tbl %>% distinct(Sample_ID_Env, Group_Env)
  unmatched_poultry <- anti_join(
    poultry_tbl, env_tbl,
    by = c("Subject_Seq_ID","Experiment")
  ) %>%
    crossing(env_samples) %>%
    mutate(
      Query_Seq_ID_Env           = "Unmatched",
      PID_Env                    = 0,
      Query_Length_AA_Env        = 0,
      Subject_Coverage_Pct_Env   = 0,
      Query_Mismatch_Pct_Env     = 0
    )
  
  poultry_samples <- poultry_tbl %>% distinct(Sample_ID_Poultry, Group_Poultry)
  unmatched_env <- anti_join(
    env_tbl, poultry_tbl,
    by = c("Subject_Seq_ID","Experiment")
  ) %>%
    crossing(poultry_samples) %>%
    mutate(
      Query_Seq_ID_Poultry         = "Unmatched",
      PID_Poultry                  = 0,
      Query_Length_AA_Poultry      = 0,
      Subject_Coverage_Pct_Poultry = 0,
      Query_Mismatch_Pct_Poultry   = 0
    )
  
  bind_rows(matched, unmatched_poultry, unmatched_env) %>%
    select(
      Sample_ID_Poultry, Group_Poultry, Experiment,
      Query_Seq_ID_Poultry, Subject_Seq_ID, PID_Poultry,
      Query_Length_AA_Poultry, Subject_Coverage_Pct_Poultry, Query_Mismatch_Pct_Poultry,
      Sample_ID_Env, Group_Env,
      Query_Seq_ID_Env, PID_Env,
      Query_Length_AA_Env, Subject_Coverage_Pct_Env, Query_Mismatch_Pct_Env
    )
}

experiment_ids <- unique(DIAMOND_Meta$Experiment)
paired_PIDs <- map(set_names(experiment_ids), pair_one_experiment)
Final_PID <- bind_rows(paired_PIDs)

Final_PID <- Final_PID %>%
  mutate(Distance = PID_Poultry - PID_Env)

Final_PID_Coverage <- Final_PID %>%
  mutate(Coverage_Distance = Subject_Coverage_Pct_Poultry - Subject_Coverage_Pct_Env)

Final_PID_Mean <- Final_PID_Coverage %>%
  mutate(Mean_Coverage_PID = (Distance + Coverage_Distance) / 2)

Final_PID <- Final_PID_Mean

Unmatched <- Final_PID %>% 
  filter(if_any(everything(), ~ grepl("Unmatched", .)))
#write.xlsx(Unmatched, "DIAMOND/Unmatched_Contigs.xlsx", rowNames = TRUE)

Recaptured_PID <- Final_PID %>% 
  filter(if_all(everything(), ~ !grepl("Unmatched", .)))
#write.xlsx(Recaptured_PID, "DIAMOND/Recaptured_PID.xlsx", rowNames = TRUE)

Recaptured_PID_Best_Matches <- Recaptured_PID %>% 
  group_by(Query_Seq_ID_Poultry) %>%                      
  mutate(
    keep = case_when(
      any(abs(Distance) <= 0) ~ abs(Distance) <= 0,     
      TRUE                     ~ abs(Distance) == min(abs(Distance))  
    )
  ) %>% 
  ungroup() %>% 
  filter(keep) %>% 
  select(-keep)                                           
#write.xlsx(Recaptured_PID_Best_Matches, "DIAMOND/Recaptured_PID_Best_Matches.xlsx", rowNames = TRUE)

Mean_Recaptured_PID_Gene_Summary <- Recaptured_PID_Best_Matches %>%                
  group_by(Subject_Seq_ID,                
           Group_Poultry,                 
           Group_Env) %>%                 
  summarise(
    meanDistance = mean(Distance, na.rm = TRUE),
    nMatches     = n(),                   
    .groups      = "drop")
#write.xlsx(Mean_Recaptured_PID_Gene_Summary, "DIAMOND/Mean_Recaptured_PID_Gene_Summary.xlsx", rowNames = TRUE)

Mean_Recaptured_PID <- Recaptured_PID %>%                
  group_by(Subject_Seq_ID,                
           Group_Poultry,                 
           Group_Env) %>%                 
  summarise(
    meanDistance = mean(Distance, na.rm = TRUE),
    nMatches     = n(),                   
    .groups      = "drop")
#write.xlsx(Mean_Recaptured_PID, "DIAMOND/Mean_Recaptured_PID.xlsx", rowNames = TRUE)

Mean_PID_Exp <- Final_PID %>%
  group_by(Subject_Seq_ID,
           Sample_ID_Poultry,
           Sample_ID_Env,
           Experiment) %>% 
  summarise(
    Group_Poultry = dplyr::first(Group_Poultry),
    Group_Env     = dplyr::first(Group_Env),
    meanDistance  = mean(Distance, na.rm = TRUE),
    nMatches      = n(),
    .groups       = "drop")
#write.xlsx(Mean_PID_Exp, "DIAMOND/Mean_PID_Experiment.xlsx", rowNames = TRUE)

Recaptured_Coverage_Best_Matches <- Recaptured_PID %>% 
  group_by(Query_Seq_ID_Poultry) %>%                      
  mutate(
    keep = case_when(
      any(abs(Coverage_Distance) <= 0) ~ abs(Coverage_Distance) <= 0,     
      TRUE                     ~ abs(Coverage_Distance) == min(abs(Coverage_Distance))  
    )
  ) %>% 
  ungroup() %>% 
  filter(keep) %>% 
  select(-keep)                                           
#write.xlsx(Recaptured_Coverage_Best_Matches, "DIAMOND/Recaptured_Coverage_Best_Matches.xlsx", rowNames = TRUE)

Mean_Recaptured_Coverage_Gene_Summary <- Recaptured_Coverage_Best_Matches %>%                
  group_by(Subject_Seq_ID,                
           Group_Poultry,                 
           Group_Env) %>%                 
  summarise(
    meanDistance = mean(Coverage_Distance, na.rm = TRUE),
    nMatches     = n(),                   
    .groups      = "drop")
#write.xlsx(Mean_Recaptured_Coverage_Gene_Summary, "DIAMOND/Mean_Recaptured_Coverage_Gene_Summary.xlsx", rowNames = TRUE)

Mean_Recaptured_Coverage <- Recaptured_PID %>%                
  group_by(Subject_Seq_ID,                
           Group_Poultry,                 
           Group_Env) %>%                 
  summarise(
    meanDistance = mean(Coverage_Distance, na.rm = TRUE),
    nMatches     = n(),                   
    .groups      = "drop")
#write.xlsx(Mean_Recaptured_Coverage, "DIAMOND/Mean_Recaptured_Coverage.xlsx", rowNames = TRUE)

Mean_Coverage_Exp <- Final_PID %>%
  group_by(Subject_Seq_ID,
           Sample_ID_Poultry,
           Sample_ID_Env,
           Experiment) %>% 
  summarise(
    Group_Poultry = dplyr::first(Group_Poultry),
    Group_Env     = dplyr::first(Group_Env),
    meanDistance  = mean(Coverage_Distance, na.rm = TRUE),
    nMatches      = n(),
    .groups       = "drop")
#write.xlsx(Mean_PID_Exp, "DIAMOND/Mean_Coverage_Experiment.xlsx", rowNames = TRUE)

#================================================
#### PID Scatter Plots - Agreement Between Seqs Summary
#================================================
Mean_PID <- read.xlsx("DIAMOND/Mean_Recaptured_PID_Gene_Summary.xlsx", rowNames = TRUE)

group_palette <- c(
  "Chicken_Throat"   = "#009E73",
  "Chicken_Cloacal"  = "#800000",
  "Duck_Throat"      = "#0072B2",
  "Duck_Cloacal"     = "#964B00",
  "Air_Slaughter"    = "#CC79A7",
  "Air_Holding"      = "#6A5ACD",
  "Air_Outside"      = "#228B22",
  "Wash_Water"       = "#F0E442",
  "Drinking_Water"   = "#56B4E9",
  "Cage"             = "#D55E00"
)

Mean_PID <- Mean_PID %>%
  mutate(
    distCat = case_when(
      abs(meanDistance) <= 1                       ~ "0–1",
      abs(meanDistance) > 1 & abs(meanDistance) <= 3 ~ "1–3",
      TRUE                                          ~ ">3"
    ),
    distCat = factor(distCat, levels = c("0–1", "1–3", ">3"))
  )

plot_df <- Mean_PID %>%
  mutate(                       
    alphaCat = case_when(
      abs(meanDistance) <= 3  ~ 1.0,
      abs(meanDistance) <= 10 ~ 0.5,
      TRUE                    ~ 0.2
    )
  )

p_scatter <- ggplot(
  plot_df,
  aes(
    x      = Subject_Seq_ID,
    y      = meanDistance,
    size   = nMatches,
    colour = distCat
  )
) +
  geom_point_rast(alpha = 0.5, raster.dpi = 300) +        
  scale_colour_manual(values = c(
    "0–1" = "#E64B35",
    "1–3" = "#F39B7F",
    ">3"  = "#4C72B0"
  )) +
  scale_size_continuous(range = c(1, 4)) +
  scale_y_continuous(limits = c(-30, 30)) +
  labs(
    x     = "Genes",
    y     = "Distance",
    title = "Sequence identity between paired poultry and ES contigs"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "grey50", fill = NA),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(colour = "black"),
    axis.text.y      = element_text(face = "bold"),
    legend.position  = "none"
  )

p_box <- ggplot(
  plot_df,
  aes(
    x      = distCat,
    y      = nMatches,
    colour = distCat
  )
) +
  geom_jitter_rast(width = 0.2, size = 1.6, alpha = 0.6, raster.dpi = 300) +
  scale_colour_manual(
    values = c(
      "0–1" = "#E64B35",
      "1–3" = "#F39B7F",
      ">3"  = "#4C72B0"
    ),
    name = "Distance Group"
  ) +
  scale_y_continuous(position = "left") +
  labs(
    x = NULL,
    y = "Number of Paired Contigs"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.border    = element_rect(colour = "grey50", fill = NA),
    axis.ticks.y    = element_line(colour = "black"),
    axis.text.y     = element_text(face = "bold"),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    legend.position = "none"
  )

bar_table <- table(plot_df$distCat)

bar_df <- data.frame(
  distCat = names(bar_table),
  Count   = as.integer(bar_table),
  stringsAsFactors = FALSE
)

bar_df$distCat <- factor(bar_df$distCat, levels = c("0–1", "1–3", ">3"))

bar <- ggplot(bar_df, aes(x = distCat, y = Count, fill = distCat)) +
  geom_col(color = "black", width = 1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c(
    "0–1" = "#E64B35",
    "1–3" = "#F39B7F",
    ">3"  = "#4C72B0"
  )) +
  labs(
    x   = "Distance Group",
    y   = "Contig Pair Frequency",
    fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(color = "black"),
    axis.text.y      = element_text(face = "bold", color = "black"),
    legend.position  = "right"
  )

combined_fig <- (p_scatter | p_box | bar) +
  plot_layout(widths = c(2, 0.9, 0.3))
print(combined_fig)

ggsave(
  filename = "PID_Distance_Paired_Contigs_Summary.pdf",
  plot     = combined_fig,
  width    = 200,
  height   = 70,
  units    = "mm"
)

#================================================
#### PID Scatter Plots - Agreement Between Seqs Summary Individual ES
#================================================
Mean_PID <- read.xlsx("DIAMOND/Mean_Recaptured_PID_Gene_Summary.xlsx",
                      rowNames = TRUE) %>%
  mutate(
    distCat = case_when(
      abs(meanDistance) <= 1 ~ "0–1",
      abs(meanDistance) <= 3 ~ "1–3",
      TRUE                  ~ ">3"
    ),
    distCat = factor(distCat, levels = c("0–1", "1–3", ">3"))
  )

CT <- subset(Mean_PID, Group_Poultry == "Chicken_Throat")
CC <- subset(Mean_PID, Group_Poultry == "Chicken_Cloacal")
DT <- subset(Mean_PID, Group_Poultry == "Duck_Throat")
DC <- subset(Mean_PID, Group_Poultry == "Duck_Cloacal")

##############################################################################
#  Function to make 3-panel figure for a given environmental group
##############################################################################
make_plots <- function(env_group, mean_tbl) {
  
  plot_df <- mean_tbl %>%
    filter(Group_Env == env_group) %>%
    mutate(
      alphaCat = case_when(
        abs(meanDistance) <=  3 ~ 1.0,
        abs(meanDistance) <= 10 ~ 0.5,
        TRUE                    ~ 0.2
      )
    )
  
  ## ── 1. Scatter plot (raster points) ─────────────────────────────────────
  p_scatter <- ggplot(
    plot_df,
    aes(
      x      = Subject_Seq_ID,
      y      = meanDistance,
      size   = nMatches,
      colour = distCat
    )
  ) +
    geom_point_rast(alpha = 0.5, raster.dpi = 300) +           # ← raster
    scale_colour_manual(values = c(
      "0–1" = "#E64B35",
      "1–3" = "#F39B7F",
      ">3"  = "#4C72B0"
    )) +
    scale_size_continuous(range = c(1, 4)) +
    scale_y_continuous(limits = c(-30, 30)) +
    labs(
      x     = "Genes",
      y     = "Distance",
      title = env_group
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(colour = "grey50", fill = NA),
      panel.grid       = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.ticks.y     = element_line(colour = "black"),
      axis.text.y      = element_text(face = "bold"),
      legend.position  = "none"
    )
  
  ## ── 2. Box / jitter plot (raster jitter points) ─────────────────────────
  p_box <- ggplot(
    plot_df,
    aes(
      x      = distCat,
      y      = nMatches,
      colour = distCat
    )
  ) +
    geom_jitter_rast(width = 0.2,
                     size  = 1.6,
                     alpha = 0.6,
                     raster.dpi = 300) +                       
    scale_colour_manual(values = c(
      "0–1" = "#E64B35",
      "1–3" = "#F39B7F",
      ">3"  = "#4C72B0"
    )) +
    scale_y_continuous(position = "left") +
    labs(
      x = NULL,
      y = "Number of Paired Contigs"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.border    = element_rect(colour = "grey50", fill = NA),
      axis.ticks.y    = element_line(colour = "black"),
      axis.text.y     = element_text(face = "bold"),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      legend.position = "none"
    )
  
  bar_df <- as.data.frame(table(plot_df$distCat))
  names(bar_df) <- c("distCat", "Count")
  bar_df$distCat <- factor(bar_df$distCat, levels = c("0–1", "1–3", ">3"))
  
  bar <- ggplot(bar_df, aes(x = distCat, y = Count, fill = distCat)) +
    geom_col(color = "black", width = 1) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c(
      "0–1" = "#E64B35",
      "1–3" = "#F39B7F",
      ">3"  = "#4C72B0"
    )) +
    labs(
      x   = "Distance Group",
      y   = "Contig Pair Frequency",
      fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_blank(),
      panel.grid       = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.ticks.y     = element_line(color = "black"),
      axis.text.y      = element_text(face = "bold"),
      legend.position  = "right"
    )
  
  (p_scatter | p_box | bar) +
    plot_layout(widths = c(4, 2, 0.7))
}

env_groups <- c("Air_Slaughter", "Air_Holding", "Air_Outside",
                "Wash_Water", "Drinking_Water", "Cage")

make_and_save <- function(poultry_tbl, outfile) {
  fig <- map(env_groups, ~ make_plots(.x, poultry_tbl)) %>%
    set_names(env_groups) %>%
    wrap_plots(ncol = 1) +
    plot_layout(guides = "collect")
  
  print(fig)
  ggsave(
    filename = outfile,
    plot     = fig,
    width    = 200,
    height   = 280,
    units    = "mm"
  )
}

make_and_save(CT, "PID_Distance_Threshold_Scatter_Chicken_Throat.pdf")
make_and_save(CC, "PID_Distance_Threshold_Scatter_Chicken_Cloacal.pdf")
make_and_save(DT, "PID_Distance_Threshold_Scatter_Duck_Throat.pdf")
make_and_save(DC, "PID_Distance_Threshold_Scatter_Duck_Cloacal.pdf")

#================================================
#### PID Boxplot - Summary of Number of Matches
#================================================
Mean_PID_exp <- read.xlsx("DIAMOND/Mean_PID_Experiment.xlsx", rowNames = TRUE)

group_palette <- c(
  "Chicken_Throat"   = "#009E73",
  "Chicken_Cloacal"  = "#800000",
  "Duck_Throat"      = "#0072B2",
  "Duck_Cloacal"     = "#964B00",
  "Air_Slaughter"    = "#CC79A7",
  "Air_Holding"      = "#6A5ACD",
  "Air_Outside"      = "#228B22",
  "Wash_Water"       = "#F0E442",
  "Drinking_Water"   = "#56B4E9",
  "Cage"             = "#D55E00"
)

Mean_PID_exp$Group_Env <- factor(Mean_PID_exp$Group_Env,
                                 levels = names(group_palette))

All_Poultry_jitter <- ggplot(
  Mean_PID_exp,
  aes(x = Group_Env, y = nMatches, colour = Group_Env)
) +
  geom_point_rast(
    position  = position_jitter(width = 0.3, height = 0),
    size      = 2,
    alpha     = 0.6,
    raster.dpi = 300
  ) +
  scale_colour_manual(values = group_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    title = "All Poultry Groups",
    x     = NULL,
    y     = "Number of Paired Contigs"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(colour = "black"),
    axis.text.y      = element_text(face = "bold", colour = "black"),
    axis.title.y     = element_text(face = "bold", colour = "black"),
    plot.title       = element_text(face = "bold", hjust = 0, colour = "black"),
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = "grey50", fill = NA, size = 1),
    legend.position  = "right",
    legend.title     = element_blank()
  )

All_Poultry_jitter
#ggsave(filename="All_Poultry_PID_nMatches.pdf", plot = All_Poultry_jitter, width=150, height=100, units="mm")

All_Poultry <- Mean_PID_exp
All_Poultry$Group_Env <- as.factor(All_Poultry$Group_Env)
Kruskall_All_Poultry_PID_nMatches <- kwAllPairsDunnTest(All_Poultry$nMatches, All_Poultry$Group_Env, p.adjust.method = "fdr")
#write.xlsx(Kruskall_All_Poultry_PID_nMatches$p.value, file ="Supplementary_Files/Kruskall_All_Poultry_PID_nMatches.xlsx", rowNames = TRUE)

#================================================
#### PID Boxplot - Number of Matches Between Each Type
#================================================
Mean_PID_exp <- read.xlsx("DIAMOND/Mean_PID_Experiment.xlsx", rowNames = TRUE)
Metadata      <- read.xlsx("Metadata.xlsx",                  rowNames = TRUE)

group_palette <- c(
  "Chicken_Throat"   = "#009E73",
  "Chicken_Cloacal"  = "#800000",
  "Duck_Throat"      = "#0072B2",
  "Duck_Cloacal"     = "#964B00",
  "Air_Slaughter"    = "#CC79A7",
  "Air_Holding"      = "#6A5ACD",
  "Air_Outside"      = "#228B22",
  "Wash_Water"       = "#F0E442",
  "Drinking_Water"   = "#56B4E9",
  "Cage"             = "#D55E00"
)

poultry_groups <- c("Chicken_Throat", "Chicken_Cloacal",
                    "Duck_Throat",    "Duck_Cloacal")

for (pg in poultry_groups) {
  
  df <- subset(Mean_PID_exp, Group_Poultry == pg)
  df$Group_Env <- factor(df$Group_Env, levels = names(group_palette))
  assign(paste0(pg, "_df"), df, envir = .GlobalEnv)
  
  plt <- ggplot(
    df,
    aes(x = Group_Env, y = nMatches, colour = Group_Env)
  ) +
    geom_point_rast(                                     
      position  = position_jitter(width = 0.2, height = 0),
      size      = 2,
      alpha     = 0.6,
      raster.dpi = 300
    ) +
    scale_colour_manual(values = group_palette, guide = "none") +
    scale_fill_manual(values = group_palette) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      title = gsub("_", " ", pg),
      x     = NULL,
      y     = "Number of Matches"
    ) +
    theme_minimal() +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.ticks.y     = element_line(colour = "black"),
      axis.text.y      = element_text(face = "bold", colour = "black"),
      axis.title.y     = element_text(face = "bold", colour = "black"),
      plot.title       = element_text(face = "bold", hjust = 0, colour = "black"),
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = "grey50", fill = NA, size = 1),
      legend.position  = "right",
      legend.title     = element_blank()
    )
  
  assign(paste0(pg, "_jitter"), plt, envir = .GlobalEnv)
}

Chicken_Throat <- Chicken_Throat_df
Chicken_Throat$Group_Env <- as.factor(Chicken_Throat$Group_Env)
Kruskall_Chicken_Throat_PID_nMatches <- kwAllPairsDunnTest(Chicken_Throat$nMatches, Chicken_Throat$Group_Env, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chicken_Throat_PID_nMatches$p.value, file ="Supplementary_Files/Kruskall_Chicken_Throat_PID_nMatches.xlsx", rowNames = TRUE)

Chicken_Cloacal <- Chicken_Cloacal_df
Chicken_Cloacal$Group_Env <- as.factor(Chicken_Cloacal$Group_Env)
Kruskall_Chicken_Cloacal_PID_nMatches <- kwAllPairsDunnTest(Chicken_Cloacal$nMatches, Chicken_Cloacal$Group_Env, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chicken_Cloacal_PID_nMatches$p.value, file ="Supplementary_Files/Kruskall_Chicken_Cloacal_PID_nMatches.xlsx", rowNames = TRUE)

Duck_Throat <- Duck_Throat_df
Duck_Throat$Group_Env <- as.factor(Duck_Throat$Group_Env)
Kruskall_Duck_Throat_PID_nMatches <- kwAllPairsDunnTest(Duck_Throat$nMatches, Duck_Throat$Group_Env, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Duck_Throat_PID_nMatches$p.value, file ="Supplementary_Files/Kruskall_Duck_Throat_PID_nMatches.xlsx", rowNames = TRUE)

Duck_Cloacal <- Duck_Cloacal_df
Duck_Cloacal$Group_Env <- as.factor(Duck_Cloacal$Group_Env)
Kruskall_Duck_Cloacal_PID_nMatches <- kwAllPairsDunnTest(Duck_Cloacal$nMatches, Duck_Cloacal$Group_Env, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Duck_Cloacal_PID_nMatches$p.value, file ="Supplementary_Files/Kruskall_Duck_Cloacal_PID_nMatches.xlsx", rowNames = TRUE)

Chicken <- Chicken_Throat_jitter / Chicken_Cloacal_jitter
ggsave(filename="Chicken_PID_nMatches.pdf", plot = Chicken, width=200, height=280, units="mm")

Duck <- Duck_Throat_jitter / Duck_Cloacal_jitter
ggsave(filename="Duck_PID_nMatches.pdf", plot = Duck, width=200, height=280, units="mm")

#================================================
#### PID Best Matches - Barplot/Distribution Boxplot Sumamry
#================================================
Recaptured_PID <- read.xlsx("DIAMOND/Recaptured_PID_Best_Matches.xlsx", rowNames = TRUE)

group_palette <- c(
  "Chicken_Throat"   = "#009E73",
  "Chicken_Cloacal"  = "#800000",
  "Duck_Throat"      = "#0072B2",
  "Duck_Cloacal"     = "#964B00",
  "Air_Slaughter"    = "#CC79A7",
  "Air_Holding"      = "#6A5ACD",
  "Air_Outside"      = "#228B22",
  "Wash_Water"       = "#F0E442",
  "Drinking_Water"   = "#56B4E9",
  "Cage"             = "#D55E00")

Recaptured_PID$Group_Env <- factor(Recaptured_PID$Group_Env,
                                   levels = names(group_palette))

p_imperfect <- Recaptured_PID %>% 
  filter(meanDistance < 100) %>% 
  ggplot(aes(y = Group_Env, x = meanDistance,
             fill = Group_Env, colour = Group_Env)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_point_rast(
    position = position_jitter(height = 0.2, width = 0),
    size     = 1,
    alpha    = 0.05,
    raster.dpi = 300
  ) +
  scale_fill_manual(values = group_palette, guide = "none") +
  scale_colour_manual(values = group_palette, guide = "none") +
  scale_x_continuous(limits = c(90, 100),
                     breaks = seq(90, 100, 2)) +
  labs(x = "PID (Env vs Poultry)", y = NULL) +
  theme_minimal() +
  theme(
    axis.text.y         = element_blank(),
    axis.ticks.y        = element_blank(),
    panel.grid.major.y  = element_blank(),
    panel.border        = element_rect(colour = "grey50", fill = NA),
    legend.position     = "right"
  )

Recaptured_PID <- Recaptured_PID %>%                
  group_by(Subject_Seq_ID,                
           Group_Poultry,                 
           Group_Env) %>%                 
  summarise(
    meanDistance = mean(PID_Env, na.rm = TRUE),
    nMatches     = n(),                   
    .groups      = "drop")

p_perfect <- Recaptured_PID %>% 
  mutate(PID_class = if_else(PID_Env == 100, "100 %", "< 100 %")) %>% 
  ggplot(aes(y = Group_Env, fill = PID_class)) +
  geom_bar(position = "fill", colour = "black", width = 0.9) +
  #geom_col(color = "black", width = 1) +
  #scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("100 %"  = "#4C72B0",
                               "< 100 %" = "#D55E00")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Proportion of contigs", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(face = "bold"),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    panel.grid.major.y = element_blank()
  )

Recaptured_PID <- Recaptured_PID %>%
  mutate(PID_class3 = case_when(
    PID_Env == 100                ~ "100 %",
    PID_Env >= 97 & PID_Env < 100 ~ "97-<100 %",
    TRUE                          ~ "<97 %"
  ))

barB_tbl <- Recaptured_PID %>% 
  group_by(Group_Env, PID_class3) %>% 
  summarise(contigs = n_distinct(Query_Seq_ID_Env), .groups = "drop") %>% 
  group_by(Group_Env) %>% 
  mutate(
    prop   = contigs / sum(contigs),
    cumprop= cumsum(prop) - prop/2
  )

pid_cols3 <- c("100 %" = "#E64B35",
               "97-<100 %" = "#F39B7F",
               "<97 %" = "#4C72B0")

p_barB <- ggplot(barB_tbl,
                 aes(y = Group_Env, x = prop, fill = PID_class3)) +
  geom_bar(stat = "identity", width = 0.9, colour = "black") +
  geom_text(aes(x = cumprop, label = sprintf("%.1f%%", 100*prop)),
            colour = "white", size = 3, fontface = "bold") +
  scale_fill_manual(values = pid_cols3) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Proportion of unique contigs", y = NULL,
       fill = NULL, title = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "right")

p_imperfect
#ggsave(filename="PID_Distribution_Summary.pdf", plot = p_imperfect, width=200, height=70, units="mm")

combined_fig <- (p_perfect / p_barB) +
  plot_layout(widths = c(1, 1))
print(combined_fig)
ggsave(filename="Proportion_PID_Contig_Barplots_Summary.pdf", plot = combined_fig, width=300, height=280, units="mm")

#================================================
#### Coverage - Raincloud Plot
#================================================
DIAMOND <- read.xlsx("DIAMOND/Recaptured_Coverage_Best_Matches.xlsx", rowNames = TRUE)

group_palette <- c(
  "Chicken_Throat"   = "#009E73",
  "Chicken_Cloacal"  = "#800000",
  "Duck_Throat"      = "#0072B2",
  "Duck_Cloacal"     = "#964B00",
  "Air_Slaughter"    = "#CC79A7",
  "Air_Holding"      = "#6A5ACD",
  "Air_Outside"      = "#228B22",
  "Wash_Water"       = "#F0E442",
  "Drinking_Water"   = "#56B4E9",
  "Cage"             = "#D55E00"
)

plot_df <- DIAMOND %>% 
  filter(!is.na(Coverage_Distance)) %>% 
  mutate(Group_Env = factor(Group_Env, levels = names(group_palette)))

p_bees <- ggplot(
  plot_df,
  aes(x = Coverage_Distance,
      y = Group_Env,
      colour = Group_Env)
) +
  geom_quasirandom_rast(                 
    width      = 0.4,
    groupOnX   = FALSE,
    alpha      = 0.25,
    size       = 1,
    raster.dpi = 300
  ) +
  scale_colour_manual(values = group_palette, guide = "none") +
  scale_x_continuous(
    limits = c(-80, 80),
    breaks = seq(-80, 80, 20),
    expand = c(0, 0)
  ) +
  labs(
    x = "Coverage Distance",
    y = NULL,
    title = "Coverage-distance points per environmental group (beeswarm)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.border       = element_rect(colour = "grey50", fill = NA),
    panel.grid.major.y = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.y        = element_text(face = "bold"),
    plot.title         = element_text(face = "bold", hjust = 0)
  )

print(p_bees)
#ggsave(filename="Coverage_Distribution_Summary.pdf", plot = p_bees, width=200, height=100, units="mm")

#================================================
#### Contig Taxonomy/Gene Distribution
#================================================
DIAMOND <- read.xlsx("DIAMOND/DIAMOND_PID_Coverage_Annotated.xlsx")
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
Gene <- read.xlsx("DIAMOND/Genes_Metadata_Final_Version.xlsx")
Contig <- read.xlsx("DIAMOND/Contig_Metadata.xlsx", rowNames = TRUE)

DIAMOND <- DIAMOND %>%
  distinct(Query_Seq_ID, .keep_all = TRUE)
DIAMOND <- DIAMOND[c(1,2,13)]
DIAMOND <- DIAMOND %>%
  left_join(Gene %>% select(Subject_Seq_ID, Genus, Species_Taxonomy, Gene_Product),
            by = "Subject_Seq_ID")
DIAMOND <- DIAMOND %>%
  left_join(Contig %>% select(Query_Seq_ID, Sample_ID, Experiment, Group),
            by = "Query_Seq_ID")

group_palette <- c(
  "Chicken_Throat"   = "#009E73",
  "Chicken_Cloacal"  = "#800000",
  "Duck_Throat"      = "#0072B2",
  "Duck_Cloacal"     = "#964B00",
  "Air_Slaughter"    = "#CC79A7",
  "Air_Holding"      = "#6A5ACD",
  "Air_Outside"      = "#228B22",
  "Wash_Water"       = "#F0E442",
  "Drinking_Water"   = "#56B4E9",
  "Cage"             = "#D55E00"
)

genus_pal <- c(
  "Avastrovirus"     = "#8DD3C7",
  "Aviadenovirus"    = "#DCF1B9",
  "Coronaviridae"    = "#E4E3C2",
  "Orthomyxoviridae" = "#C4B4CF",
  "Orthoreovirus"    = "#EE8B86",
  "Orthoretrovirinae"= "#BD98A2",
  "Paramyxoviridae"  = "#F0B36D",
  "Parvovirinae"     = "#D0CD66",
  "Picobirnavirus"   = "#C8D88E",
  "Rotavirus"        = "#FCCDE5",
  "Other"            = "#99B1BC"
)

DIAMOND2 <- DIAMOND %>%
  filter(!is.na(Genus), !is.na(Group)) %>%
  mutate(
    Genus = if_else(Genus %in% names(genus_pal), Genus, "Other"),
    Genus = factor(Genus,
                   levels = c(names(genus_pal)[names(genus_pal) != "Other"], "Other")),
    Group = factor(Group, levels = names(group_palette))
  )

genus_group_counts1 <- DIAMOND2 %>%
  dplyr::count(Group, Genus)

p_stack1 <- ggplot(genus_group_counts1,
                   aes(x = Group, y = n, fill = Genus)) +
  coord_flip() +
  geom_col(width = 0.7, colour = "black") +
  scale_fill_manual(values = genus_pal, name = "Genus") +
  scale_x_discrete(
    limits = names(group_palette),
    labels = gsub("_", " ", names(group_palette)))+
  labs(x = "Sample Group", y = "Number of Contigs", title = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.border    = element_rect(colour = "grey50", fill = NA),
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x    = element_blank(),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black"),
    plot.title      = element_text(face = "bold", hjust = 0),
    legend.position = "bottom"
  )

genus_group_counts2 <- DIAMOND2 %>%
  dplyr::count(Genus, Group) %>%
  mutate(Genus = fct_reorder(Genus, n, .fun = sum, .desc = TRUE))

p_stack2 <- ggplot(genus_group_counts2,
                   aes(x = Genus, y = n, fill = Group)) +
  coord_flip() +
  geom_col(width = 0.7, colour = "black") +
  scale_fill_manual(values = group_palette, name = "Group") +
  labs(x = "Genus", y = "Number of Contigs", title = NULL) +
  #scale_x_discrete(expand = c(0, 0)) +       
  scale_y_continuous(expand = c(0, 0)) +     
  theme_minimal(base_size = 11) +
  theme(
    panel.border    = element_rect(colour = "grey50", fill = NA),
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x    = element_blank(),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black"),
    plot.title      = element_text(face = "bold", hjust = 0),
    legend.position = "bottom"
  )

print(p_stack1)
#ggsave(filename="Matched_Contig_Taxonomy_Distribution.pdf", plot = p_stack1, width=200, height=100, units="mm")

print(p_stack2)
#ggsave(filename="Matched_Contig_Sample_Distribution.pdf", plot = p_stack2, width=200, height=100, units="mm")

p1 <- p_stack1 +
  guides(fill = guide_legend(
    nrow        = 2,           
    byrow       = TRUE,        
    keywidth    = unit(5, "mm"),
    keyheight   = unit(5, "mm"),
    label.position = "right",
    label.hjust     = 0
  )) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.key.size   = unit(3, "mm"),
    legend.spacing.x  = unit(1, "mm"),
    legend.text       = element_text(size = 8),
    legend.title      = element_blank(),
    legend.margin     = margin(t = 0, r = 0, b = 0, l = 0)
  )

print(p1)
#ggsave(filename="Legend_Matched_Contig_Taxonomy_Distribution.pdf", plot = p1, width=200, height=100, units="mm")

#================================================
#### Matched versus Unmatched
#================================================
Unmatched <- read.xlsx("DIAMOND/Unmatched_Contigs.xlsx")
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
Gene <- read.xlsx("DIAMOND/Genes_Metadata_Final_Version.xlsx")

Unmatched <- Unmatched %>%
  left_join(Gene %>% select(Subject_Seq_ID, Genus, Species_Taxonomy, Gene_Product),
            by = "Subject_Seq_ID")

Poultry <- Unmatched[c(2,3,5,21,22,23)]
Poultry <- Poultry %>%
  filter(!if_any(everything(), ~ str_detect(.x, "Unmatched")))
colnames(Poultry)[colnames(Poultry) == "Sample_ID_Poultry"]    <- "Sample_ID"
colnames(Poultry)[colnames(Poultry) == "Group_Poultry"]        <- "Group"
colnames(Poultry)[colnames(Poultry) == "Query_Seq_ID_Poultry"] <- "Query_Seq_ID"
Poultry <- Poultry %>%
  distinct(Query_Seq_ID, .keep_all = TRUE)

Env <- Unmatched[c(11,12,13,21,22,23)]
Env <- Env %>%
  filter(!if_any(everything(), ~ str_detect(.x, "Unmatched")))
colnames(Env)[colnames(Env) == "Sample_ID_Env"]    <- "Sample_ID"
colnames(Env)[colnames(Env) == "Group_Env"]        <- "Group"
colnames(Env)[colnames(Env) == "Query_Seq_ID_Env"] <- "Query_Seq_ID"
Env <- Env %>%
  distinct(Query_Seq_ID, .keep_all = TRUE)

Contigs <- rbind(Poultry, Env)

group_palette <- c(
  "Chicken_Throat"   = "#009E73",
  "Chicken_Cloacal"  = "#800000",
  "Duck_Throat"      = "#0072B2",
  "Duck_Cloacal"     = "#964B00",
  "Air_Slaughter"    = "#CC79A7",
  "Air_Holding"      = "#6A5ACD",
  "Air_Outside"      = "#228B22",
  "Wash_Water"       = "#F0E442",
  "Drinking_Water"   = "#56B4E9",
  "Cage"             = "#D55E00"
)

genus_pal <- c(
  "Avastrovirus"      = "#8DD3C7",
  "Aviadenovirus"     = "#DCF1B9",
  "Coronaviridae"     = "#E4E3C2",
  "Orthomyxoviridae"  = "#C4B4CF",
  "Orthoreovirus"     = "#EE8B86",
  "Orthoretrovirinae" = "#BD98A2",
  "Paramyxoviridae"   = "#F0B36D",
  "Parvovirinae"      = "#D0CD66",
  "Picobirnavirus"    = "#C8D88E",
  "Rotavirus"         = "#FCCDE5",
  "Other"             = "#99B1BC"
)

Contigs_clean <- Contigs %>%
  filter(!is.na(Genus), !is.na(Group)) %>%
  mutate(
    Genus = if_else(Genus %in% names(genus_pal), Genus, "Other"),
    Genus = factor(Genus,
                   levels = c(setdiff(names(genus_pal), "Other"), "Other")),
    Group = factor(Group, levels = names(group_palette))
  )

counts_genus_group <- dplyr::count(Contigs_clean, Genus, Group)
counts_group_genus <- dplyr::count(Contigs_clean, Group, Genus)

counts_group_gene <- Contigs %>%
  filter(!is.na(Gene_Product), !is.na(Group)) %>%
  mutate(
    Gene_Product = fct_lump_n(Gene_Product, n = 19, other_level = "Other"),
    Group        = factor(Group, levels = names(group_palette))
  ) %>%
  dplyr::count(Group, Gene_Product) %>%
  mutate(Gene_Product = fct_reorder(Gene_Product, n, .fun = sum, .desc = TRUE))

gene_levels <- levels(counts_group_gene$Gene_Product)
gene_pal <- setNames(
  ggthemes::tableau_color_pal("Tableau 20")(length(gene_levels)),
  gene_levels
)

base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "grey50", fill = NA),
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black"),
    plot.title   = element_text(face = "bold", hjust = 0),
    legend.position = "bottom",
  )

p1 <- ggplot(counts_genus_group,
             aes(x = Genus, y = n, fill = Group)) +
  coord_flip()+
  geom_col(width = 0.7, colour = "black") +
  scale_fill_manual(values = group_palette) +
  labs(x = "Genus", y = "Contigs", title = NULL) +
  base_theme

p2 <- ggplot(counts_group_genus,
             aes(x = Group, y = n, fill = Genus)) +
  coord_flip()+
  geom_col(width = 0.7, colour = "black") +
  scale_fill_manual(values = genus_pal) +
  scale_x_discrete(
    limits = names(group_palette),
    labels = gsub("_", " ", names(group_palette))
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Group", y = "Contigs", title = NULL) +
  base_theme

p3 <- ggplot(counts_group_gene,
             aes(x = Group, y = n, fill = Gene_Product)) +
  coord_flip()+
  geom_col(width = 0.7, colour = "black") +
  scale_fill_manual(values = gene_pal, name = "Gene Product") +
  scale_x_discrete(
    limits = names(group_palette),
    labels = gsub("_", " ", names(group_palette))
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Group", y = "Contigs", title = NULL) +
  base_theme

p3_stack <- p3 +
  guides(fill = guide_legend(
    nrow        = 4,           
    byrow       = TRUE,        
    keywidth    = unit(5, "mm"),
    keyheight   = unit(5, "mm"),
    label.position = "right",
    label.hjust     = 0
  )) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.key.size   = unit(3, "mm"),
    legend.spacing.x  = unit(1, "mm"),
    legend.text       = element_text(size = 8),
    legend.title      = element_blank(),
    legend.margin     = margin(t = 0, r = 0, b = 0, l = 0)
  )

print(p1)
#ggsave(filename="Unmatched_Gene_Sample_Taxonomy_Sumamry.pdf", plot = p1, width=200, height=100, units="mm")

print(p2)
#ggsave(filename="Unmatched_Gene_Sample_Taxonomy_Breakdown.pdf", plot = p2, width=200, height=100, units="mm")

print(p3_stack)
#ggsave(filename="Unmatched_Gene_Sample_Sumamry.pdf", plot = p3_stack, width=200, height=100, units="mm")
