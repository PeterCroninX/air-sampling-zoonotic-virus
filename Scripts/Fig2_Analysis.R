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
#### Beta-Diversity - Jaccard - Poultry Viruses
#================================================
Jaccard_Li <- read.xlsx("Jaccard.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

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

fill_scale <- scale_fill_manual(values = group_palette)
color_scale <- scale_color_manual(values = group_palette)

calc_centroid_distances <- function(df, group_name) {
  coords <- subset(df, Group == group_name)
  centroid <- data.frame(
    A1 = median(coords$A1),
    A2 = median(coords$A2)
  )
  all_coords <- df[, c("A1", "A2")]
  all_coords <- rbind(centroid, all_coords)
  dists <- as.matrix(vegdist(all_coords, method = "euclidean"))[1, -1]
  return(data.frame(Sample = rownames(df), Distance = dists))
}

centroid_df_list <- list()
groups <- c("Chicken_Throat", "Chicken_Cloacal", "Duck_Throat", "Duck_Cloacal")

for (grp in groups) {
  tmp <- calc_centroid_distances(Jaccard_Li, grp)
  colnames(tmp)[2] <- paste0("Distance_", grp)
  centroid_df_list[[grp]] <- tmp
}

Distance_Median_Centroids <- Reduce(function(x, y) merge(x, y, by = "Sample"), centroid_df_list)
Distance_Median_Centroids <- merge(Distance_Median_Centroids, Metadata, by.x = "Sample", by.y = "row.names")

Distance_Median_Centroids$Location_Grouping <- factor(
  Distance_Median_Centroids$Location_Grouping,
  levels = names(group_palette), ordered = TRUE
)

plot_distance_boxplot <- function(df, group_col, avian_group_label, species_name, y_limit) {
  exclude_groups <- setdiff(c("Chicken_Throat", "Chicken_Cloacal", "Duck_Throat", "Duck_Cloacal"), avian_group_label)
  
  df_filtered <- df[, c("Location_Grouping", group_col)]
  colnames(df_filtered)[2] <- "Distance"
  df_filtered$Location_Grouping <- as.character(df$Location_Grouping)
  df_filtered <- df_filtered[!(df_filtered$Location_Grouping %in% exclude_groups), ]
  
  group_means <- df_filtered %>%
    group_by(Location_Grouping) %>%
    summarise(mean_distance = mean(Distance, na.rm = TRUE)) %>%
    arrange(ifelse(Location_Grouping == avian_group_label, -Inf, mean_distance))
  
  df_filtered$Location_Grouping <- factor(df_filtered$Location_Grouping,
                                          levels = group_means$Location_Grouping)
  
  ggplot(df_filtered, aes(x = Location_Grouping, y = Distance)) +
    geom_boxplot(
      aes(fill = Location_Grouping),
      outlier.shape = NA,
      color = "black",
      width = 0.6
    ) +
    geom_jitter(
      aes(color = Location_Grouping),
      width = 0.2,
      size = 2,
      alpha = 0.6
    ) +
    fill_scale +
    color_scale +
    scale_y_continuous(limits = y_limit, expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      title = species_name,
      subtitle = paste0(gsub("_", " ", sub(".*_", "", group_col)), " Sample"),
      x = NULL,
      y = "Distance From Median Centroid"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.text.y = element_text(face = "bold", color = "black"),
      axis.title.y = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, color = "black", size = 20),
      plot.subtitle = element_text(hjust = 0, face = "bold", size = 14),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, size = 1),
      legend.position = "right"
    )
}

Kruskall_Chicken_Throat <- kwAllPairsDunnTest(Distance_Median_Centroids$Distance_Chicken_Throat, Distance_Median_Centroids$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chicken_Throat$p.value, file ="Supplementary_Files/Supplementary_Median_Centroid_Kruskall_Chicken_Throat.xlsx", rowNames = TRUE)

Kruskall_Chicken_Cloacal <- kwAllPairsDunnTest(Distance_Median_Centroids$Distance_Chicken_Cloacal, Distance_Median_Centroids$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Chicken_Cloacal$p.value, file ="Supplementary_Files/Supplementary_Median_Centroid_Kruskall_Chicken_Cloacal.xlsx", rowNames = TRUE)

Kruskall_Duck_Throat <- kwAllPairsDunnTest(Distance_Median_Centroids$Distance_Duck_Throat, Distance_Median_Centroids$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Duck_Throat$p.value, file ="Supplementary_Files/Supplementary_Median_Centroid_Kruskall_Duck_Throat.xlsx", rowNames = TRUE)

Kruskall_Duck_Cloacal <- kwAllPairsDunnTest(Distance_Median_Centroids$Distance_Duck_Cloacal, Distance_Median_Centroids$Location_Grouping, p.adjust.method = "fdr")
#write.xlsx(Kruskall_Duck_Cloacal$p.value, file ="Supplementary_Files/Supplementary_Median_Centroid_Kruskall_Duck_Cloacal.xlsx", rowNames = TRUE)

Chicken_Throat_Boxplot <- plot_distance_boxplot(Distance_Median_Centroids, "Distance_Chicken_Throat", "Chicken_Throat", "", c(0, 0.5))
Chicken_Cloacal_Boxplot <- plot_distance_boxplot(Distance_Median_Centroids, "Distance_Chicken_Cloacal", "Chicken_Cloacal", "", c(0, 0.5))
Duck_Throat_Boxplot <- plot_distance_boxplot(Distance_Median_Centroids, "Distance_Duck_Throat", "Duck_Throat", "", c(0, 0.8))
Duck_Cloacal_Boxplot <- plot_distance_boxplot(Distance_Median_Centroids, "Distance_Duck_Cloacal", "Duck_Cloacal", "", c(0, 0.65))

Final = Chicken_Throat_Boxplot / Chicken_Cloacal_Boxplot / Duck_Throat_Boxplot / Duck_Cloacal_Boxplot
ggsave(filename="Median_Centroid_Distance_Boxplots.pdf", plot = Final, width=200, height=400, units="mm")

#================================================
#### Beta-Diversity - Spearman Correlations - Centroid
#================================================
Correlations <- read.xlsx("Spearman_Correlations_Virus_Abund_Centroid.xlsx",
                          rowNames = TRUE) |>
  t() |>
  as.data.frame()

Detected   <- read.xlsx("Avian_Species_Detected_List.xlsx",   rowNames = TRUE)
Undetected <- read.xlsx("Avian_Species_Undetected_List.xlsx", rowNames = TRUE)
Detected_Taxa   <- colnames(Detected)
Undetected_Taxa <- colnames(Undetected)
Distance_Correlations <- bind_rows(
  Correlations[, Detected_Taxa]   |> t() |> as.data.frame() |> mutate(Group = "Detected"),
  Correlations[, Undetected_Taxa] |> t() |> as.data.frame() |> mutate(Group = "Undetected")
) |> 
  mutate(Group = factor(Group, levels = c("Detected", "Undetected")))

wilcox.test(Distance_Correlations$Chicken_Throat_Median_Centroid  ~ Distance_Correlations$Group)
wilcox.test(Distance_Correlations$Chicken_Cloacal_Median_Centroid ~ Distance_Correlations$Group)
wilcox.test(Distance_Correlations$Duck_Throat_Median_Centroid     ~ Distance_Correlations$Group)
wilcox.test(Distance_Correlations$Duck_Cloacal_Median_Centroid    ~ Distance_Correlations$Group)

plot_cfg <- list(
  Chicken_Throat_Median_Centroid  = list(title = "Chicken – Throat" , ylim = c(-0.4,  0.4)),
  Chicken_Cloacal_Median_Centroid = list(title = "Chicken – Cloacal", ylim = c(-0.6,  0.4)),
  Duck_Throat_Median_Centroid     = list(title = "Duck – Throat"    , ylim = c(-0.5,  0.75)),
  Duck_Cloacal_Median_Centroid    = list(title = "Duck – Cloacal"   , ylim = c(-0.4,  0.75))
)

palette_det <- c(Detected = "#F8766D", Undetected = "#00BFC4")

make_box <- function(df, metric, cfg) {
  ggplot(df, aes(Group, .data[[metric]], fill = Group, colour = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.4, size = 1, colour = "black") +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = palette_det, name = "") +
    scale_colour_manual(values = palette_det, guide = "none") +
    scale_y_continuous(limits = cfg$ylim) +
    labs(title = cfg$title, x = NULL, y = "Spearman Distance") +
    theme_minimal(base_size = 14) +
    theme(
      panel.border   = element_rect(colour = "grey50", fill = NA, size = 1),
      axis.line.x    = element_line(colour = "grey50", size = 1),
      axis.line.y    = element_line(colour = "grey50", size = 1),
      axis.ticks.x   = element_blank(),
      axis.text.x    = element_blank(),
      axis.text.y    = element_text(face = "bold", colour = "black"),
      axis.title.y   = element_text(face = "bold", colour = "black"),
      plot.title     = element_text(hjust = 0, face = "bold.italic", size = 16),
      legend.position = "bottom",
      legend.text    = element_text(face = "bold", colour = "black")
    )
}

box_plots <- purrr::imap(plot_cfg,
                         ~ make_box(Distance_Correlations, .y, .x))
combined_2x2 <- wrap_plots(box_plots, ncol = 2)
combined_2x2                           
ggsave(filename="Spearman_Distance_Centroid_Virus_Abund_Boxplot.pdf", plot = combined_2x2 , width=200, height=200, units="mm")

#================================================
#### Source Tracking
#================================================
Species <- read.xlsx("Raw_Species_Count.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)

Species <- as.matrix(Species_Count)
Common.Samples <- intersect(rownames(Metadata), rownames(Species))
Species <- Species[Common.Samples,]
Metadata <- Metadata[Common.Samples,]
if(length(Common.Samples) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}
Train <- which(Metadata$Source_Sink=='Source')
Test <- which(Metadata$Source_Sink=='Sink')
Envs <- Metadata$Location_Grouping
if(is.element('Description',colnames(Metadata))) Desc <- Metadata$Location_Grouping

source("/Users/peterc/NUS Dropbox/Peter Brian Cronin/01_Cambodia_Gates_Project/LBM/Rebuttal_March_2025/SourceTracker.r")

Tune.Results <- tune.st(Species[Train,], Envs[Train])
Alpha1 <- Tune.Results$best.alpha1
Alpha2 <- Tune.Results$best.alpha2

ST <- sourcetracker(Species[Train,], Envs[Train])
Results <- predict(ST,Species[Test,], alpha1=Alpha1, alpha2=Alpha2)
Results.Train <- predict(ST, alpha1=Alpha1, alpha2=Alpha2)
Source_Tracking_Results <- as.data.frame(Results$proportions)

ST <- read.xlsx("Source_Tracking_Results.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Source_Tracking_Metadata.xlsx", rowNames = TRUE)

ST <- ST[c(11,1,2,3,4,5)]
ST_Melt <- melt(ST)

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

ST_Ridge_Plot <- ggplot(ST_Melt, aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges(color = "black", linewidth = 0.7, scale = 1) +
  scale_fill_manual(values = group_palette) +
  scale_color_manual(values = group_palette) +
  scale_y_discrete(expand = c(0.01, 0.01)) +
  theme_minimal() +
  theme(
    panel.border     = element_rect(color = "grey50", fill = NA, size = 1),
    panel.background = element_rect(fill = "white"),
    axis.line.x      = element_blank(), 
    axis.line.y      = element_blank(),
    axis.text.x      = element_text(size = 10, face = "bold", color = "black",  
                                    margin = ggplot2::margin(t = 2)),
    axis.text.y      = element_blank(),
    axis.ticks.x     = element_line(color = "grey50", size = 1),
    axis.ticks.y     = element_blank(),
    axis.title.y     = element_blank(),
    axis.title.x     = element_blank(),
    plot.title       = element_text(hjust = 0, size = 14, face = "bold"),
    legend.position  = "right",
    plot.margin      = ggplot2::margin(2, 5, 2, 5))
ST_Ridge_Plot

ggsave(filename="Source_Tracking_Ridge_Plot.pdf", plot = ST_Ridge_Plot, width=120, height=150, units="mm")

#================================================
#### Detected/Undetected Species RPM
#================================================
Species_RPM_Filt <- read.xlsx("Species_RPM_Filt.xlsx", rowNames = TRUE)
Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
Avian_Species <- read.xlsx("Avian_Species_Count.xlsx", rowNames = TRUE)
Avian_Species_Undetected_List <- read.xlsx("Avian_Species_Undetected_List.xlsx", rowNames = TRUE)
Avian_Species_Detected_List <- read.xlsx("Avian_Species_Detected_List.xlsx", rowNames = TRUE)

Avian_RPM <- Species_RPM_Filt[names(Species_RPM_Filt) %in% names(Avian_Species)]
Avian_RPM[Avian_Species == 0] <- 0
Group <- Metadata[c(3)]
Avian_RPM <- cbind(Group, Avian_RPM)
Avian_RPM <- subset(Avian_RPM, Group=="Avian")
Avian_RPM <- Avian_RPM[-c(1)]
Detected_Taxa <- colnames(Avian_Species_Detected_List)
Undetected_Taxa <- colnames(Avian_Species_Undetected_List)
Detected_RPM <- Avian_RPM[Detected_Taxa]
Undetected_RPM <- Avian_RPM[Undetected_Taxa]
Detected_RPM$Group <- "Detected"
Undetected_RPM$Group <- "Undetected"
Detected_Melt <- melt(Detected_RPM, id.vars = "Group")
Undetected_Melt <- melt(Undetected_RPM, id.vars = "Group")
Detection_Rate_RPM <- rbind(Detected_Melt, Undetected_Melt)

Wilcox_Test <- wilcox.test(Detection_Rate_RPM$value ~ Detection_Rate_RPM$Group)

Detection_Rate_RPM$Group<- factor(Detection_Rate_RPM$Group,
                                  levels = c("Detected", "Undetected"),ordered = TRUE)
Detection_Rate_RPM$value <- as.numeric(Detection_Rate_RPM$value)
group_palette <- c(Detected   = "#F8766D",
                   Undetected = "#00BFC4")

Detection_Rate_RPM_Boxplot <- Detection_Rate_RPM %>%                
  ggplot(aes(x = Group, y = value, fill = Group, colour = Group)) +
  geom_boxplot(outlier.shape = NA, colour = "black", width = 0.4) +               
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +                
  scale_fill_manual(values = group_palette) +
  scale_colour_manual(values = group_palette) +
  scale_y_continuous(limits = c(10000, 1000000)) +      
  labs(
    title    = "P<0.0001***",
    subtitle = NULL,
    x        = NULL,
    y        = "RPM"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_blank(),                        
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(colour = "grey50"),
    axis.line.x      = element_line(colour = "grey50"),        
    axis.line.y      = element_line(colour = "grey50"),
    axis.text.y      = element_text(face = "bold", color = "black"),
    axis.title.x     = element_text(face = "bold", color = "black"),
    axis.title.y     = element_text(face = "bold", color = "black"),
    plot.title       = element_text(hjust = 0, face = "bold", color = "black"),
    panel.border     = element_rect(colour = "grey50", fill = NA, linewidth = 1),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(face = "bold", color = "black"))
Detection_Rate_RPM_Boxplot

#================================================
#### Detected/Undetected Species - Times Detected
#================================================
Species <- read.xlsx("Avian_Time_Appeared.xlsx", rowNames = TRUE)

Wilcox_Test <- wilcox.test(Species$Total ~ Species$Group)
Species$Group<- factor(Species$Group,
                       levels = c("Detected", "Undetected"),ordered = TRUE)
Species$Total <- as.numeric(Species$Total)
group_palette <- c(Detected   = "#F8766D",
                   Undetected = "#00BFC4")

Times_Boxplot <- Species %>%                
  ggplot(aes(x = Group, y = Total, fill = Group, colour = Group)) +
  geom_boxplot(outlier.shape = NA, 
               colour = "black", 
               width = 0.4) +                      
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +                
  scale_fill_manual(values = group_palette) +
  scale_colour_manual(values = group_palette) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +      
  labs(
    title    = "P=0.0001***",
    subtitle = NULL,
    x        = NULL,
    y        = "Detection Frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_blank(),    
    axis.text.y      = element_text(face = "bold", color = "black"),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_line(colour = "grey50"),
    axis.line.x      = element_line(colour = "grey50"),                  
    axis.line.y      = element_line(colour = "grey50"),
    axis.title.x     = element_text(face = "bold", color = "black"),
    axis.title.y     = element_text(face = "bold", color = "black"),
    plot.title       = element_text(hjust = 0, face = "bold", color = "black"),
    panel.border     = element_rect(colour = "grey50", fill = NA, linewidth = 1),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(face = "bold", color = "black"))
Times_Boxplot

Combined = Detection_Rate_RPM_Boxplot | Times_Boxplot
ggsave(filename="Detected_Undetected_RPM_FrequencyBoxplots.pdf", plot = Combined, width=200, height=140, units="mm")
