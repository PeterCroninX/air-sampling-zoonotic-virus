#================================================
#### Packages
#================================================
library(openxlsx)    
library(treeio)      
library(ape)         
library(ggtree)     
library(ggplot2)     
library(dplyr)       
library(patchwork)

#================================================
#### H5 Tree
#================================================
meta <- read.xlsx("Metadata/H5_Tree_Metadata.xlsx")
tf <- read.newick("Final_Tree/H5_Final.tre")

tips_in_tree <- tf$tip.label
meta_filtered <- meta[ meta$Tip_Labels %in% tips_in_tree, ]
meta_filtered <- meta_filtered[ match(tips_in_tree, meta_filtered$Tip_Labels), ]
rownames(meta_filtered) <- meta_filtered$Tip_Labels
length(tips_in_tree)          
nrow(meta_filtered)           
setdiff(tips_in_tree, meta_filtered$Tip_Labels) 

tree_df <- fortify(tf) %>% 
  filter(isTip)
rownames(tree_df) <- tree_df$label
rownames(meta_filtered) <- meta_filtered$Tip_Labels
meta_matched <- meta_filtered[ rownames(tree_df), ]
H5_Metadata <- cbind(tree_df, meta_matched)
Clean_Tree <- ggtree(tf, color = "red", size = 2) %<+% H5_Metadata +
  geom_tiplab(aes(color = LBM %in% drop), size = 3)
H5_TreeData <- as.treedata(Clean_Tree)

H5 <- ggtree(H5_TreeData, size = 0.6) +
  xlim(0, 0.2) +
  geom_treescale() +
  geom_tree(linewidth = 0.6) +
  geom_tippoint(aes(
    color = factor(Clade),
    size  = LBM,
    shape = factor(Type)
  ),
  alpha = 1    
  ) +
  scale_color_discrete(name = "Clade") +
  scale_size_continuous(name = "LBM", range = c(1, 1)) +
  scale_shape_discrete(name = "Type") +
  theme(legend.position = "none")

print(H5)
#ggsave(filename="Unmatched_Gene_Sample_Sumamry.pdf", plot = p3_stack, width=200, height=100, units="mm")

tf <- read.newick("Final_Tree/H5_Clade1.tre")
ggtree(tf)
phy <- if (inherits(tf, "treedata")) tf@phylo else tf
tips <- phy$tip.label
tips_df <- data.frame(Tip_Label = tips, stringsAsFactors = FALSE)
#write.xlsx(tips_df, file= "Metadata/H5_Clade1_Metadata.xlsx")

tf   <- read.tree("Final_Tree/H5_Clade1.tre")
meta <- read.xlsx("Metadata/H5_Clade1_Metadata.xlsx", colNames = TRUE)

tips <- tf$tip.label
meta_filtered <- meta %>% 
  filter(Tip_Label %in% tips) %>% 
  arrange(factor(Tip_Label, levels = tips))
tree_df <- fortify(tf) %>% 
  filter(isTip)

H5_Metadata <- left_join(
  tree_df,
  meta_filtered,
  by = c("label" = "Tip_Label")
)

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

type_levels  <- unique(H5_Metadata$Type)
shape_values <- setNames(
  c(17, 16, 15)[seq_along(type_levels)], 
  type_levels
)

H5_TreeData <- as.treedata(
  ggtree(tf) %<+% H5_Metadata
)

H5_Clade1 <- ggtree(H5_TreeData, size = 1) +
  geom_treescale() +
  geom_tree(linewidth = 1) +
  geom_tippoint(aes(
    color = Group,
    shape = Type
  ), size = 5, alpha = 1) +
  scale_color_manual(name = "Group", values = group_palette) +
  scale_shape_manual(name = "Type", values = shape_values) +
  theme(legend.position = "bottom")

print(H5_Clade1)

heat <- meta[c(2,4)]
rownames(heat) <- heat$Tip_Label
heat$Tip_Label  <- NULL
H5_Clade1 <- gheatmap(
  p      = H5_Clade1,      
  data   = heat, 
  colnames = FALSE, 
  color  = "black",
  offset = 0.02,           
  width  = 0.1
) + 
  hexpand(0.9)
print(H5_Clade1)

tf <- read.newick("Final_Tree/H5_Clade2.tre")
ggtree(tf)
phy <- if (inherits(tf, "treedata")) tf@phylo else tf
tips <- phy$tip.label
tips_df <- data.frame(Tip_Label = tips, stringsAsFactors = FALSE)
#write.xlsx(tips_df, file= "Metadata/H5_Clade2_Metadata.xlsx")

tf   <- read.tree("Final_Tree/H5_Clade2.tre")
meta <- read.xlsx("Metadata/H5_Clade2_Metadata.xlsx", colNames = TRUE)

tips <- tf$tip.label
meta_filtered <- meta %>% 
  filter(Tip_Label %in% tips) %>% 
  arrange(factor(Tip_Label, levels = tips))
tree_df <- fortify(tf) %>% 
  filter(isTip)

H5_Metadata <- left_join(
  tree_df,
  meta_filtered,
  by = c("label" = "Tip_Label")
)

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

type_levels  <- unique(H5_Metadata$Type)
shape_values <- setNames(
  c(17, 16, 15)[seq_along(type_levels)], 
  type_levels
)

H5_TreeData <- as.treedata(
  ggtree(tf) %<+% H5_Metadata
)

H5_Clade2 <- ggtree(H5_TreeData, size = 1) +
  geom_treescale() +
  geom_tree(linewidth = 1) +
  geom_tippoint(aes(
    color = Group,
    shape = Type
  ), size = 5, alpha = 1) +
  scale_color_manual(name = "Group", values = group_palette) +
  scale_shape_manual(name = "Type", values = shape_values) +
  theme(legend.position = "bottom")

print(H5_Clade2)

heat <- meta[c(2,4)]
rownames(heat) <- heat$Tip_Label
heat$Tip_Label  <- NULL
H5_Clade2 <- gheatmap(
  p      = H5_Clade2,      
  data   = heat, 
  colnames = FALSE, 
  color  = "black",
  offset = 0.02,           
  width  = 0.1
) + 
  hexpand(0.9)
print(H5_Clade2)

combined <- H5_Clade1 / H5_Clade2
print(combined)
#ggsave(filename="FiguresH5_Clades.pdf", plot = combined, width=400, height=280, units="mm")

#================================================
#### H6 Tree
#================================================
meta <- read.xlsx("Metadata/H6_Tree_Metadata.xlsx")
tf <- read.newick("Final_Tree/H6_Final.tre")

tips_in_tree <- tf$tip.label
meta_filtered <- meta[ meta$Tip_Labels %in% tips_in_tree, ]
meta_filtered <- meta_filtered[ match(tips_in_tree, meta_filtered$Tip_Labels), ]
rownames(meta_filtered) <- meta_filtered$Tip_Labels
length(tips_in_tree)          
nrow(meta_filtered)           
setdiff(tips_in_tree, meta_filtered$Tip_Labels) 

tree_df <- fortify(tf) %>% 
  filter(isTip)
rownames(tree_df) <- tree_df$label
rownames(meta_filtered) <- meta_filtered$Tip_Labels
meta_matched <- meta_filtered[ rownames(tree_df), ]
H6_Metadata <- cbind(tree_df, meta_matched)
Clean_Tree <- ggtree(tf, color = "red", size = 2) %<+% H6_Metadata +
  geom_tiplab(aes(color = LBM %in% drop), size = 3)
H6_TreeData <- as.treedata(Clean_Tree)

H6 <- ggtree(H6_TreeData, size = 0.6) +
  #xlim(0, 0.2) +
  geom_treescale() +
  geom_tree(linewidth = 0.6) +
  geom_tippoint(aes(
    color = factor(Clade),
    size  = LBM,
    shape = factor(Type)
  ),
  alpha = 1    
  ) +
  scale_color_discrete(name = "Clade") +
  scale_size_continuous(name = "LBM", range = c(1, 1)) +
  scale_shape_discrete(name = "Type") +
  theme(legend.position = "right")

print(H6)
#ggsave(filename="Unmatched_Gene_Sample_Sumamry.pdf", plot = p3_stack, width=200, height=100, units="mm")

tf <- read.newick("Final_Tree/H6_Clade1.tre")
ggtree(tf)
phy <- if (inherits(tf, "treedata")) tf@phylo else tf
tips <- phy$tip.label
tips_df <- data.frame(Tip_Label = tips, stringsAsFactors = FALSE)
#write.xlsx(tips_df, file= "Metadata/H6_Clade1_Metadata.xlsx")

tf   <- read.tree("Final_Tree/H6_Clade1.tre")
meta <- read.xlsx("Metadata/H6_Clade1_Metadata.xlsx", colNames = TRUE)

tips <- tf$tip.label
meta_filtered <- meta %>% 
  filter(Tip_Label %in% tips) %>% 
  arrange(factor(Tip_Label, levels = tips))
tree_df <- fortify(tf) %>% 
  filter(isTip)

H6_Metadata <- left_join(
  tree_df,
  meta_filtered,
  by = c("label" = "Tip_Label")
)

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

type_levels  <- unique(H6_Metadata$Type)
shape_values <- setNames(
  c(17, 16, 15)[seq_along(type_levels)], 
  type_levels
)

H6_TreeData <- as.treedata(
  ggtree(tf) %<+% H6_Metadata
)

H6_Clade1 <- ggtree(H6_TreeData, size = 1) +
  geom_treescale() +
  geom_tree(linewidth = 1) +
  geom_tippoint(aes(
    color = Group,
    shape = Type
  ), size = 5, alpha = 1) +
  scale_color_manual(name = "Group", values = group_palette) +
  scale_shape_manual(name = "Type", values = shape_values) +
  theme(legend.position = "bottom")

print(H6_Clade1)

heat <- meta[c(2,4)]
rownames(heat) <- heat$Tip_Label
heat$Tip_Label  <- NULL
H6_Clade1 <- gheatmap(
  p      = H6_Clade1,      
  data   = heat, 
  colnames = FALSE, 
  color  = "black",
  offset = 0.02,           
  width  = 0.1
) + 
  hexpand(0.9)
print(H6_Clade1)

tf <- read.newick("Final_Tree/H6_Clade2.tre")
ggtree(tf)
phy <- if (inherits(tf, "treedata")) tf@phylo else tf
tips <- phy$tip.label
tips_df <- data.frame(Tip_Label = tips, stringsAsFactors = FALSE)
#write.xlsx(tips_df, file= "Metadata/H6_Clade2_Metadata.xlsx")

tf   <- read.tree("Final_Tree/H6_Clade2.tre")
meta <- read.xlsx("Metadata/H6_Clade2_Metadata.xlsx", colNames = TRUE)

tips <- tf$tip.label
meta_filtered <- meta %>% 
  filter(Tip_Label %in% tips) %>% 
  arrange(factor(Tip_Label, levels = tips))
tree_df <- fortify(tf) %>% 
  filter(isTip)

H6_Metadata <- left_join(
  tree_df,
  meta_filtered,
  by = c("label" = "Tip_Label")
)

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

type_levels  <- unique(H6_Metadata$Type)
shape_values <- setNames(
  c(17, 16, 15)[seq_along(type_levels)], 
  type_levels
)

H6_TreeData <- as.treedata(
  ggtree(tf) %<+% H6_Metadata
)

H6_Clade2 <- ggtree(H6_TreeData, size = 1) +
  geom_treescale() +
  geom_tree(linewidth = 1) +
  geom_tippoint(aes(
    color = Group,
    shape = Type
  ), size = 5, alpha = 1) +
  scale_color_manual(name = "Group", values = group_palette) +
  scale_shape_manual(name = "Type", values = shape_values) +
  theme(legend.position = "bottom")

print(H6_Clade2)

heat <- meta[c(2,4)]
rownames(heat) <- heat$Tip_Label
heat$Tip_Label  <- NULL
H6_Clade2 <- gheatmap(
  p      = H6_Clade2,      
  data   = heat, 
  colnames = FALSE, 
  color  = "black",
  offset = 0.02,           
  width  = 0.1
) + 
  hexpand(0.9)
print(H6_Clade2)

combined <- H6_Clade1 / H6_Clade2
print(combined)
#ggsave(filename="FiguresH6_Clades.pdf", plot = combined, width=400, height=280, units="mm")

#================================================
#### H7 Tree
#================================================
meta <- read.xlsx("Metadata/H7_Tree_Metadata.xlsx")
tf <- read.newick("Final_Tree/H7_Final.tre")

tips_in_tree <- tf$tip.label
meta_filtered <- meta[ meta$Tip_Labels %in% tips_in_tree, ]
meta_filtered <- meta_filtered[ match(tips_in_tree, meta_filtered$Tip_Labels), ]
rownames(meta_filtered) <- meta_filtered$Tip_Labels
length(tips_in_tree)          
nrow(meta_filtered)           
setdiff(tips_in_tree, meta_filtered$Tip_Labels) 

tree_df <- fortify(tf) %>% 
  filter(isTip)
rownames(tree_df) <- tree_df$label
rownames(meta_filtered) <- meta_filtered$Tip_Labels
meta_matched <- meta_filtered[ rownames(tree_df), ]
H7_Metadata <- cbind(tree_df, meta_matched)
Clean_Tree <- ggtree(tf, color = "red", size = 2) %<+% H7_Metadata +
  geom_tiplab(aes(color = LBM %in% drop), size = 3)
H7_TreeData <- as.treedata(Clean_Tree)

H7 <- ggtree(H7_TreeData, size = 0.6) +
  #xlim(0, 0.2) +
  geom_treescale() +
  geom_tree(linewidth = 0.6) +
  geom_tippoint(aes(
    color = factor(Clade),
    size  = LBM,
    shape = factor(Type)
  ),
  alpha = 1    
  ) +
  scale_color_discrete(name = "Clade") +
  scale_size_continuous(name = "LBM", range = c(1.5, 1.5)) +
  scale_shape_discrete(name = "Type") +
  theme(legend.position = "right")

print(H7)
#ggsave(filename="Figures/H7_Final_Tree.pdf", plot = H7, width=200, height=280, units="mm")

#================================================
#### H9 Tree
#================================================
meta <- read.xlsx("Metadata/H9_Tree_Metadata.xlsx")
tf <- read.newick("Final_Tree/H9_Final.tre")

tips_in_tree <- tf$tip.label
meta_filtered <- meta[ meta$Tip_Labels %in% tips_in_tree, ]
meta_filtered <- meta_filtered[ match(tips_in_tree, meta_filtered$Tip_Labels), ]
rownames(meta_filtered) <- meta_filtered$Tip_Labels
length(tips_in_tree)          
nrow(meta_filtered)           
setdiff(tips_in_tree, meta_filtered$Tip_Labels) 

tree_df <- fortify(tf) %>% 
  filter(isTip)
rownames(tree_df) <- tree_df$label
rownames(meta_filtered) <- meta_filtered$Tip_Labels
meta_matched <- meta_filtered[ rownames(tree_df), ]
H9_Metadata <- cbind(tree_df, meta_matched)
Clean_Tree <- ggtree(tf, color = "red", size = 2) %<+% H9_Metadata +
  geom_tiplab(aes(color = LBM %in% drop), size = 3)
H9_TreeData <- as.treedata(Clean_Tree)

H9 <- ggtree(H9_TreeData, size = 0.6) +
  xlim(0, 0.2) +
  geom_treescale() +
  geom_tree(linewidth = 0.6) +
  geom_tippoint(aes(
    color = factor(Clade),
    size  = LBM,
    shape = factor(Type)
  ),
  alpha = 1    
  ) +
  scale_color_discrete(name = "Clade") +
  scale_size_continuous(name = "LBM", range = c(1, 1)) +
  scale_shape_discrete(name = "Type") +
  theme(legend.position = "none")

print(H9)

combined <- H5 | H6 | H9
#ggsave(filename="H5_H6_H9_Trees.pdf", plot = combined, width=200, height=100, units="mm")

tf <- read.newick("Final_Tree/H9_Clade1.tre")
ggtree(tf)
phy <- if (inherits(tf, "treedata")) tf@phylo else tf
tips <- phy$tip.label
tips_df <- data.frame(Tip_Label = tips, stringsAsFactors = FALSE)
#write.xlsx(tips_df, file= "Metadata/H9_Clade1_Metadata.xlsx")

tf   <- read.tree("Final_Tree/H9_Clade1.tre")
meta <- read.xlsx("Metadata/H9_Clade1_Metadata.xlsx", colNames = TRUE)

tips <- tf$tip.label
meta_filtered <- meta %>% 
  filter(Tip_Label %in% tips) %>% 
  arrange(factor(Tip_Label, levels = tips))
tree_df <- fortify(tf) %>% 
  filter(isTip)

H9_Metadata <- left_join(
  tree_df,
  meta_filtered,
  by = c("label" = "Tip_Label")
)

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

type_levels  <- unique(H9_Metadata$Type)
shape_values <- setNames(
  c(17, 16, 15)[seq_along(type_levels)], 
  type_levels
)

H9_TreeData <- as.treedata(
  ggtree(tf) %<+% H9_Metadata
)

H9_Clade1 <- ggtree(H9_TreeData, size = 1) +
  geom_treescale() +
  geom_tree(linewidth = 1) +
  geom_tippoint(aes(
    color = Group,
    shape = Type
  ), size = 5, alpha = 1) +
  scale_color_manual(name = "Group", values = group_palette) +
  scale_shape_manual(name = "Type", values = shape_values) +
  theme(legend.position = "bottom")

print(H9_Clade1)

heat <- meta[c(2,4)]
rownames(heat) <- heat$Tip_Label
heat$Tip_Label  <- NULL
H9_Clade1 <- gheatmap(
  p      = H9_Clade1,      
  data   = heat, 
  colnames = FALSE, 
  color  = "black",
  offset = 0.02,           
  width  = 0.1
) + 
  hexpand(0.9)
print(H9_Clade1)

tf <- read.newick("Final_Tree/H9_Clade2.tre")
ggtree(tf)
phy <- if (inherits(tf, "treedata")) tf@phylo else tf
tips <- phy$tip.label
tips_df <- data.frame(Tip_Label = tips, stringsAsFactors = FALSE)
#write.xlsx(tips_df, file= "Metadata/H9_Clade2_Metadata.xlsx")

tf   <- read.tree("Final_Tree/H9_Clade2.tre")
meta <- read.xlsx("Metadata/H9_Clade2_Metadata.xlsx", colNames = TRUE)

tips <- tf$tip.label
meta_filtered <- meta %>% 
  filter(Tip_Label %in% tips) %>% 
  arrange(factor(Tip_Label, levels = tips))
tree_df <- fortify(tf) %>% 
  filter(isTip)

H9_Metadata <- left_join(
  tree_df,
  meta_filtered,
  by = c("label" = "Tip_Label")
)

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

type_levels  <- unique(H9_Metadata$Type)
shape_values <- setNames(
  c(17, 16, 15)[seq_along(type_levels)], 
  type_levels
)

H9_TreeData <- as.treedata(
  ggtree(tf) %<+% H9_Metadata
)

H9_Clade2 <- ggtree(H9_TreeData, size = 1) +
  geom_treescale() +
  geom_tree(linewidth = 1) +
  geom_tippoint(aes(
    color = Group,
    shape = Type
  ), size = 5, alpha = 1) +
  scale_color_manual(name = "Group", values = group_palette) +
  scale_shape_manual(name = "Type", values = shape_values) +
  theme(legend.position = "bottom")

print(H9_Clade2)

heat <- meta[c(2,4)]
rownames(heat) <- heat$Tip_Label
heat$Tip_Label  <- NULL
H9_Clade2 <- gheatmap(
  p      = H9_Clade2,      
  data   = heat, 
  colnames = FALSE, 
  color  = "black",
  offset = 0.02,           
  width  = 0.1
) + 
  hexpand(0.9)
print(H9_Clade2)

combined <- H9_Clade1 / H9_Clade2
print(combined)
#ggsave(filename="Figures/H9_Clades.pdf", plot = combined, width=400, height=280, units="mm")

#================================================
#### N1 Tree
#================================================
meta <- read.xlsx("Metadata/N1_Tree_Metadata.xlsx")
tf <- read.newick("Final_Tree/N1_Final.tre")

tips_in_tree <- tf$tip.label
meta_filtered <- meta[ meta$Tip_Labels %in% tips_in_tree, ]
meta_filtered <- meta_filtered[ match(tips_in_tree, meta_filtered$Tip_Labels), ]
rownames(meta_filtered) <- meta_filtered$Tip_Labels
length(tips_in_tree)          
nrow(meta_filtered)           
setdiff(tips_in_tree, meta_filtered$Tip_Labels) 

tree_df <- fortify(tf) %>% 
  filter(isTip)
rownames(tree_df) <- tree_df$label
rownames(meta_filtered) <- meta_filtered$Tip_Labels
meta_matched <- meta_filtered[ rownames(tree_df), ]
N1_Metadata <- cbind(tree_df, meta_matched)
Clean_Tree <- ggtree(tf, color = "red", size = 2) %<+% N1_Metadata +
  geom_tiplab(aes(color = LBM %in% drop), size = 3)
N1_TreeData <- as.treedata(Clean_Tree)

N1 <- ggtree(N1_TreeData, size = 0.6) +
  xlim(0, 0.2) +
  geom_treescale() +
  geom_tree(linewidth = 0.6) +
  geom_tippoint(aes(
    color = factor(Clade),
    size  = LBM,
    shape = factor(Type)
  ),
  alpha = 1    
  ) +
  scale_color_discrete(name = "Clade") +
  scale_size_continuous(name = "LBM", range = c(1, 1)) +
  scale_shape_discrete(name = "Type") +
  theme(legend.position = "none")

print(N1)
#ggsave(filename="Unmatched_Gene_Sample_Sumamry.pdf", plot = p3_stack, width=200, height=100, units="mm")

tf <- read.newick("Final_Tree/N1_Clade1.tre")
ggtree(tf)
phy <- if (inherits(tf, "treedata")) tf@phylo else tf
tips <- phy$tip.label
tips_df <- data.frame(Tip_Label = tips, stringsAsFactors = FALSE)
#write.xlsx(tips_df, file= "Metadata/N1_Clade1_Metadata.xlsx")

tf   <- read.tree("Final_Tree/N1_Clade1.tre")
meta <- read.xlsx("Metadata/N1_Clade1_Metadata.xlsx", colNames = TRUE)

tips <- tf$tip.label
meta_filtered <- meta %>% 
  filter(Tip_Label %in% tips) %>% 
  arrange(factor(Tip_Label, levels = tips))
tree_df <- fortify(tf) %>% 
  filter(isTip)

N1_Metadata <- left_join(
  tree_df,
  meta_filtered,
  by = c("label" = "Tip_Label")
)

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

type_levels  <- unique(N1_Metadata$Type)
shape_values <- setNames(
  c(17, 16, 15)[seq_along(type_levels)], 
  type_levels
)

N1_TreeData <- as.treedata(
  ggtree(tf) %<+% N1_Metadata
)

N1_Clade1 <- ggtree(N1_TreeData, size = 1) +
  geom_treescale() +
  geom_tree(linewidth = 1) +
  geom_tippoint(aes(
    color = Group,
    shape = Type
  ), size = 5, alpha = 1) +
  scale_color_manual(name = "Group", values = group_palette) +
  scale_shape_manual(name = "Type", values = shape_values) +
  theme(legend.position = "bottom")

print(N1_Clade1)

heat <- meta[c(2,4)]
rownames(heat) <- heat$Tip_Label
heat$Tip_Label  <- NULL
N1_Clade1 <- gheatmap(
  p      = N1_Clade1,      
  data   = heat, 
  colnames = FALSE, 
  color  = "black",
  offset = 0.02,           
  width  = 0.1
) + 
  hexpand(0.9)
print(N1_Clade1)

tf <- read.newick("Final_Tree/N1_Clade2.tre")
ggtree(tf)
phy <- if (inherits(tf, "treedata")) tf@phylo else tf
tips <- phy$tip.label
tips_df <- data.frame(Tip_Label = tips, stringsAsFactors = FALSE)
#write.xlsx(tips_df, file= "Metadata/N1_Clade2_Metadata.xlsx")

tf   <- read.tree("Final_Tree/N1_Clade2.tre")
meta <- read.xlsx("Metadata/N1_Clade2_Metadata.xlsx", colNames = TRUE)

tips <- tf$tip.label
meta_filtered <- meta %>% 
  filter(Tip_Label %in% tips) %>% 
  arrange(factor(Tip_Label, levels = tips))
tree_df <- fortify(tf) %>% 
  filter(isTip)

N1_Metadata <- left_join(
  tree_df,
  meta_filtered,
  by = c("label" = "Tip_Label")
)

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

type_levels  <- unique(N1_Metadata$Type)
shape_values <- setNames(
  c(17, 16, 15)[seq_along(type_levels)], 
  type_levels
)

N1_TreeData <- as.treedata(
  ggtree(tf) %<+% N1_Metadata
)

N1_Clade2 <- ggtree(N1_TreeData, size = 1) +
  geom_treescale() +
  geom_tree(linewidth = 1) +
  geom_tippoint(aes(
    color = Group,
    shape = Type
  ), size = 5, alpha = 1) +
  scale_color_manual(name = "Group", values = group_palette) +
  scale_shape_manual(name = "Type", values = shape_values) +
  theme(legend.position = "bottom")

print(N1_Clade2)

heat <- meta[c(2,4)]
rownames(heat) <- heat$Tip_Label
heat$Tip_Label  <- NULL
N1_Clade2 <- gheatmap(
  p      = N1_Clade2,      
  data   = heat, 
  colnames = FALSE, 
  color  = "black",
  offset = 0.02,           
  width  = 0.1
) + 
  hexpand(0.9)
print(N1_Clade2)

combined <- N1_Clade1 / N1_Clade2
print(combined)
ggsave(filename="Figures/N1_Clades.pdf", plot = combined, width=400, height=280, units="mm")

#================================================
#### N2 Tree
#================================================
meta <- read.xlsx("Metadata/N2_Tree_Metadata.xlsx")
tf <- read.newick("Final_Tree/N2_Final.tre")

tips_in_tree <- tf$tip.label
meta_filtered <- meta[ meta$Tip_Labels %in% tips_in_tree, ]
meta_filtered <- meta_filtered[ match(tips_in_tree, meta_filtered$Tip_Labels), ]
rownames(meta_filtered) <- meta_filtered$Tip_Labels
length(tips_in_tree)          
nrow(meta_filtered)           
setdiff(tips_in_tree, meta_filtered$Tip_Labels) 

tree_df <- fortify(tf) %>% 
  filter(isTip)
rownames(tree_df) <- tree_df$label
rownames(meta_filtered) <- meta_filtered$Tip_Labels
meta_matched <- meta_filtered[ rownames(tree_df), ]
N2_Metadata <- cbind(tree_df, meta_matched)
Clean_Tree <- ggtree(tf, color = "red", size = 2) %<+% N2_Metadata +
  geom_tiplab(aes(color = LBM %in% drop), size = 3)
N2_TreeData <- as.treedata(Clean_Tree)

N2 <- ggtree(N2_TreeData, size = 0.6) +
  #xlim(0, 0.2) +
  geom_treescale() +
  geom_tree(linewidth = 0.6) +
  geom_tippoint(aes(
    color = factor(Clade),
    size  = LBM,
    shape = factor(Type)
  ),
  alpha = 1    
  ) +
  scale_color_discrete(name = "Clade") +
  scale_size_continuous(name = "LBM", range = c(1, 1)) +
  scale_shape_discrete(name = "Type") +
  theme(legend.position = "none")

print(N2)
#ggsave(filename="Unmatched_Gene_Sample_Sumamry.pdf", plot = p3_stack, width=200, height=100, units="mm")

tf <- read.newick("Final_Tree/N2_Clade1.tre")
ggtree(tf)
phy <- if (inherits(tf, "treedata")) tf@phylo else tf
tips <- phy$tip.label
tips_df <- data.frame(Tip_Label = tips, stringsAsFactors = FALSE)
#write.xlsx(tips_df, file= "Metadata/N2_Clade1_Metadata.xlsx")

tf   <- read.tree("Final_Tree/N2_Clade1.tre")
meta <- read.xlsx("Metadata/N2_Clade1_Metadata.xlsx", colNames = TRUE)

tips <- tf$tip.label
meta_filtered <- meta %>% 
  filter(Tip_Label %in% tips) %>% 
  arrange(factor(Tip_Label, levels = tips))
tree_df <- fortify(tf) %>% 
  filter(isTip)

N2_Metadata <- left_join(
  tree_df,
  meta_filtered,
  by = c("label" = "Tip_Label")
)

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

type_levels  <- unique(N2_Metadata$Type)
shape_values <- setNames(
  c(17, 16, 15)[seq_along(type_levels)], 
  type_levels
)

N2_TreeData <- as.treedata(
  ggtree(tf) %<+% N2_Metadata
)

N2_Clade1 <- ggtree(N2_TreeData, size = 1) +
  geom_treescale() +
  geom_tree(linewidth = 1) +
  geom_tippoint(aes(
    color = Group,
    shape = Type
  ), size = 5, alpha = 1) +
  scale_color_manual(name = "Group", values = group_palette) +
  scale_shape_manual(name = "Type", values = shape_values) +
  theme(legend.position = "bottom")

print(N2_Clade1)

heat <- meta[c(2,4)]
rownames(heat) <- heat$Tip_Label
heat$Tip_Label  <- NULL
N2_Clade1 <- gheatmap(
  p      = N2_Clade1,      
  data   = heat, 
  colnames = FALSE, 
  color  = "black",
  offset = 0.02,           
  width  = 0.1
) + 
  hexpand(0.9)
print(N2_Clade1)

tf <- read.newick("Final_Tree/N2_Clade2.tre")
ggtree(tf)
phy <- if (inherits(tf, "treedata")) tf@phylo else tf
tips <- phy$tip.label
tips_df <- data.frame(Tip_Label = tips, stringsAsFactors = FALSE)
#write.xlsx(tips_df, file= "Metadata/N2_Clade2_Metadata.xlsx")

tf   <- read.tree("Final_Tree/N2_Clade2.tre")
meta <- read.xlsx("Metadata/N2_Clade2_Metadata.xlsx", colNames = TRUE)

tips <- tf$tip.label
meta_filtered <- meta %>% 
  filter(Tip_Label %in% tips) %>% 
  arrange(factor(Tip_Label, levels = tips))
tree_df <- fortify(tf) %>% 
  filter(isTip)

N2_Metadata <- left_join(
  tree_df,
  meta_filtered,
  by = c("label" = "Tip_Label")
)

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

type_levels  <- unique(N2_Metadata$Type)
shape_values <- setNames(
  c(17, 16, 15)[seq_along(type_levels)], 
  type_levels
)

N2_TreeData <- as.treedata(
  ggtree(tf) %<+% N2_Metadata
)

N2_Clade2 <- ggtree(N2_TreeData, size = 1) +
  geom_treescale() +
  geom_tree(linewidth = 1) +
  geom_tippoint(aes(
    color = Group,
    shape = Type
  ), size = 5, alpha = 1) +
  scale_color_manual(name = "Group", values = group_palette) +
  scale_shape_manual(name = "Type", values = shape_values) +
  theme(legend.position = "bottom")

print(N2_Clade2)

heat <- meta[c(2,4)]
rownames(heat) <- heat$Tip_Label
heat$Tip_Label  <- NULL
N2_Clade2 <- gheatmap(
  p      = N2_Clade2,      
  data   = heat, 
  colnames = FALSE, 
  color  = "black",
  offset = 0.02,           
  width  = 0.1
) + 
  hexpand(0.9)
print(N2_Clade2)

combined <- N2_Clade1 / N2_Clade2
print(combined)
#ggsave(filename="FiguresN2_Clades.pdf", plot = combined, width=400, height=280, units="mm")

#================================================
#### N6 Tree
#================================================
meta <- read.xlsx("Metadata/N6_Tree_Metadata.xlsx")
tf <- read.newick("Final_Tree/N6_Final.tre")

tips_in_tree <- tf$tip.label
meta_filtered <- meta[ meta$Tip_Labels %in% tips_in_tree, ]
meta_filtered <- meta_filtered[ match(tips_in_tree, meta_filtered$Tip_Labels), ]
rownames(meta_filtered) <- meta_filtered$Tip_Labels
length(tips_in_tree)          
nrow(meta_filtered)           
setdiff(tips_in_tree, meta_filtered$Tip_Labels) 

tree_df <- fortify(tf) %>% 
  filter(isTip)
rownames(tree_df) <- tree_df$label
rownames(meta_filtered) <- meta_filtered$Tip_Labels
meta_matched <- meta_filtered[ rownames(tree_df), ]
N6_Metadata <- cbind(tree_df, meta_matched)
Clean_Tree <- ggtree(tf, color = "red", size = 2) %<+% N6_Metadata +
  geom_tiplab(aes(color = LBM %in% drop), size = 3)
N6_TreeData <- as.treedata(Clean_Tree)

N6 <- ggtree(N6_TreeData, size = 0.6) +
  #xlim(0, 0.2) +
  geom_treescale() +
  geom_tree(linewidth = 0.6) +
  geom_tippoint(aes(
    color = factor(Clade),
    size  = LBM,
    shape = factor(Type)
  ),
  alpha = 1    
  ) +
  scale_color_discrete(name = "Clade") +
  scale_size_continuous(name = "LBM", range = c(1, 1)) +
  scale_shape_discrete(name = "Type") +
  theme(legend.position = "none")

print(N6)
#ggsave(filename="Unmatched_Gene_Sample_Sumamry.pdf", plot = p3_stack, width=200, height=100, units="mm")

combined <- N1 | N2 | N6
ggsave(filename="N1_N2_N6_Tree.pdf", plot = combined, width=200, height=100, units="mm")

