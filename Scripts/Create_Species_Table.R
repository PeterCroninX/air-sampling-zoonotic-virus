#================================================
#### Packages
#================================================
library(tidyverse)
library(openxlsx)

#================================================
#### Species/Table Preperation (CZID Output Manipulation)
#================================================
myFiles <- list.files(path = "CZID_Files",
                      pattern = "\\.csv$",
                      full.names = TRUE)

Species <- myFiles %>% 
  set_names() %>% 
  map_dfr(.f = read_delim,
          delim = ";",
          .id = "file_name")

Species_Count <- str_split_fixed(Species$`tax_id,tax_level,genus_tax_id,name,common_name,category,is_phage,agg_score,max_z_score,nt_z_score,nt_rpm,nt_count,nt_contigs,nt_contig_r,nt_percent_identity,nt_alignment_length,nt_e_value,nt_bg_mean,nt_bg_stdev,nt_bg_mean_mass_normalized,nt_bg_stdev_mass_normalized,nr_z_score,nr_rpm,nr_count,nr_contigs,nr_contig_r,nr_percent_identity,nr_alignment_length,nr_e_value,nr_bg_mean,nr_bg_stdev,nr_bg_mean_mass_normalized,nr_bg_stdev_mass_normalized,species_tax_ids,known_pathogen`, ',', 35)
Species_Count <- as.data.frame(Species_Count)
Samples <- str_split_fixed(Species$`file_name`, '/', 9)
Samples <- as.data.frame(Samples)
Samples <- Samples[c(9)]
names(Samples)[names(Samples) == 'V9'] <- 'Sample_Name'
Species_Count <- cbind(Samples, Species_Count)

Species_Type <- Species_Count[c(1,3,5,7,8,12,13,14,15,24,25,26,27)]
Species_List <- subset(Species_Type, V2=="1")
Virus_List <- subset(Species_List, V6=="viruses")
Clean_Virus_List <- subset(Virus_List, V7=="false")
NT_Virus <- Clean_Virus_List[c(1,3,6,7,8,9)]
names(NT_Virus)[names(NT_Virus) == 'V4'] <- 'Species'
names(NT_Virus)[names(NT_Virus) == 'V11'] <- 'RPM'
names(NT_Virus)[names(NT_Virus) == 'V12'] <- 'Read_Count'
names(NT_Virus)[names(NT_Virus) == 'V13'] <- 'Contig_Number'
names(NT_Virus)[names(NT_Virus) == 'V14'] <- 'Contig_Read_Count'
Final_NT_Virus <- NT_Virus[!(is.na(NT_Virus$Contig_Number) | NT_Virus$Contig_Number==""), ]
NR_Virus <- Clean_Virus_List[c(1,3,10,11,12,13)]
names(NR_Virus)[names(NR_Virus) == 'V4'] <- 'Species'
names(NR_Virus)[names(NR_Virus) == 'V23'] <- 'RPM'
names(NR_Virus)[names(NR_Virus) == 'V24'] <- 'Read_Count'
names(NR_Virus)[names(NR_Virus) == 'V25'] <- 'Contig_Number'
names(NR_Virus)[names(NR_Virus) == 'V26'] <- 'Contig_Read_Count'
Final_NR_Virus <- NR_Virus[!(is.na(NR_Virus$Contig_Number) | NR_Virus$Contig_Number==""), ]
Combined_Virus <- rbind(Final_NT_Virus, Final_NR_Virus)
Combined_Virus$Names <- paste(Combined_Virus$Sample_Name, "_", Combined_Virus$Species)
Combined_Virus$Contig_Number <- as.numeric(Combined_Virus$Contig_Number)
Ordered_Virus <- Combined_Virus[order(Combined_Virus$Names, -abs(Combined_Virus$Contig_Number) ), ]
Clean_Combined_Virus <- Ordered_Virus[ !duplicated(Ordered_Virus$Names), ]
Virus_Count <- Clean_Combined_Virus[c(1,2,6)]
Virus_Count$Contig_Read_Count <- as.numeric(Virus_Count$Contig_Read_Count)
Virus_Count <- aggregate(Contig_Read_Count ~ Species + Sample_Name, Virus_Count, sum)
Virus_Count <- reshape(Virus_Count, idvar = "Sample_Name", timevar = "Species", direction = "wide")
Virus_Count[is.na(Virus_Count)] = 0
#write.xlsx(Virus_Count, file = "Negative_Control_Virus_Count.xlsx", rowNames = TRUE)
