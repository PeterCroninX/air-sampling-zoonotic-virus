#================================================
#### Packages
#================================================
library(sf)       
library(ggplot2)

#================================================
#### Cambodia Map - Supplementary Figure
#================================================
cambodia_shp <- read_sf("Cambodia.json")
cambodia2_shp <- read_sf("Cambodia_Level2.json")
cambodia_shp$colorSource <- "Other" 
cambodia_shp[cambodia_shp$NAME_1 == "Takêv",]$colorSource <- "Takeo"
cambodia_shp[cambodia_shp$NAME_1 == "PhnomPenh",]$colorSource <- "Phnom Penh"

takeo <- subset(cambodia2_shp, NAME_1 == "Takêv")
takeo$colorSource <- "Other" 
takeo[takeo$NAME_2 == "DounKaev",]$colorSource <- "DounKaev"

pp <- subset(cambodia2_shp, NAME_1 == "PhnomPenh")
pp$colorSource <- "Other" 
pp[pp$NAME_2 == "PrampirMeakkakra",]$colorSource <- "Prampi Makara"

rec_map <- ggplot(cambodia_shp) + 
  geom_sf(aes(fill = colorSource), color = NA) + 
  geom_sf(fill = NA, color = "white", lwd = 1) +
  labs(fill = "Provinces") +
  scale_fill_manual(
    values = c(
      "Phnom Penh" = "#F8766D",  
      "Takeo" = "#00BFC4",  
      "Other" = "grey80"
    )
  ) +
  theme_void()
rec_map

takeo_map <- ggplot(takeo) + 
  geom_sf(aes(fill = colorSource), color = NA) + 
  geom_sf(fill = NA, color = "white", lwd = 1) +
  labs(fill = "Provinces") +
  scale_fill_manual(
    values = c(
      "DounKaev" = "#7CAE00",  
      "Other" = "#00BFC4"
    )
  ) +
  theme_void()
takeo_map

pp_map <- ggplot(pp) + 
  geom_sf(aes(fill = colorSource), color = NA) + 
  geom_sf(fill = NA, color = "white", lwd = 1) +
  labs(fill = "Provinces") +
  scale_fill_manual(
    values = c(
      "Prampi Makara" = "#C77CFF",  
      "Other" = "#F8766D"
    )
  ) +
  theme_void()
pp_map

ggsave(filename="Cambodia_Map.pdf", plot = rec_map, width=200, height=280, units="mm")
ggsave(filename="Takeo_Map.pdf", plot = takeo_map, width=100, height=140, units="mm")
ggsave(filename="Phnom_Penh_Map.pdf", plot = pp_map, width=100, height=140, units="mm")


