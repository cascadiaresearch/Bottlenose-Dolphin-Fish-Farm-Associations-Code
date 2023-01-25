##### Code for figures and analyses used in Harnish et al. paper on bottlenose associations with a fish farm #####

## code for group size analysis ##

library(tidyverse)
library(stats)
library(ggplot2)

data <- read.csv("Cascadia/Fish Farm/CRCvsIDsgrpsize.csv")

# start with CRC group sizes #
CRCgroups <- filter(data, Source=="CRC")
CRCgroups <- filter(CRCgroups, CRCGrp_Size!="NA")

# test whether the data is normal #
CRCgroupsfarm <- filter(CRCgroups, AtFarm=="Farm")
CRCgroupsnofarm <- filter(CRCgroups, AtFarm=="no")
shapiro.test(CRCgroupsfarm$CRCGrp_Size)
shapiro.test(CRCgroupsnofarm$CRCGrp_Size)

# not all data is normal, so we're going to have to do a Mann-Whitney U Test #
res.mannwhitney <- wilcox.test(CRCGrp_Size ~ AtFarm, data=CRCgroups)
# results have W = 181, p=0.1005 with alternative that location shift is not equal to 0 #
# means that our results are not significant #



## code to create a discovery curve from photo-ID data, using numbers already compiled.## 

library(ggplot2)

disc <- read.csv("Data/AppendixS4_Discoverycurvecounts_27Aug2021.csv")

ggplot() + 
  geom_line(data=disc, color="black", aes(x=Total_IDs, y=Total_uniqueIDs), color="black") +
  theme_classic(base_size = 14) + 
  scale_x_continuous(expand=c(0,0), limits = c(0,105), breaks=c(0,25,50,75,100)) + 
  scale_y_continuous(expand=c(0,0), limits=c(0,105), breaks=c(0,25,50,75,100)) + 
  geom_point(data=disc, color="black", aes(x=Total_IDs, y=Total_uniqueIDs), size=1) + 
  geom_abline(slope=1, color="black", linetype="dashed") + 
  xlab("Cumulative number of identifications") + 
  ylab("Cumulative number of individuals") + 
  theme(axis.text = element_text(colour="black")) + 
  coord_fixed(ratio=1) + 
  geom_text(label="2018", aes(x=22, y=9)) +
  geom_text(label="2019", aes(x=70, y=25)) +
  geom_text(label="2020", aes(x=71, y=32)) +
  geom_text(label="2021", aes(x=100, y=39)) +
  geom_text(label="2010", aes(x=12, y=3)) + 
  geom_text(label="1:1 Reference trendline", aes(x=30, y=86), color="black", size=4) + 
  geom_text(label="Discovery curve of farm associates", aes(x=41, y=93), size=4) +
  geom_text(label="Encounters", aes(x=20, y=100), size=4)+
  geom_point(shape=16, size=2, aes(x=5, y=100)) +
  geom_segment(aes(x=2, xend=7, y=93, yend=93)) +
  geom_segment(aes(x=2, xend=7, y=86, yend=86), linetype="dashed", color="black")
 
  
ggsave("Discoverycurve_labeled_consbio_v2.jpg", width=5, dpi=600)

## code to create maps to view Cascadia Research Collective Hawaii Island effort tracklines ##
## and bottlenose dolphin sightings ##

library(sf)
library(tidyverse)
library(marmap)
library(ggplot2)
library(lubridate)
library(rgdal)
library(ggspatial)
library(metR)
library(raster)
library(sp)


# map themes, labels

label_df <- data.frame(label = c(paste0("Hawai\u02BBi")), lat=c(19.6), lon = c(-155.5))
label_sf <- st_as_sf(label_df, coords = c("lon", "lat"), crs = 4326)

farm <- data.frame(lat=c(19.743), lon=c(-156.061))
farm_sf <- st_as_sf(farm, coords = c("lon", "lat"), crs = 4326)

theme_map <- function() {theme_bw() + 
    theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 1.25), 
          axis.text = element_text(colour = 'black'), 
          plot.title = element_text(colour = 'black', face = 'bold'))}

coast <- st_read("Shapefiles", layer="Coastline")
coastr <- st_transform(coast, crs=4326)
plot(coastr)

# Get effort data

HI_eff <- read.csv("Effort data with on off effort_2021Decv1.csv")

prefarm_HI_eff <- dplyr::filter(HI_eff, Year_ < 2006)
postfarm_HI_eff <- dplyr::filter(HI_eff, Year_ %in% 2006:2021)

prefarm_HI_eff_sf <- sf::st_as_sf(prefarm_HI_eff, coords=c("Long_dd", "Lat_dd"), crs=4326)
postfarm_HI_eff_sf <- sf::st_as_sf(postfarm_HI_eff, coords=c("Long_dd", "Lat_dd"), crs=4326)

prefarm_tracks_sf <- prefarm_HI_eff_sf %>%
  group_by(Date_, Vessel) %>%
  summarise(do_union=F)%>%
  st_cast("LINESTRING")

postfarm_tracks_sf <- postfarm_HI_eff_sf %>%
  group_by(Date_, Vessel) %>%
  summarise(do_union=F)%>%
  st_cast("LINESTRING")

# get bottlenose sightings #

tt_sights <- read.csv("Bottlenose_crcsightings.csv")
tt_sights_sf <- sf::st_as_sf(tt_sights, coords=c("longitude", "latitude"), crs=4326)

prefarm_tt_sights <- dplyr::filter(tt_sights, Year < 2006)
prefarm_tt_sights_sf <- sf::st_as_sf(prefarm_tt_sights, coords=c("longitude", "latitude"), 
                                     crs=4326)

postfarm_tt_sights <- dplyr::filter(tt_sights, Year > 2005)
postfarm_tt_sights_sf <- sf::st_as_sf(postfarm_tt_sights, coords=c("longitude", "latitude"), 
                                      crs=4326)

# farm location #

farm <- data.frame(Long_dd=c(-156.061), Lat_dd=c(19.743))
farm_sf <- st_as_sf(farm, coords=c("Long_dd", "Lat_dd"), crs=4326)

# effort maps #
# pre-farm effort and sightings plot #

ggplot() + 
  geom_sf(data=coastr, fill="white", color="darkgrey", lwd=0.4) + 
  geom_sf(data=prefarm_tracks_sf, lwd=0.4, alpha=0.8, color="red") +
  geom_sf(data=prefarm_tt_sights_sf, shape=16, color="black", size=1.5) +
  geom_sf(data=farm_sf, shape=17, color="yellow", size=2) +
  geom_sf_text(data = label_sf, aes(label = label, fontface = "bold"), size = 4.5) +
  coord_sf(xlim=c(-156.75, -154.75), ylim=c(18.8, 20.3), crs=4326) +
  annotation_north_arrow(location = "tr", which_north="true", 
                         style = north_arrow_fancy_orienteering()) + 
  annotation_scale(location="bl", text_cex=unit(1, "cm")) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), 
        axis.text = element_text(size=11)) + 
  theme_map() + xlab("") + ylab("")

ggsave("HI_effortandsightingsmap_prefarmsightings_v2.jpg")

# post-farm effort and sightings plot #

ggplot() + 
  geom_sf(data=coastr, fill="white", color="darkgrey", lwd=0.4) + 
  geom_sf(data=postfarm_tracks_sf, lwd=0.4, alpha=0.8, color="red") +
  geom_sf(data=postfarm_tt_sights_sf, shape=16, color="black", size=1.5) +
  geom_sf(data=farm_sf, shape=17, color="yellow", size=2) +
  geom_sf_text(data = label_sf, aes(label = label, fontface = "bold"), size = 4.5) +
  coord_sf(xlim=c(-156.75, -154.75), ylim=c(18.8, 20.3), crs=4326) +
  annotation_north_arrow(location = "tr", which_north="true", 
                         style = north_arrow_fancy_orienteering()) + 
  annotation_scale(location="bl", text_cex=unit(1, "cm")) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), 
        axis.text = element_text(size=11)) + 
  theme_map() + xlab("") + ylab("")

ggsave("HI_effortandsightingsmap_postfarmsightings_v2.jpg")

