library(rnaturalearth)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(sf)
library(mapview)
library(dplyr)
library(terra)
setwd("~/Desktop/BANS_adaptive_units/")

# Read pops and meta
pops <- read.table("analysis/03.GEA/BANS_sample_list_pop.tsv", sep="\t", header=FALSE) %>%
  rename(BGP_ID = V1, Pop = V3) %>%
  dplyr::select(BGP_ID, Pop)

meta <- read.delim("data/spatial_files/BANS.initial.worldclim.txt") %>%
  dplyr::select(BGP_ID, Long, Lat)

pops <- left_join(meta, pops, by="BGP_ID")

# Define cluster groups
West_clusters    <- c("Cluster_1", "Cluster_5")
Midwest_clusters <- c("Cluster_4", "Cluster_8", "Cluster_12", "Cluster_2", "Cluster_3", "Cluster_7")
East_clusters    <- c("Cluster_11", "Cluster_14", "Cluster_13", "Cluster_6", "Cluster_9", "Cluster_10")

# Assign region-based colors
W.cols      <- setNames(brewer.pal(n = length(West_clusters),    "Reds"),   West_clusters)
Midwest.cols<- setNames(brewer.pal(n = length(Midwest_clusters), "Blues"),  Midwest_clusters)
E.cols      <- setNames(brewer.pal(n = length(East_clusters),    "Greens"), East_clusters)

# Combine
cluster_colors <- c(W.cols, Midwest.cols, E.cols)

# Map data
map <- ne_states(country = c("United States of America", "Canada", "Mexico"), 
                 returnclass = 'sf')

BANS_range <- rast("data/spatial_files/banswa_abundance_seasonal_breeding_mean_2023.tif") %>% 
  project(crs(map))

ext.a <- ext(map)
BANS_range.crop <- crop(BANS_range, ext.a)
BANS_range.mask <- mask(BANS_range.crop, vect(map))
BANS_range.crop.xy <- as.data.frame(BANS_range.mask, xy=TRUE) %>% na.omit()

region_order <- c(
  "Cluster_1", "Cluster_5", 
  "Cluster_4", "Cluster_8", "Cluster_12", "Cluster_2", "Cluster_3", "Cluster_7",
  "Cluster_11", "Cluster_14", "Cluster_13", "Cluster_6", "Cluster_9", "Cluster_10"
)

# Set Pop as a factor with this order
pops <- pops %>%
  mutate(Pop = factor(Pop, levels = region_order))

# Now the legend will follow this order
p <- ggplot() +
  geom_sf(data = map, fill = "white", color = "black", size = 0.1) +
  geom_raster(data = BANS_range.crop.xy, aes(x=x, y=y), fill = "grey", alpha = 0.5) +
  geom_point(data = pops, aes(x = Long, y = Lat, color = Pop), 
             size = 1.5, position = position_jitter(width = 0.1, height = 0.1)) +
  scale_colour_manual(values = cluster_colors, name = "Spatial Clusters") +
  coord_sf(xlim = c(-140, -65), ylim = c(15, 50)) + 
  theme_minimal() +
  labs(title = "All clusters with ≥4 individuals",
       x = "Longitude",
       y = "Latitude")

p

pops <- pops %>%
  mutate(
    Region = case_when(
      Pop %in% West_clusters    ~ "West",
      Pop %in% Midwest_clusters ~ "Midwest",
      Pop %in% East_clusters    ~ "East"
    ),
    Pop = factor(Pop, levels = region_order)
  )

# Make a named vector of colors in the right order
cluster_colors <- c(W.cols, Midwest.cols, E.cols)

# Relabel with Region headers
labels_with_groups <- c(
  "West" = "--- West ---",
  setNames(West_clusters, West_clusters),
  "Midwest" = "--- Midwest ---",
  setNames(Midwest_clusters, Midwest_clusters),
  "East" = "--- East ---",
  setNames(East_clusters, East_clusters)
)

p <- ggplot() +
  geom_sf(data = map, fill = "white", color = "black", size = 0.1) +
  geom_raster(data = BANS_range.crop.xy, aes(x=x, y=y), fill = "grey", alpha = 0.5) +
  geom_point(data = pops, aes(x = Long, y = Lat, color = Pop), 
             size = 1.5, position = position_jitter(width = 0.1, height = 0.1)) +
  scale_colour_manual(
    values = cluster_colors,
    breaks = names(labels_with_groups),
    labels = labels_with_groups,
    name = "Spatial Clusters"
  ) +
  coord_sf(xlim = c(-170, -50), ylim = c(14.5462794760001, 83.1165225280001)) +
  theme_minimal() +
  labs(title = "All clusters with ≥4 individuals",
       x = "Longitude", y = "Latitude")

pdf("sample_plot_bycluster.pdf")
print(p)
dev.off()

ggsave("sample_plot.png")
#ggplotly(p, tooltip = "text")

############ Updated Sample Map Sep 2025
pops <- read.csv("data/BANS_sample_data.csv") %>% 
  dplyr::select(BGP.ID, Country_code, State, CityTown, Lat, Long) %>%
  rename(BGP_ID=BGP.ID)
pops$seq <- "Y"

to.seq <- read.csv("data/BANS_to be sequenced.csv") %>%
  dplyr::select(BGP_ID, Country_code, State, CityTown, Lat, Long)
to.seq$seq <- "N"

samples <- rbind(pops, to.seq)
samples.sf <- st_as_sf(samples, coords=c("Long", "Lat"), crs=4326)
library(mapview)

mapview(samples.sf)

# Define cluster groups
West_clusters    <- c("AK", "Yukon")
Midwest_clusters <- c("Northwest Territories", "British Columbia", "Alberta", "SK", "WA", "MT", "WY", "Manitoba","ND")
East_clusters    <- c("Ontario", "Quebec", "Newfoundland", "New Brunswick", "Prince Edward Island", "Nova Scotia")

# Assign region-based colors
W.cols      <- setNames(brewer.pal(n = length(West_clusters),    "Reds"),   West_clusters)
Midwest.cols<- setNames(brewer.pal(n = length(Midwest_clusters), "Blues"),  Midwest_clusters)
E.cols      <- setNames(brewer.pal(n = length(East_clusters),    "Greens"), East_clusters)

# Combine
cluster_colors <- c(W.cols, Midwest.cols, E.cols)

#getting base map data
map <- ne_states(country = c("United States of America", "Canada", "Mexico"), 
                 returnclass = 'sf')

BANS_range <- st_read("data/spatial_files/banswa_range_2023/banswa_range_2023.gpkg", quiet = TRUE) %>%
  filter(season == "breeding")

BANS_mask_sf <- st_intersection(BANS_range, map)
plot(BANS_mask_sf["season"])

region_order <- c(West_clusters, Midwest_clusters, East_clusters)

# Set Pop as a factor with this order
samples <- samples %>%
  mutate(Pop = factor(State, levels = region_order))

# Now the legend will follow this order

p.2 <- ggplot() +
  geom_sf(data = map, fill = "white", color = "black", size = 0.1) +
  geom_sf(data = BANS_mask_sf, fill = "grey", alpha = 0.5) +
  geom_point(data = samples %>% distinct(Long, Lat, State, seq), 
             aes(x = Long, y = Lat, 
                 fill = State,   # Area colors
                 shape = seq),   # Sequenced status
             color = "black", stroke = 0.4, size = 3) +
  scale_fill_manual(values = cluster_colors, name = "Area") +
  scale_shape_manual(values = c("Y" = 21, "N" = 24), name = "Sequenced?") +
  coord_sf(xlim = c(-170, -55), ylim = c(25, 80)) +
  theme_minimal() +
  labs(title = "All Sample Sites",
       x = "Longitude",
       y = "Latitude")

p.2

pdf("all_sample_plot_bycluster.pdf")
print(p.2)
dev.off()


