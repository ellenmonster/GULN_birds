### PLOT MAPS
library(tidyverse)
library(bbsBayes2)
library(leaflet)


### Read in map parts
park_bounds_list <- readRDS(here::here("Data_out", "park_bounds_list.RDS"))
bcr_subset <- readRDS(here::here("bbs_maps", "bcr_subset.RDS"))
# ### Pull in the GULN sf
sf_list <- readRDS(here::here("Data_out", "sf_list.RDS"))
guln_bound <- sf_list[["GULN"]]

## Read in trend data
all_trends <- readRDS("~/Desktop/NPS/NPS_GULN_BIRDS/RESULTS/COMBINED_TREND_RESULTS.RDS") %>% dplyr::rename(strata_name = subunit)
all_trends %<>% dplyr::filter(source == "BBS BCR" | estim_mean_per_surv > 0.2)

### USER ENTER
map_spec <- "NOCA"
full_name <- "Northern cardinal"


bbs_trend <- bcr_subset %>% 
  dplyr::left_join(all_trends %>% dplyr::filter(spec == map_spec), by = "strata_name")

### Create point map for GULN parks
df_locs <- readRDS("~/Desktop/NPS/NPS_GULN_BIRDS/Data_out/df_locs.RDS")
temp <- df_locs
temp$unit_code[temp$unit_code == "SAAN"] <- "SAAN-MI"
temp$unit_code[temp$location_name %in% c("SAAN01", "SAAN02", "SAAN03", "SAAN04", "SAAN05", "SAAN06")] <- "SAAN-RA"
temp %<>% dplyr::group_by(unit_code) %>%
  dplyr::rename(strata_name = unit_code) %>%
  dplyr::summarize(lat = mean(latitude),
                   long = mean(longitude))

guln_trend <- temp %>%
  dplyr::left_join(all_trends %>% dplyr::filter(spec == map_spec), by = "strata_name")

bump_df <-data.frame(rbind(
  c("PAAL", -.5, 1),
  c("SAAN-MI", 1, -1.2),
  c("SAAN-RA", -.25, 1),
  c("BITH", 0.5, -1.2),
  c("VICK", 0.5,-1.2), 
  c("JELA", .5, -1),
  c("GUIS-FL", -.25, 1),
  c("GUIS-MS", .5, -1.2))
)
names(bump_df) <- c("strata_name", "hjust", "vjust")
bump_df$hjust <- as.numeric(bump_df$hjust)
bump_df$vjust <- as.numeric(bump_df$vjust)

guln_trend %<>%
  dplyr::left_join(bump_df, by = "strata_name")
  
guln_parks_sf <- sf::st_as_sf(guln_trend, coords = c("long", "lat"), crs = 4326)

p <-ggplot() +
   geom_sf(data = bbs_trend, aes(fill = annual_perc_change), alpha = 0.7) +
   geom_sf_text(data = bbs_trend, aes(label = strata_name), color = "gray") +
  geom_sf(data = guln_parks_sf, aes(fill = annual_perc_change), size = 10, alpha = 0.7, shape = 21, color = "black") +
  geom_sf_text(data = guln_parks_sf, aes(label = strata_name, hjust = hjust, vjust = vjust)) +
  scale_fill_gradient2(low = "#b2182b", high = "#2166ac", na.value = "black") +
  labs(fill = "Annual % change", title = paste0("Estimated trends for ", full_name, " (2010 to 2022)"), x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.width = unit(1.5, "cm"),
        plot.title=element_text(face = "bold", hjust=0.5))

p
saveRDS(p, here::here("RESULTS", "trend_maps", paste0("map_", map_spec, ".RDS")))

