### bbsBayes - This script uses the bbsBayes package to estimate trend by BCR (bird conservation region), then creates maps/figures to compare BBS trends with trends estimated from our bird survey data, by species

### NOTES ----
# > Their spatial model looks promising
# > I CAN'T FIGURE OUT HOW TO CHECK MODEL FIT BEYOND RHAT
# Look into the first difference model, which Link and Sauer 2020 say may be better
# The models are overdispersed Poissons, and better with heavy tails (t-distribution)

### ABOUT THE FUNCTIONS ----
# > Note that BBS is favoring the first-differenced model over the slopes one now (actually going toward the spatial in the future, probably). This code uses the slopes model, but can be changed to use first-differenced.
# > Overdispersed Poisson generates a count-level random effect to account for extra-Poisson variance, following the approaches used in most of the official BBS analyses (note that NB is the default b/c uses less memory--so need to specify overdispersed Poisson explicitly).
# The generate_trends() function produces estimates of the overall percent change in the population between the first and last years of the trend period, by default. This calculation may be easier to interpret than an average annual rate of change. These percent change estimates have associated uncertainty bounds, and so can be helpful for deriving statements such as “between 2008 and 2018, the population has declined by 20 percent, but that estimate is relatively uncertain and the true decline may be as little as 2 percent or as much as 50 percent”. 
# The generate_trends() function can optionally calculate the posterior conditional probability that a population has changed more or less than some user-defined threshold(s), using the prob_decrease and prob_increase arguments. The calculate the conditional probability a species population has decreased, set the argument prob_decrease = 0 . These values can be useful for deriving statements such as “the first- difference spatial model suggests that it is extremely likely that the Scissor-tailed Flycatcher population monitored by the BBS has decreased between 1966 and 2021 (prob_decrease_0_percent > 0.999), and that there is an approximate 21% probability that the species has decreased by at least 50% in that same time period (prob_decrease_50-percent = 0.20775)”.
# generate_trends(slope = TRUE...) - The slope of this line is transformed into an average annual percent change across the time-period of interest. 
# generate_indices() --Default is "n" which for all models represents an index of the estimated annual relative abundance, scaled as the expected mean count averaged across all BBS routes and observers in a given stratum and year. 

### ABOUT THE OUTPUTS ----
# trend - Estimated median annual percent change over the trend time-period according to end point comparison of annual indices for the start_year and the end_year
# slope_trend - Estimated median annual percent change over the trend time-period, according to the slope of a linear regression through the log-transformed annual indices. (Only if slope = TRUE)
# rel_abundance - Mean annual index value across all years. An estimate of the average relative abundance of the species in the region. Can be interpreted as the predicted average count of the species in an average year on an average route by an average observer, for the years, routes, and observers in the existing data


### START CODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

### Load packages ----
pkgs <- c("bbsBayes2", "tidyverse", "sf", "patchwork", "bayesplot", "sf", "httr", "here", "magrittr")
installed_pkgs <- pkgs %in% installed.packages()
if (length(pkgs[!installed_pkgs]) > 0) install.packages(pkgs[!installed_pkgs], repos = "https://cloud.r-project.org" , dep=TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

options(mc.cores = 4, brms.backend = "cmdstanr")

### FUNCTIONS
FuncGetPark <- function(unit_code) {
  # Get park unit boundary data
  unitBoundaryURL <- paste0("https://services1.arcgis.com/fBc8EJBxQRMcHlei/ArcGIS/rest/services/IMD_Units_Generic_areas_of_analysis_(AOAs)_-_IMD_BND_ALL_UNITS_AOA_nad_py_view/FeatureServer/0/query?where=UNITCODE+%3D+%27", unit_code, "%27&objectIds=&time=&geometry=&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&distance=&units=esriSRUnit_Meter&relationParam=&outFields=*&returnGeometry=true&maxAllowableOffset=&geometryPrecision=&outSR=4326&gdbVersion=&returnDistinctValues=false&returnIdsOnly=false&returnCountOnly=false&returnExtentOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&resultOffset=&resultRecordCount=&f=geojson") # save output as WGS84
  
  if(http_status(GET(unitBoundaryURL))$category!="Success") { # if not a valid service call or if the web service is down, abort
    showModal(
      urlModal(unitBoundaryURL, title = "Data Retrieval Error", subtitle = paste0("Error retrieving Park unit boundary data from IRMA. The message from the web service is: `", http_status(GET(unitBoundaryURL))$message, "`.\n\nTo confirm this is a problem with the web service (and not the dashboard), enter the above URL in a browser and see if data successfully downloads. If an error results, email IRMA support (irma@nps.gov) for assistance."))
    )
  }
  
  imported_dat <- tryCatch(sf::st_read(unitBoundaryURL), error=function(e) print("Error retrieving Park unit boundary data")) # return error message if problems still arise with downloading from web services
  
  if(!"sf" %in% class(imported_dat)) {
    showModal(
      urlModal(unitBoundaryURL, title = "Data Retrieval Error", subtitle = paste0("Could not successfully retrieve Park unit boundary data from IRMA. To see if this is a problem with the web service (and not the dashboard), enter the above URL in a browser and see if data successfully downloads with boundary information in geojson format. If necessary, email IRMA support (irma@nps.gov) for assistance."))
    )
  }
  
  park_bound <- sf::st_transform(imported_dat, crs = "+proj=longlat +datum=WGS84") # convert to WGS84
  
  
  # Check for problems with self-intersections, etc. These problems should be fixed in the updated version of LandscapeDynamics
  if(sf::st_is_valid(park_bound) == FALSE) {
    park_bound <- sf::st_make_valid(park_bound)
  }
  return(park_bound)
}

### ONCE AND DONE ----
## Grab the most current version of bbs data
bbsBayes2::fetch_bbs_data()

## Get park boundaries selected parks 
park_bounds_list <- list()
for (i in c("BITH", "GUIS", "JELA", "PAAL", "SAAN", "VICK")) { # <<<<<<<<<<<<<<<< ENTER PARK CODES HERE
  park_bounds_list[[i]] <- FuncGetPark(i)
}
saveRDS(park_bounds_list, here::here("Data_out", "park_bounds_list.RDS")) # <<<<<<<<<<<<<<<<<< Saves the .RDS to a 'Data_out' subfolder of the working directory

# ### Pull in the GULN sf
# sf_list <- readRDS(here::here("Data_out", "sf_list.RDS"))
# guln_bound <- sf_list[["GULN"]]
# 
# ### Clip BCR's to GULN region
# 
# bcr_map_orig <- load_map("bcr") # get bbs stratum map
# bcr_map <- sf::st_transform(bcr_map_orig, crs = st_crs(guln_bound)) #convert to same coord. system as guln sf
#
# bcr_clipped <- sf::st_intersection(x = bcr_map, y = guln_bound) %>% dplyr::filter(strata_name != "BCR20") # BCR20 is only a sliver in the clipped bounds, so omit. This clips the bcrs to the GULN boundaries.
# bcr_subset <- bcr_map[guln_bound,] %>% dplyr::filter(strata_name %in% c("BCR20", "BCR21", "BCR36", "BCR37", "BCR25", "BCR26", "BCR27")) # This gives the subset of BCRs that include the parks
# saveRDS(bcr_subset, here::here("bbs_out", "bcr_subset.RDS"))
# 
# # This is the map to use for plotting results
# 
# (base_bbs_map <- ggplot() +
#   geom_sf(data = bcr_subset, aes(fill = strata_name), alpha = 0.3) +
#   geom_sf_text(data = bcr_subset, aes(label = strata_name)) +
#   # geom_sf(data = guln_bound, color = "blue", fill = NA) +
#   geom_sf(data = park_bounds_list[["BITH"]]) +
#   geom_sf_label(data = park_bounds_list[["BITH"]], label = "BITH") +
#   geom_sf(data = park_bounds_list[["JELA"]]) +
#   geom_sf_label(data = park_bounds_list[["JELA"]], label = "JELA") +
#   geom_sf(data = park_bounds_list[["GUIS"]]) +
#   geom_sf_label(data = park_bounds_list[["GUIS"]], label = "GUIS") +
#   geom_sf(data = park_bounds_list[["SAAN"]]) +
#   geom_sf_label(data = park_bounds_list[["SAAN"]], label = "SAAN") +
#   geom_sf(data = park_bounds_list[["VICK"]]) +
#   geom_sf_label(data = park_bounds_list[["VICK"]], label = "VICK") +
#   geom_sf(data = park_bounds_list[["PAAL"]]) +
#   geom_sf_label(data = park_bounds_list[["PAAL"]], label = "PAAL") +
#   scale_fill_viridis_d(guide = "none", na.value = "darkgray") +
#   theme_minimal())
# 
# saveRDS(base_bbs_map, here::here("bbs_out", "base_bbs_map.RDS"))

### START MODELS HERE ----

### FUNCTIONS ----

FuncRunBBSMod <- function(mod_dat, mod_var) {
  mp <- prepare_model(mod_dat,
                      model = "first_diff",
                      model_variant = mod_var,
                      use_pois = TRUE,
                      heavy_tailed = TRUE)
  
  m <- run_model(mp,
                 iter_warmup = 2000,
                 iter_sampling = 2000, 
                 thin = 4,
                 chains = 4)
  return(m)
}

### Species to generate maps for ----
spec_df <- data.frame(
  spec_code = c("NOCA", "CARW", "NOMO", "RBWO", "TUTI", "WEVI", "BLJA", "MODO", "AMCR", "GTGR", "YBCU", "EATO"),
  spec_name = c("Northern Cardinal", "Carolina Wren", "Northern Mockingbird", "Red-bellied Woodpecker", "Tufted Titmouse", "White-eyed Vireo", "Blue Jay", "Mourning Dove", "American Crow", "Great-tailed Grackle", "Yellow-billed Cuckoo", "Eastern Towhee")
)

first_yr <- 2010
last_yr <- 2022

# Load files
base_bbs_map <- readRDS(here::here("bbs_maps", "base_bbs_map.RDS"))
bcr_subset <- readRDS(here::here("bbs_maps", "bcr_subset.RDS"))

for(i in 1:nrow(spec_df)) {
  spec_code = spec_df$spec_code[i]
  spec_name = spec_df$spec_name[i]
# Subset to species and selected bcr's
s <- stratify(by = "bcr", 
              strata_custom = bcr_subset, 
              species = spec_name) # subset to species and selected bcr's

# which routes fall within the subset area?
s_routes <- s$routes_strata %>%
  dplyr::filter(!is.na(route)) %>%
  sf::st_as_sf(., coords = c("longitude","latitude"))
s_routes <- st_set_crs(s_routes, st_crs(bcr_subset)) # I mapped these out, and determined it's WGS84
s_routes_subset <- sf::st_intersection(x = s_routes, y = bcr_subset) # So I want to subset the data to include only these routes

# Check that only have routes that fall in clipped area. These may change by species b/c some routes will be included or excluded from analyses depending on the species
(bbs_map_routes <- base_bbs_map +
  # geom_sf(data = s_routes_subset, size = 1) +
  geom_sf(data = s_routes, shape = 1, size = 1, color = "darkgray", alpha = 0.6) +
  labs(title = "Park units and bird conservation regions (BCRs) in GULN region",
       subtitle = "Breeding bird survey routes are shown as gray circles"))
saveRDS(bbs_map_routes, here::here("bbs_out", "bbs_map_routes.RDS"))

# # Keep the route_data_id for which the route lat-long fall within the clipped area
# route_data_keep <- unique(s_routes_subset$route_data_id)
# 
# # Now subset the bbs data to include only the routes that fall within the clipped area
# s_clipped <- s
# s_clipped$meta_strata <- bcr_clipped
# s_clipped$birds_strata %<>% 
#   dplyr::filter(route_data_id %in% route_data_keep)
# s_clipped$routes_strata %<>% 
#   dplyr::filter(route_data_id %in% route_data_keep)

p <- prepare_data(s,
                  min_year = first_yr,
                  max_year = last_yr,
                  min_n_routes = 3,
                  min_max_route_years = 3,
                  min_mean_route_years = 1) # keeping the default mins

## for spatial...
sp <- prepare_spatial(p, bcr_subset, add_map = bcr_subset)
sp$spatial_data$map

m <- FuncRunBBSMod(mod_dat = sp, mod_var = "spatial") # prepare data and run model
names(m)
saveRDS(m, here::here("bbs_mods", paste0("bbs_GULN_", spec_code, "_firstdiff_spatial_fattail.RDS")))

# ..... end of spatial model fit
}


# ## for nonspatial hierarchical...
# m <- FuncRunBBSMod(mod_dat = p, mod_var = "hier")
# names(m)
# saveRDS(m, here::here("bbs_mods", paste0("bbs_GULN_", spec_code, "_slope_hier_fattail.RDS")))
# # ...... end of nonspatial hierarchical

park_bounds_list <- readRDS(here::here("Data_out", "park_bounds_list.RDS"))
bcr_subset <- readRDS(here::here("bbs_maps", "bcr_subset.RDS"))
# ### Pull in the GULN sf
sf_list <- readRDS(here::here("Data_out", "sf_list.RDS")) # <<<<<<<<<<<<<<<<<<< PICK UP FROM HERE. WHAT SCRIPT CREATES THIS?
guln_bound <- sf_list[["GULN"]]

### Check model and plot results ----
# get_summary(m) %>% View() # Check Rhats

i <- generate_indices(model_output = m, hpdi = TRUE) # Calculate annual indices of relative abundance by year for different regions.

(plot_ind <- plot_indices(indices = i, 
                  add_observed_means = TRUE, # optional argument to show raw observed mean
                  add_number_routes = TRUE))

trend_est <- generate_trends(i, quantiles = c(0.025, 0.05, 0.95, 0.975), slope = TRUE, hpdi = TRUE, prob_decrease = c(0,25)) # without 'slope = TRUE' we only get the end-point trends. Asking for probability the pop has declined, and prob. it has declined by at least 10% (e.g., from 100 to 90 or fewer)
# View(trend_est[["trends"]]) # Estimated median annual percent change over the trend time-period. The raw data column "count" gives the total number of detections of that species on that route for that year

### Create final trend map----

# Append trend results to clipped bcr sf

trend_bit <- trend_est[["trends"]] %>%
  dplyr::select(strata_name = region, start_year, end_year, mean_n_routes, annual_perc_change = slope_trend, Low95CI = slope_trend_q_0.025, High95CI = slope_trend_q_0.975, prob_decline = prob_decrease_0_percent, prob_decline_tot25perc = prob_decrease_25_percent) %>%
  dplyr::mutate(
    annual_perc_change = round(annual_perc_change, 2),
    Low95CI = round(Low95CI, 2),
    High95CI = round(High95CI, 2),
    prob_decline = round(100*prob_decline, 1),
    prob_decline_tot25perc = round(100*prob_decline_tot25perc, 1),
    hover_text = paste0("<span style='font-size:12px; font-weight:bold;'>", "Annual % change [95%CI]: ", annual_perc_change, "% [", Low95CI, ", ", High95CI, "]</span><br>Year range: ", start_year, " to ", end_year, "<br>Avg # routes/yr: ", round(mean_n_routes, 1), "<br>Probability of decline: ", prob_decline, "%<br>Probabilty of decline >25%: ", prob_decline_tot25perc, "%"))

bbs_trend <- bcr_subset %>% 
  dplyr::left_join(trend_bit, by = "strata_name")
  
# bbs_trend = sf::st_cast(bbs_trend, "MULTIPOLYGON") # this is required to convert geom_sf map to ggplotly

# Plot results
trend_map <- ggplot() +
   geom_sf(data = bbs_trend, aes(fill = annual_perc_change, text = hover_text), alpha = 0.3) +
   geom_sf_text(data = bbs_trend, aes(label = strata_name)) +
   # geom_sf(data = guln_bound, fill = NA) +
   geom_sf(data = park_bounds_list[["BITH"]]) +
   geom_sf_label(data = park_bounds_list[["BITH"]], label = "BITH") +
   geom_sf(data = park_bounds_list[["JELA"]]) +
   geom_sf_label(data = park_bounds_list[["JELA"]], label = "JELA") +
   geom_sf(data = park_bounds_list[["GUIS"]]) +
   geom_sf_label(data = park_bounds_list[["GUIS"]], label = "GUIS") +
   geom_sf(data = park_bounds_list[["SAAN"]]) +
   geom_sf_label(data = park_bounds_list[["SAAN"]], label = "SAAN") +
   geom_sf(data = park_bounds_list[["VICK"]]) +
   geom_sf_label(data = park_bounds_list[["VICK"]], label = "VICK") +
   geom_sf(data = park_bounds_list[["PAAL"]]) +
   geom_sf_label(data = park_bounds_list[["PAAL"]], label = "PAAL") +
  scale_fill_gradient2(low = "red", high = "blue", na.value = "black") +
  labs(fill = "Annual % change", title = paste0("Estimated trends for ", spec_name, " (", first_yr, " - ", last_yr, ")"), x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.width = unit(1.5, "cm"),
        plot.title=element_text(face = "bold", hjust=0.5))

trend_map


plotly::ggplotly(trend_map, tooltip = "text")
  
## AFTER I GET TREND ESTIMATES FOR PARKS, ADD THEM TO THIS MAP AND ALSO CREATE A FACET_WRAP PLOT OF TRENDS W/PARKS FIRST, AND ALSO SHOW A POINT W/ERROR BAR GRAPH SHOWING TREND AND 95% CI FOR ALL 
