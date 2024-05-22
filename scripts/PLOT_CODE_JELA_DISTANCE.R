### This script generates bar plots of # and % of bird detections in each distance bin, by species


rm(list=ls())

library(tidyverse)
# library(ggeffects)
# library(ubms)
# library(patchwork)
library(here)

source(here::here("scripts", "ANALYSIS_STARTUP.R"))

# Generate plots for these bird species
birds_predict <- read_csv(here::here("Data_in", "jela_predictions.csv"))

df_finaldat <- readRDS(here::here("Data_out", "df_finaldat.RDS")) # This is the actual count data, with time bin, distance bin, etc. for each obs. DOES NOT INCLUDE THE ZERO-COUNTS.
df_full_obs <- readRDS(here::here("Data_out", "df_full_obs.RDS"))# For each survey event, this is the total number of individuals detected per species. Includes the zero's

# Get the JELA data for the species of interest
temp_dat <- df_finaldat %>% 
  dplyr::filter(unit_code == "JELA" & species_code %in% birds_predict$spec) %>%
  dplyr::filter(!species_code %in% c("ACFL", "GCFL", "REVI", "YBCU")) %>% # actually, too few data for these. Total counts are: ACFL = 13, GCFL=65, REVI=96, YBCU=70
  dplyr::select(location_name, species_code, count, distance_bin_id, event_date, yr)

# Grab the yr_visit info per survey event, from df_full_obs. Note that 2018 is when surveys began including ALL LOCS.
assign_yr_visit <- df_full_obs %>%
  dplyr::filter(unit_code == "JELA") %>%
  dplyr::select(location_name, event_date, yr_visit) %>%
  dplyr::distinct() 

plot_dat <- temp_dat %>% 
  dplyr::left_join(assign_yr_visit, by = c("location_name", "event_date")) %>%
  dplyr::arrange(species_code, location_name, yr_visit)

# Initial plot
ggplot(plot_dat, aes(x = yr_visit, y = count, fill = distance_bin_id)) +
  # geom_col(position = position_fill(reverse = TRUE)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  ylab("Number of detections in each distance bin") +
  xlab("Year-visit") +
  labs(title = "At JELA, number of bird detections in each distance bin (each facet is a species, each row is a survey round)",
       subtitle = "Bin distance intervals are: 0 = 0-25m, 1 = 25-50m, 2 = 50-100m, 3 = 100m+; starting 2018 (red line), a survey round included all sites") +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  geom_vline(xintercept = 9.5, color = "red") +
  scale_fill_manual(values = c("0" = "lightblue", "1" = "lightblue4", "2" = "steelblue3", "3" = "black")) +
  theme_bw() +
  facet_wrap(facets = vars(species_code), scales = "free_y")

  
  # template of all surveys <<<<<<<<<<<<<<<<< PICK UP FROM HERE. REPEAT FOR EACH SPECIES.
  all_surveys <- sub_dat %>%
    dplyr::filter(unit_code %in% park) %>%
    dplyr::select(unit_code, subunit, location_name, event_date, yr, yr_visit) %>%
    distinct()
  
  # Now actually filter for species
  sub_dat %<>%
    dplyr::filter(unit_code == park & species_code == spec) %>%
    dplyr::select(unit_code, subunit, location_name, event_date, distance_bin_id, count)
  
  # Modify and sum the distance bin counts
  summarytools::freq(sub_dat$distance_bin_id) # << EXAMINE TO SEE HOW MUCH DATA WOULD BE LOST IF TRUNCATE THE LAST DISTANCE CLASS
  
  # Drop the last distance bin
  sub_dat %<>% 
    dplyr::filter(distance_bin_id != "3") %>%
    droplevels(.)
  
  # Arrange the y-counts and fill in zero's
  dist_df <- sub_dat %>%
    dplyr::group_by(unit_code, subunit, location_name, event_date, distance_bin_id) %>%
    dplyr::summarize(sum_indiv = sum(count)) %>%
    dplyr::arrange(as.integer(distance_bin_id)) %>%
    tidyr::pivot_wider(names_from = distance_bin_id, values_from = sum_indiv, names_prefix = "distbin_") %>%
    dplyr::right_join(all_surveys, by = c("unit_code", "subunit", "location_name", "event_date"))
  
  dist_cols <- names(dist_df)[stringr::str_detect(names(dist_df), "distbin_")]
  dist_df[, dist_cols][is.na(dist_df[, dist_cols])] <- 0 # add in the zero-counts
  
  dist_df$yr <- lubridate::year(dist_df$event_date)
  
  summary(dist_df)
  
  # Add covariates data
  survey_dat_keep <- df_survey_conditions %>%
    dplyr::select(location_name, event_date, researcher, julian_prop, hrs_since_rise, weather_wind)
  loc_dat_keep <- df_locs %>%
    dplyr::select(location_name, hab_type_100, hab_type_200, prop_understory_cover_50, prop_understory_cover_100, prop_understory_cover_200, understory_cover_sd_50, understory_cover_sd_100, understory_cover_sd_200)
  
  # Add to the count data
  dat_Nmix <- dist_df %>% 
    dplyr::left_join(loc_dat_keep, by = c("location_name")) %>%
    dplyr::left_join(survey_dat_keep, by = c("location_name", "event_date"))
  
  # Remove data we will not use
  if(park == "VICK") {
    
    ## For VICK, only using Daniel;s data and only data starting 2012 so can avoid messy changes of researchers and also the very odd timing of 5 surveys (the yr_visits had long timespans that overlapped each other) in 2010
    dat_Nmix %<>% 
      dplyr::filter(researcher == "Twedt, Daniel" & yr >= 2012)
  }
  
  if(park == "BITH") {
    ## REMOVE first three years of data b/c each of these first three years had a different researcher and strong researcher impact
    dat_Nmix %<>% filter(yr >= 2017)
  }
  
  if(park == "GUIS") {
    dat_Nmix %<>%
      dplyr::filter(!researcher %in% c("Walker,  Jake", "Sculley,  Mike"))
  }
  
  # Combine covariates and calculate scaled covariates
  med_yr <- median(range(dat_Nmix$yr, na.rm = TRUE), na.rm = TRUE)
  dat_Nmix$yr_sc <- scale(dat_Nmix$yr, center = TRUE, scale = TRUE)# THIS IS MEAN-CENTERED, NOT MEDIAN-CENTERED!
  dat_Nmix$hrs_since_rise_sc <- scale(dat_Nmix$hrs_since_rise, center = TRUE, scale = TRUE)
  dat_Nmix %<>%
    dplyr::mutate(yr_c = yr - med_yr, # centered on median
                  yr_f = as.factor(yr),
                  weather_wind_comb = case_when(
                    weather_wind %in% c("0_calm", "1_smoke_drifts", "2_light_breeze") ~ "calm",
                    weather_wind %in% c( "3_constant_breeze", "4_branches_move", "5_trees_sway") ~ "windy"
                  ),
                  hab_type_200_comb = case_when(
                    hab_type_200 == "forest" ~ "forest",
                    hab_type_200 != "forest" ~ "not_forest"
                  ),
                  hab_type_100_comb = case_when(
                    hab_type_100 == "forest" ~ "forest",
                    hab_type_100 != "forest" ~ "not_forest"
                  )) %>%
    dplyr::mutate(survey_id = paste(unit_code, location_name, event_date)) 
  
  dat_Nmix$julian_prop_c = scale(dat_Nmix$julian_prop, center = TRUE, scale = FALSE)
  
  dat_Nmix %<>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::select(-weather_wind)
  
  # Identify the first year for each researcher
  first_yr <- dat_Nmix %>%
    dplyr::select(researcher, yr) %>%
    dplyr::group_by(researcher) %>%
    dplyr::summarize(first_yr = min(yr, na.rm = TRUE))
  
  dat_Nmix$first_yr = 0
  
  for(i in 1:nrow(first_yr)) {
    dat_Nmix$first_yr[dat_Nmix$researcher == first_yr$researcher[i] & dat_Nmix$yr == first_yr$first_yr[i]] <- 1
  }
  dat_Nmix$first_yr <- as.factor(dat_Nmix$first_yr)
  
  # if(park %in% c("BITH", "PAAL")) {
  #   dat_Nmix %<>% dplyr::select(-hab_type_100, -hab_type_100_comb, -hab_type_200, -hab_type_200_comb) 
  # }
  
  
  dat_Nmix %<>% 
    ungroup(.) %>% 
    droplevels() %>%
    dplyr::arrange(location_name, yr)

