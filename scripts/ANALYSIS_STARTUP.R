### ANALYSES START-UP 
## CHECK ON THE ISSUE OF ADDING ZERO'S--I THINK I ALREADY HAVE IT IN THE OBS
## ADD THE HRS-SINCE-RISE AND THE UNDERSTORY INFO. FOR ANALYSES, and 200-M FOREST
## COMBINE GUIS FOR ANALYSES BUT HAVE THE SUBUNITS AS FIXED EFFECTS
## REALLY, WE DON'T NEED TO USE THE COMMON SPECIES THRESHOLDS FOR THIS
## ADD OBSERVER AS COVARIATE??

library(tidyverse)
library(summarytools)
library(magrittr)

### Read in formatted data----
df_locs <- readRDS(here::here("Data_out", "df_locs.RDS")) # site locations
df_survey_conditions <- readRDS(here::here("Data_out", "df_survey_conditions.RDS")) # event covariate data
df_full_obs <- readRDS(here::here("Data_out", "df_full_obs.RDS"))# For each survey event, this is the total number of individuals detected per species. Includes the zero's
df_finaldat <- readRDS(here::here("Data_out", "df_finaldat.RDS")) # This is the actual count data, with time bin, distance bin, etc. for each obs. DOES NOT INCLUDE THE ZERO-COUNTS.

### Functions ----
FuncSubsetAnalysisData <- function(dat) {
  # Remove data we will not use
  if(park == "VICK") {
    
    ## For VICK, only using Daniel;s data and only data starting 2012 so can avoid messy changes of researchers and also the very odd timing of 5 surveys (the yr_visits had long timespans that overlapped each other) in 2010
    dat %<>% 
      dplyr::filter(researcher == "Twedt, Daniel" & yr >= 2012 & yr_visit != "2023_3")
  }
  
  if(park == "BITH") {
    ## REMOVE first three years of data b/c each of these first three years had a different researcher and strong researcher impact
    dat %<>% dplyr::filter(yr >= 2017)
    dat %<>% dplyr::filter(yr_visit != "2018_2")
  }
  
  if(park == "GUIS") {
    dat %<>%
      dplyr::filter(!researcher %in% c("Walker,  Jake", "Sculley,  Mike"))
  }
  
  if(park == "PAAL") { # NOTE: I only started doing this for models I ran in May 2024
    dat %<>%
      dplyr::filter(researcher == "Pruitt,  Kenneth" & yr >= 2013)
  }
  
  return(dat)
}

FuncFormatUnmarkedDist <- function(park, spec, single_visit) { # formerly called 'FuncFormatFinalDat'
  # Format the df_finaldat for use in DISTANCE SAMPLING.
  # Add a column to indicate a researcher's first year of survey
  # Combine cov. levels
  # Merge in the location data
  sub_dat <- df_finaldat
  
  # Add subunits for GUIS & SAAN-- for other parks, unit name is same as subunit
  sub_dat$subunit <- sub_dat$unit_code
  sub_dat$unit_code[sub_dat$subunit %in% c("GUIS-FL", "GUIS-MS")] <- "GUIS"
  sub_dat$subunit[sub_dat$unit_code == "SAAN"] <- "MI"
  sub_dat$subunit[sub_dat$location_name %in% c("SAAN01", "SAAN02", "SAAN03", "SAAN04", "SAAN05", "SAAN06")] <- "RA"
  
  # Grab the yr_visit info per survey event, from df_full_obs
  assign_yr_visit <- df_full_obs %>%
    dplyr::select(location_name, event_date, yr_visit) %>%
    dplyr::distinct() 
  
  sub_dat %<>% 
    dplyr::left_join(assign_yr_visit, by = c("location_name", "event_date")) %>%
    rowwise() %>%
    dplyr::mutate(visit = str_split(yr_visit, pattern = "_")[[1]][2]) %>%
    dplyr::arrange(location_name, yr_visit)
  
  # If only keeping first visit, subset the data
  if(single_visit) {
    sub_dat %<>% dplyr::filter(visit == "1")
  }
  
  # template of all surveys
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
    dplyr::select(location_name, perc_forest_100, perc_forest_200, perc_opendev_100, perc_opendev_200, hab_type_100, hab_type_200, prop_understory_cover_50, prop_understory_cover_100, prop_understory_cover_200, understory_cover_sd_50, understory_cover_sd_100, understory_cover_sd_200)
  
  # Add to the count data
  dat_Nmix <- dist_df %>% 
    dplyr::left_join(loc_dat_keep, by = c("location_name")) %>%
    dplyr::left_join(survey_dat_keep, by = c("location_name", "event_date"))
  
  # Remove data we will not use
  dat_Nmix <- FuncSubsetAnalysisData(dat_Nmix)
  # if(park == "VICK") {
  #   
  #   ## For VICK, only using Daniel;s data and only data starting 2012 so can avoid messy changes of researchers and also the very odd timing of 5 surveys (the yr_visits had long timespans that overlapped each other) in 2010
  #   dat_Nmix %<>% 
  #     dplyr::filter(researcher == "Twedt, Daniel" & yr >= 2012 & yr_visit != "2023_3")
  # }
  # 
  # if(park == "BITH") {
  #   ## REMOVE first three years of data b/c each of these first three years had a different researcher and strong researcher impact
  #   dat_Nmix %<>% filter(yr >= 2017)
  #   dat_Nmix %<>% dplyr::filter(yr_visit != "2018_2")
  # }
  # 
  # if(park == "GUIS") {
  #   dat_Nmix %<>%
  #     dplyr::filter(!researcher %in% c("Walker,  Jake", "Sculley,  Mike"))
  # }
  # 
  # if(park == "PAAL") { # NOTE: I only started doing this for models I ran in May 2024
  #   obs_dat %<>%
  #     dplyr::filter(researcher == "Pruitt,  Kenneth" & yr >= 2013)
  # }
  
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
  
  return(dat_Nmix)
  }

FuncFormatFullObs <- function(park, spec, limit_100m = TRUE) {
  # Format the df_full_obs for use
  # Add a column to indicate a researcher's first year of survey
  # Combine cov. levels
  # Merge in the location data
  
  obs_dat <- df_full_obs
  
  # Add subunits for GUIS & SAAN-- for other parks, unit name is same as subunit
  obs_dat$subunit <- obs_dat$unit_code
  obs_dat$unit_code[obs_dat$subunit %in% c("GUIS-FL", "GUIS-MS")] <- "GUIS"
  obs_dat$subunit[obs_dat$unit_code == "SAAN"] <- "MI"
  obs_dat$subunit[obs_dat$location_name %in% c("SAAN01", "SAAN02", "SAAN03", "SAAN04", "SAAN05", "SAAN06")] <- "RA"
  
  obs_dat %<>% 
    dplyr::filter(unit_code == park & species_code == spec)
  
  if(limit_100m == TRUE) { # recalc sum_indiv if limiting to 100m
    obs_dat$sum_indiv <- NULL
    
    # Drop the last distance bin and get new sums
    new_sums <- df_finaldat %>%
      dplyr::filter(distance_bin_id != "3") %>%
      droplevels(.) %>%
      dplyr::group_by(location_name, event_date, species_code) %>%
      dplyr::summarize(sum_indiv = sum(count))
    
    # Now put the new sums in obs_dat
    obs_dat %<>% 
      dplyr::left_join(new_sums, by = c("location_name", "event_date", "species_code"))
    
    obs_dat$sum_indiv[is.na(obs_dat$sum_indiv)] <- 0
    
  }
  
  # Subset the data and add covariates
  survey_dat_keep <- df_survey_conditions %>%
    dplyr::select(location_name, event_date, researcher, julian_prop, hrs_since_rise, weather_wind)
  loc_dat_keep <- df_locs %>%
    dplyr::select(location_name, perc_forest_100, perc_forest_200, perc_opendev_100, perc_opendev_200, hab_type_100, hab_type_200, prop_understory_cover_50, prop_understory_cover_100, prop_understory_cover_200, understory_cover_sd_50, understory_cover_sd_100, understory_cover_sd_200)
  
  # Add and combine covariates as needed
  obs_dat %<>%
    dplyr::left_join(loc_dat_keep, by = c("location_name")) %>%
    dplyr::left_join(survey_dat_keep, by = c("location_name", "event_date"))
  
  obs_dat <- FuncSubsetAnalysisData(obs_dat)
  # if(park == "VICK") {
  #   
  #   # For VICK, only using Daniel;s data and only data starting 2012 so can avoid messy changes of researchers and also the very odd timing of 5 surveys (the yr_visits had long timespans that overlapped each other) in 2010
  #   obs_dat %<>% 
  #     dplyr::filter(researcher == "Twedt, Daniel" & yr >= 2012 & yr_visit != "2023_3")
  # }
  # 
  # if(park == "BITH") {
  #   ## REMOVE first three years of data b/c each of these first three years had a different researcher and strong researcher impact
  #   obs_dat %<>% dplyr::filter(yr >= 2017)
  #   obs_dat %<>% dplyr::filter(yr_visit != "2018_2")
  #   
  # }
  # 
  # if(park == "GUIS") {
  #   obs_dat %<>%
  #     dplyr::filter(!researcher %in% c("Walker,  Jake", "Sculley,  Mike"))
  # }
  # 
  # if(park == "PAAL") {
  #   obs_dat %<>%
  #     dplyr::filter(researcher == "Pruitt,  Kenneth" & yr >= 2013)
  # }
  
  mod_dat <- obs_dat %>%
    dplyr::mutate(yr_c = yr - median(range(obs_dat$yr, na.rm = TRUE), na.rm = TRUE), # centered on median
                  yr_c_f = as.factor(yr_c),
                  julian_prop_c = scale(julian_prop, center = TRUE, scale = FALSE),
                  hrs_since_rise_c = scale(hrs_since_rise, center = TRUE, scale = FALSE),
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
    dplyr::mutate(survey_id = paste(unit_code, location_name, event_date)) %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::select(-weather_wind, -park_breeding)
  
  # Identify the first year for each researcher
  first_yr <- mod_dat %>%
    dplyr::select(researcher, yr) %>%
    dplyr::group_by(researcher) %>%
    dplyr::summarize(first_yr = min(yr, na.rm = TRUE))
  
  mod_dat$first_yr = 0
  
  for(i in 1:nrow(first_yr)) {
    mod_dat$first_yr[mod_dat$researcher == first_yr$researcher[i] & mod_dat$yr == first_yr$first_yr[i]] <- 1
  }
  mod_dat$first_yr <- as.factor(mod_dat$first_yr)
  
  if(park == "PAAL") {
    mod_dat %<>% dplyr::select(-hab_type_100)
  }
  
  mod_dat %<>%
    droplevels()
  
  return(mod_dat)
}

FuncFormatGDist <- function(park, spec) {
  
  
}
# Function to automatically name list elements
# actually can use tibble::lst() instead
listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}

# Below is my function for generating an AIC summary table from a list of unmarked models. The unmarked function does not work with all its data classes.
FuncFitList <- function(mod_list) {
  
  temp_aic <- lapply(mod_list, FUN = function(x) 
    c(round(x@AIC, 2), deparse(x@formula, width.cutoff = 500))) %>% 
    as.data.frame(.) %>%
    t(.) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(., var = "model") %>%
    dplyr::rename(AIC = V1, formula = V2) %>%
    dplyr::mutate(AIC = as.numeric(AIC))
  
  temp_aictab <- akaike.weights(temp_aic$AIC) %>% 
    as.data.frame(.) %>%
    dplyr::mutate(across(where(is.numeric), round, 3))
  out <- cbind(temp_aic, temp_aictab) %>%
    dplyr::arrange(deltaAIC) %>%
    dplyr::select(-rel.LL)
  
  return(out)
}

# Define new fitstats function
#   (introduced in AHM1 Section 7.9.3)
fitstats2 <- function(fm) {
  observed <- unmarked::getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  n.obs <- apply(observed,1,sum,na.rm=TRUE)
  n.pred <- apply(expected,1,sum,na.rm=TRUE)
  sse <- sum(resids^2,na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
  freeTuke.n<- sum((sqrt(n.obs)-sqrt(n.pred))^2,na.rm=TRUE)
  sse.n <- sum( (n.obs -n.pred)^2,na.rm=TRUE)
  chisq.n <- sum((n.obs - n.pred)^2 / expected,na.rm=TRUE)
  
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke,
           SSE.n = sse.n, Chisq.n = chisq.n, freemanTukey.n=freeTuke.n)
  return(out)
}

FuncPlotPredict <- function(func_name, mod, newdata, is_gdr = FALSE, incl_rem = FALSE, incl_re = FALSE) {
  abund <- predict(mod, type="lambda", newdata = newdata, appendData = FALSE) %>% 
    dplyr::rename(abund = Predicted, abund_SE = SE, abund_low = lower, abund_high = upper)
  abund <- round(abund, 3)
  
  if(single_visit == FALSE) {
    phi <- predict(mod, type="phi", newdata = newdata, appendData = FALSE) %>% 
      dplyr::rename(phi = Predicted, phi_SE = SE, phi_low = lower, phi_high = upper)
    phi <- round(phi, 3)
  }

  
  if(is_gdr) {
    rem <- unmarked::predict(mod, type="rem", newdata = newdata, appendData = FALSE) %>% 
      dplyr::rename(rem = Predicted, rem_SE = SE, rem_low = lower, rem_high = upper)
    rem <- round(rem, 3)
    
    predict_dat <- cbind(newdata, abund, phi, rem)
  } else {
    if(single_visit) {
      predict_dat <- cbind(newdata, abund)
    } else {
      predict_dat <- cbind(newdata, abund, phi)
    }
  }
  
  p <- ggplot(predict_dat, aes(x = yr, y = abund)) 
  
  if(incl_re) {
    p <- p + geom_line(aes(color = location_name))
  } else {
      p <- p + 
        geom_point() +
        geom_errorbar(aes(ymin=abund_low, ymax=abund_high), width=0) 
    }
      
  p <- p + 
    labs(title = paste0(func_name, " predictions: ", save_prefix), subtitle = toString(mod@formula[-1]), x = "Year", y = "Birds/point") + 
    expand_limits(y = c(0.00, NA)) +
    theme_bw() # Can add a facet_wrap
  return(list(predict_dat = predict_dat, p = p))
}

# FuncPlotBUP <- function(func_name, mod, dat) {
#   re <- ranef(mod)
#   
#   # BUP
#   temp_bup <- data.frame(bup(re, stat = "mean"), confint(re))
#   names(temp_bup) <- c("bup_abund", "bup_abund_CIlow", "bup_abund_CIhigh")
#   head(temp_bup)
#   
#   # PREDICT
#   temp_pred <- predict(mod, type = "lambda", nsims = 1000) %>%
#     dplyr::select(pred_abund = Predicted, pred_abund_CIlow = lower, pred_abund_CIhigh = upper) %>%
#     round(., 2)
#   head(temp_pred)
#   
#   # actual counts
#   actual <- dat %>% 
#     dplyr::group_by(location_name, yr, visit) %>% 
#     dplyr::summarize(
#       tot_obs = sum(tot_count),
#       obs_dens = round(tot_obs/dens_denom, 2)
#     )
#   head(actual)
#   
#   bup_dat <- cbind(
#     mod@data@siteCovs,
#     temp_bup, 
#     temp_pred
#   )
#   
#   # Remember that these are forced to share the same trend estimate
#   p <- ggplot(bup_dat, aes(x = yr, y = bup_abund)) +
#     labs(title = paste0(func_name, " abundances: ", save_prefix), subtitle = toString(mod@formula[-1]), x = "Year", y = "Birds/point-survey") +
#     geom_line() + 
#     geom_point(size = 0.8) + 
#     geom_ribbon(aes(x = yr, ymin = bup_abund_CIlow, ymax = bup_abund_CIhigh), alpha = 0.3) + 
#     geom_point(data = actual, aes(x = yr, y = obs_dens), size = 0.8, color = "red") +
#     theme_bw() + 
#     facet_wrap("location_name")
#   return(list(bup_dat = bup_dat, p = p))
# }
# 
# FuncPlotBUPGDR <- function(func_name, mod, dat) {
#   dens_denom <- (pi * max_bin_bound^2)/10000 # this is the number of hectares for each point plot survey area
#   re <- ranef(mod)
#   
#   # BUP
#   temp_bup <- data.frame(bup(re, stat = "mean")[,2], confint(re)[,,2])
#   names(temp_bup) <- c("bup_dens", "bup_dens_CIlow", "bup_dens_CIhigh")
#   temp_bup <- round(temp_bup/dens_denom, 2) # <<<<<<<<<< CONFIRM I CAN SIMPLY DO THIS WITH THE CI'S!
#   head(temp_bup)
#   
#   # PREDICT
#   temp_pred <- predict(mod, type = "lambda", nsims = 1000) %>%
#     dplyr::select(pred_dens = Predicted, pred_dens_CIlow = lower, pred_dens_CIhigh = upper) %>%
#     round(., 2)
#   head(temp_pred)
#   
#   # actual counts
#   actual <- dat %>% 
#     dplyr::group_by(location_name, yr, visit) %>% 
#     dplyr::summarize(
#       tot_obs = sum(tot_count),
#       obs_dens = round(tot_obs/dens_denom, 2)
#     )
#   head(actual)
#   
#   bup_dat <- cbind(
#     mod@data@siteCovs,
#     temp_bup, 
#     temp_pred
#   )
#   
#   # Remember that these are forced to share the same trend estimate
#   p <- ggplot(bup_dat, aes(x = yr, y = bup_dens)) +
#     labs(title = paste0(func_name, " densities: ", save_prefix), subtitle = toString(mod@formula[-1]), x = "Year", y = "Birds/ha") +
#     geom_line() + 
#     geom_point(size = 0.8) + 
#     geom_ribbon(aes(x = yr, ymin = bup_dens_CIlow, ymax = bup_dens_CIhigh, fill = hab_type_200), alpha = 0.3) + 
#     geom_point(data = actual, aes(x = yr, y = obs_dens), size = 0.8, color = "red") +
#     theme_bw() + 
#     facet_wrap("location_name")
#   return(list(bup_dat = bup_dat, p = p))
# }

FuncMatrix_Nmix <- function(dat, num_surv, cov_vec, response_mat = FALSE, stack = FALSE, incl_site_yr = FALSE, gdistsamp = FALSE) {
  # Function to create a list of matrices of observation covariates, for N-mixture analysis. Can also be used to generate the single matrix of the response value.
  #
  # Args:
  #   dat:  The data frame with columns corresponding to elements of cov_vec. Data should be ordered by site.
  #   num_surv: The maximum number of secondary surveys
  #   cov_vec:  Pass in a named vector of the observation covariates (for the response matrix this should be the response variable, not covariates)
  #   response_mat: TRUE if function should generate a response matrix
  #   stack: TRUE if the data should be stacked for repeated measures
  #   incl_site_yr: TRUE to include columns with site and year. Useful for checking the data, but then should be removed for analyses
  #   gdistsamp: TRUE if formatting for use with gdistsamp()
  #
  # Returns:
  #   List of matrices of covariates, or a response matrix, for N-mixture analysis
  #
  if(response_mat) {
    out <- data.frame()
    temp <- dat %>%
      dplyr::select(location_name, yr_visit, all_of(cov_vec)) %>%
      dplyr::arrange(location_name, yr_visit) %>%
      tidyr::pivot_wider(names_from = yr_visit, values_from = all_of(cov_vec), names_sort = TRUE)
    temp2 <- as.matrix(temp[,2:ncol(temp)])
    rownames(temp2) <- as.character(temp$location_name)
    
    if(stack) {
      if(gdistsamp) {
        temp3 <- FuncDistStack(temp2, surv_yrs = sort(unique(na.omit(dat$yr))), surveys_per_year=num_surv, incl_site_yr)
        temp3 %<>%
          dplyr::arrange(location_name, loc_yr)
      } else {
        temp3 <- FuncStack(temp2, nyears=length(unique(na.omit(dat$yr))), surveys_per_year=num_surv, incl_site_yr)
      }
    }
    
    out <- temp3
    
    } else {
      out <- list()
  for(i in cov_vec) {
    cat(i)
    temp <- dat %>%
      dplyr::select(location_name, yr_visit, all_of(i)) %>%
      dplyr::arrange(location_name, yr_visit) %>%
      tidyr::pivot_wider(names_from = yr_visit, values_from = all_of(i), names_sort = TRUE)
    temp2 <- as.matrix(temp[,2:ncol(temp)])
    rownames(temp2) <- as.character(temp$location_name)
    
    if(stack) {
      if(gdistsamp) {
        temp3 <- FuncDistStack(temp2, surv_yrs = sort(unique(na.omit(dat$yr))), surveys_per_year=num_surv, incl_site_yr)
        temp3 %<>%
          dplyr::arrange(location_name, loc_yr)
      } else {
        temp3 <- FuncStack(temp2, nyears=length(unique(na.omit(dat$yr))), surveys_per_year=num_surv, incl_site_yr)
      }
    }

    out[[i]] <- temp3
  }
    }
  
  
    
  return(out)
}

# FuncSubsetDat_XXX <- function(dat, spec, unit, max_bin_bound, stats_meth = "dist_time") {
#   # Use this for distance and removal models, where you require presence of species at the site to estimate detection and removal functions. Do not use for GLMMs, where you need the full zero-data incorporated (including sites where the species was never detected)
#   temp <- dat %>% 
#     dplyr::filter(unit_code %in% unit) %>% 
#     dplyr::left_join(df_locs[c("location_name", "hab_type_200")], by = "location_name") %>%
#     dplyr::left_join(df_survey_conditions[c("location_name", "event_date", "weather_wind", "weather_temperature_cs", "weather_sky_revised_num", "weather_noise_num", "weather_noise", "julian_prop")], by = c("location_name", "event_date")) %>%
#     FuncComb(., meth = stats_meth)
#   
#   # Template of sampling events at all sites where the species was detected at least once (sites where species were never detected can probably be considered structural zero's and should be included in the GLMM part of model) with their event covariate info. Fill in 0-count surveys.
#   surveys <- temp %>%
#     dplyr::select(-species_code, -distance_bin_id, -distance_bin, -time_bin_id, -count) %>%
#     distinct() 
#   
#   # These are locations at which species was detected at least once
#   bird_at_loc <- temp %>%
#     dplyr::filter(species_code == spec) %>%
#     dplyr::select(location_name) %>%
#     distinct()
#   
#   # Template should only include those locations where species was detected at least once, and should include all survey events at those locations (incl. events species not detected)
#   
#   template <- bird_at_loc %>%
#     dplyr::left_join(surveys, by = "location_name")
#   
#   dat <- temp %>%
#     dplyr::filter(species_code == spec) %>%
#     dplyr::select(-species_code) %>%
#     dplyr::right_join(template)
#   dat$count[is.na(dat$count)] <- 0
#   
#   if(!is.na(max_bin_bound) & max_bin_bound == 100) {
#     dat %<>% 
#       dplyr::filter(distance_bin_id != "3") %>%
#       droplevels(.)
#   }
#   return(dat)
# }

# This is a function that catches warnings and attaches them as attribute                                                                     
CatchWarn <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
    })
  attr(value, "val_warn") <- warn
  attr(value, "val_err") <- err
  return(value)
}

# Convert data to stacked format for N-mixture analysis of multi-year data
FuncStack <- function(input_df, nyears, surveys_per_year, incl_site_yr){
  inds <- split(1:(nyears*surveys_per_year), rep(1:nyears, each=surveys_per_year)) # for each survey year, indicates the corresponding columns of data to pull
  split_df <- lapply(1:nyears, function(i){
    out <- input_df[,inds[[i]]] %>%
      as.data.frame() 
    if(incl_site_yr) {
      out %<>%
      dplyr::mutate(
        location_name = rownames(input_df),
        yr_c = i)
    }
    names(out)[1:surveys_per_year] <- paste0("visit",1:surveys_per_year)
    out
  })
  stack_df <- do.call("rbind", split_df)
  
  if(incl_site_yr) {
    stack_df$location_name <- as.factor(stack_df$location_name)
    stack_df$yr_c <- as.factor(stack_df$yr_c)
  }
  stack_df
}

# Convert data to stacked format for N-mixture analysis of distance data (gdistance)
FuncDistStack <- function(input_df, surv_yrs, surveys_per_year, incl_site_yr){
  # Grab each distance class one-by-one
  split_df <- substack_df <- list()
  
  for(colstring in c("distbin_0", "distbin_1", "distbin_2")) { # separately for each distance bin
    grab_cols <- sub_df <- inds <- out <- stack_df <- NULL
    
    grab_cols <- colnames(input_df)[str_detect(colnames(input_df), colstring)]
    sub_df <- input_df[,grab_cols] %>% as.data.frame

    
  inds <- split(1:(length(surv_yrs)*surveys_per_year), rep(1:length(surv_yrs), each=surveys_per_year)) # for each survey year, indicates the corresponding columns of data to pull
  split_df <- lapply(1:length(surv_yrs), function(i){
    out <- sub_df[,inds[[i]]] %>% as.data.frame()
    if(incl_site_yr) {
      out %<>%
        dplyr::mutate(
          location_name = rownames(sub_df),
          loc_yr = paste(location_name, surv_yrs[i], sep = "_"))
    }
    names(out)[1:surveys_per_year] <- paste0("visit",1:surveys_per_year)
    out
  })
  stack_df <- do.call("rbind", split_df)
  
  stack_df
  if(surveys_per_year==1) {
    names(stack_df)[1] <- paste0(colstring, "_", "visit1")
  }
  if(surveys_per_year==2) {
    names(stack_df)[1:2] <- paste0(colstring, "_", c("visit1", "visit2"))
  }
  
  substack_df[[colstring]] <- stack_df
  }
  
  out_df <- dplyr::full_join(substack_df[[1]], substack_df[[2]]) %>%
    dplyr::full_join(substack_df[[3]]) 
  
  if(surveys_per_year==1) {
    out_df %<>%
    dplyr::select(distbin_0_visit1, distbin_1_visit1, distbin_2_visit1, everything())
  }
  
  if(surveys_per_year==2) {
    out_df %<>%
      dplyr::select(distbin_0_visit1, distbin_1_visit1, distbin_2_visit1, distbin_0_visit2, distbin_1_visit2, distbin_2_visit2, everything())
  }
  return(out_df)
}

