### N-MIXTURE GDISTREMOVAL in unmarked

## >>>>>>>> REMEMBER THAT FOR FINAL ANALYSES NEED TO USE PARAMETRIC BOOTSTRAP ON GDISTSAMP AND GMULTMIX MODEL OUTPUTS B/C COULD NOT INCLUDE RANDOM EFFECTS!

### Unmarked gdistsamp()
# https://groups.google.com/g/unmarked/c/hcGgk6VxHaA
# https://darinjmcneil.weebly.com/hierarchical-distance-models.html
# https://groups.google.com/g/unmarked/c/OesoTLrawgE

### Follow this one for post-model diagnostics and predictions
# https://cran.r-project.org/web/packages/unmarked/vignettes/distsamp.html

### CHECK ON THIS:
# Do I need to run nonparametric bootstrap to get more appropriate CI's?

### WHAT TO DO:

# gdistremoval() to incorporate repeat visits, distance sampling and removal sampling
# > separate01_100_5min (both surveys)
# > separate01_100 (10min, both surveys)

rm(list = ls())

### Install development version of 'unmarked' from GitHub
remotes::install_github("rbchan/unmarked", ref="845a8e7")

### Load packages ----
library(tidyverse)
library(magrittr)
library(unmarked)
library(qpcR)
library(plotly)

source("TEST_STARTUP.R")

### USER ENTRY ----
df <- read_csv("common_spec.csv")

i = 10
(spec <- df$species_code[i])
(unit <- df$unit_code[i])
# spec = "TUTI"
# unit = "JELA"
max_bin_bound <- 100 # MUST enter some value--cannot use Inf
combine_bins <- FALSE # TRUE to combine first and second bin intervals
max_time <- 10
pair_time_int <- TRUE # TRUE to pair time intervals in removal sampling, i.e., five 2-min intervals
single_visit <- FALSE # TRUE if only using the first survey for each site-year (no TE)
  
save_prefix <- paste0(spec, "_", unit, "_", ifelse(combine_bins, "combine01", "separate01"), "_", ifelse(is.na(max_bin_bound), "unlim", paste0(max_bin_bound, "m")), "_", max_time, "min", ifelse(pair_time_int, "PAIRED", ""), ifelse(single_visit, "_1visit", "_allvisits"))

## Subset the observation raw data ----
df_finaldat <- readRDS("~/NPS PROJECTS_CURRENT/GULN_birds/R/Data_out/df_finaldat.RDS")

df_full_obs <- readRDS("~/NPS PROJECTS_CURRENT/GULN_birds/R/Data_out/df_full_obs.RDS")

# Add yr_visit and separate visit cols to df_finaldat, then restrict to max 2 visits

# Grab the yr_visit info per survey event, from df_full_obs
assign_yr_visit <- df_full_obs %>%
  dplyr::select(location_name, event_date, yr_visit) %>%
  dplyr::distinct() %>%
  dplyr::mutate(visit = as.numeric(sub(".*_", "", yr_visit)))

df_finaldat %<>% 
  dplyr::left_join(assign_yr_visit, by = c("location_name", "event_date")) %>%
  dplyr::filter(visit <= ifelse(single_visit, 1, 2)) # for VICK, there were 5 survey visits in 2010. I'm only using the first two to be consistent with the remaining years of survey

# These are all the survey events for the selected unit. Add the habitat and survey covariates.
all_surveys_cov <- df_finaldat %>%
  dplyr::filter(unit_code %in% unit) %>%
  dplyr::select(unit_code, location_name, event_date, yr, visit, yr_visit) %>% 
  distinct() %>%
  dplyr::left_join(df_locs[c("location_name", "hab_type")], by = "location_name") %>%
  dplyr::left_join(df_survey_conditions[c("location_name", "event_date", "weather_wind", "weather_temperature_cs", "weather_sky_revised_num", "weather_noise_num", "weather_noise", "julian_prop", "time")], by = c("location_name", "event_date")) %>%
  dplyr::mutate(yr = lubridate::year(event_date)) %>%
  FuncComb(., meth = "dist_time") %>%
  dplyr::select(survey_id, unit_code, location_name, event_date, start_frac_time, yr, visit, yr_visit, yr_base, julian_prop, hab_type, hab_type_comb, weather_wind_comb) %>%
  dplyr::filter(complete.cases(.))

# Get all the distance, time, and count data for a specified park AND SPECIES.
sub_dat <- df_finaldat %>%
  dplyr::filter(unit_code %in% unit & species_code == spec) %>%
  dplyr::select(unit_code, location_name, event_date, distance_bin_id, time_bin_id, count) %>%
  dplyr::mutate(
    distance_bin_id = as.integer(distance_bin_id),
    time_bin_id = as.integer(time_bin_id)+1) %>% # relabel time bins so they reflect the end of the time interval, e.g., 1 = min 0-1
  dplyr::filter(time_bin_id <= max_time)

# Truncate to 100m distance if needed
if(max_bin_bound == 100) {
  sub_dat %<>% 
    dplyr::filter(distance_bin_id != "3") %>%
    droplevels(.)
}

# Combine the first two distance bins if needed
if(combine_bins == TRUE) {
  sub_dat$distance_bin_id[sub_dat$distance_bin_id == 0] <- 1
  dist_breaks_km <- unique(c(0, 0.05, 0.1, max_bin_bound/1000))
} else {# if combining distance bins 0 and 1, then distance bin 0 is reclassified as 1
  dist_breaks_km <- unique(c(0, 0.025, 0.05, 0.1, max_bin_bound/1000))
}

# Pair the time intervals if needed
if(pair_time_int) {
  sub_dat$time_bin_id <- sub_dat$time_bin_id + sub_dat$time_bin_id%%2
}

# Add in all the surveys that were conducted, where species may not have been detected
dat_gmix <- sub_dat %>%
  dplyr::full_join(all_surveys_cov)
  
# These are all the actual count records, with time bin and distance bin information. Excludes survey events with incomplete survey covariate information. Includes true survey events with no detections (NA for count, distance_bin_id, time_bin_id)--survey occurred, but zero count
dat_gmix %<>%
  dplyr::filter(complete.cases(dat_gmix[setdiff(names(dat_gmix), c("distance_bin_id", "time_bin_id", "count"))]))
 dat_gmix$count[is.na(dat_gmix$count)] <- 0
 
# Now we need to add NA's for missing second visits within a year--should be NA's, not zero, for these, b/c the visits didn't actually happen
template <- all_surveys_cov %>% dplyr::select(location_name, yr) %>% distinct()
 
if(single_visit) {
  template_full <- template %>%
    dplyr::mutate(visit = 1)
} else {
  template_full <- merge(template, 1:2) 
}

template_full %<>% 
  dplyr::mutate(loc_yr = paste0(location_name, "_", yr))
names(template_full) <- c("location_name", "yr", "visit", "loc_yr")
 
### Starting here, modifying the code specifically for gdistsamp(), gmultimix(), and gdistremoval()----

### Response variable for distance (yDistance)
## Columns for y_st should be: "dist1visit1", "dist2visit1", "dist3visit1", dist1visit2", "dist2visit2", dist3visit2"...Each year gets its own row, but replicate visits within a year (considered the number of primary sessions, which is 2) are in a single row.

dat_gdist <- dat_gmix %>%
   dplyr::group_by(location_name, yr, visit, distance_bin_id) %>%
   dplyr::summarise(tot_count = sum(count)) %>%
  dplyr::ungroup()
 
# This includes the real 0-counts
temp_gdist_y <- dat_gdist %>%
   tidyr::pivot_wider(names_from = distance_bin_id, values_from = tot_count, names_sort = TRUE)
temp_gdist_y[is.na(temp_gdist_y)] <- 0
temp_gdist_y$`NA` <- NULL

# Now add in the NA-visits
temp_gdist_y %<>% 
  dplyr::full_join(template_full, by = c("location_name", "yr", "visit")) %>%
  dplyr::select(loc_yr, location_name, yr, visit, everything())

# Separate out visit1 from visit2, then cbind the dataframes
temp_gdist_y1 <- temp_gdist_y %>%
  dplyr::filter(visit == 1) %>%
  dplyr::select(-location_name, -yr, -visit)
colnames(temp_gdist_y1)[2:ncol(temp_gdist_y1)] <- paste0("visit1_dist", colnames(temp_gdist_y1)[2:ncol(temp_gdist_y1)])
dim(temp_gdist_y1)
head(temp_gdist_y1)

if(!single_visit) {
  temp_gdist_y2 <- temp_gdist_y %>%
    dplyr::filter(visit == 2)%>%
    dplyr::select(-location_name, -yr, -visit)
  colnames(temp_gdist_y2)[2:ncol(temp_gdist_y2)] <- paste0("visit2_dist", colnames(temp_gdist_y2)[2:ncol(temp_gdist_y2)])
  dim(temp_gdist_y2)
  head(temp_gdist_y2)
  
  yDistance <- full_join(temp_gdist_y1, temp_gdist_y2, by = c("loc_yr"))
} else {
  yDistance <- temp_gdist_y1
}

yDistance %<>%
  dplyr::arrange(loc_yr) %>%
  column_to_rownames(., var = "loc_yr") %>%
  as.matrix()
dim(yDistance)
head(yDistance)

### Response variable for removal sampling (yRemoval)
## Columns for y_st should be: "visit1t1", "visit1t2", "visit1t3"... "visit2t1", "visit2t2", "visit2t3"

dat_gremoval <- dat_gmix %>%
  dplyr::group_by(location_name, yr, visit, time_bin_id) %>%
  dplyr::summarise(tot_count = sum(count)) %>%
  dplyr::ungroup()

# This includes the real 0-counts
temp_gremoval_y <- dat_gremoval %>%
  tidyr::pivot_wider(names_from = time_bin_id, values_from = tot_count, names_sort = TRUE)
temp_gremoval_y[is.na(temp_gremoval_y)] <- 0
temp_gremoval_y$`NA` <- NULL

# Now add in the NA-visits
temp_gremoval_y %<>% 
  dplyr::full_join(template_full, by = c("location_name", "yr", "visit")) %>%
  dplyr::select(loc_yr, location_name, yr, visit, everything())

# Separate out visit1 from visit2, then cbind the dataframes
temp_gremoval_y1 <- temp_gremoval_y %>%
  dplyr::filter(visit == 1) %>%
  dplyr::select(-location_name, -yr, -visit)
colnames(temp_gremoval_y1)[2:ncol(temp_gremoval_y1)] <- paste0("visit1_t", colnames(temp_gremoval_y1)[2:ncol(temp_gremoval_y1)])
dim(temp_gremoval_y1)
head(temp_gremoval_y1)

if(!single_visit) {
  temp_gremoval_y2 <- temp_gremoval_y %>%
    dplyr::filter(visit == 2)%>%
    dplyr::select(-location_name, -yr, -visit)
  colnames(temp_gremoval_y2)[2:ncol(temp_gremoval_y2)] <- paste0("visit2_dist", colnames(temp_gremoval_y2)[2:ncol(temp_gremoval_y2)])
  dim(temp_gremoval_y2)
  head(temp_gremoval_y2)
  
  yRemoval <- full_join(temp_gremoval_y1, temp_gremoval_y2, by = c("loc_yr"))
} else {
  yRemoval <- temp_gremoval_y1
}

yRemoval %<>%
  dplyr::arrange(loc_yr) %>%
  column_to_rownames(., var = "loc_yr") %>%
  as.matrix()
dim(yRemoval)
head(yRemoval)

### `Site` covariates (includes habitat and year trend)
siteCovs <- dat_gmix %>%
  dplyr::mutate(loc_yr = paste0(location_name, "_", yr)) %>%
  dplyr::select(loc_yr, location_name, yr, hab_type, hab_type_comb) %>%
  distinct() %>%
  dplyr::right_join(data.frame(loc_yr = rownames(yDistance))) %>%
  dplyr::arrange(loc_yr)
siteCovs$yr_sc <- scale(siteCovs$yr) %>% as.numeric() # for year trend
siteCovs$yr_f <- as.factor(siteCovs$yr)
dim(siteCovs)
head(siteCovs)

### `Observation` covariates (julian date, start time, and wind)--should have MXT rows, where T is the number of site visits
yearlySiteCovs <- dat_gmix %>%
  dplyr::mutate(loc_yr = paste0(location_name, "_", yr)) %>%
  dplyr::select(loc_yr, visit, julian_prop, start_frac_time, weather_wind_comb) %>%
  distinct() %>%
  dplyr::right_join(template_full[c("loc_yr", "visit")]) %>% 
  dplyr::arrange(loc_yr, visit) # site-yr major, visit minor order

dim(yearlySiteCovs)
head(yearlySiteCovs)

# gdistremoval can't handle NA's in covs. These NA's are in here because no survey was conducted (e.g., no second survey). Set these values to mean julian_prop, mean start_frac_time, and wind = "calm". Should not affect results because there are no counts associated with these events since the survey did not "exist" 
yearlySiteCovs_gdistrem <- yearlySiteCovs
yearlySiteCovs_gdistrem$julian_prop[is.na(yearlySiteCovs_gdistrem$julian_prop)] <- mean(yearlySiteCovs_gdistrem$julian_prop, na.rm = TRUE)
yearlySiteCovs_gdistrem$start_frac_time[is.na(yearlySiteCovs_gdistrem$start_frac_time)] <- mean(yearlySiteCovs_gdistrem$start_frac_time, na.rm = TRUE)
yearlySiteCovs_gdistrem$weather_wind_comb[is.na(yearlySiteCovs_gdistrem$weather_wind_comb)] <- "calm"

# For gdistsamp() distance sampling
umf_gds <- unmarkedFrameGDS(
  y = yDistance,
  siteCovs = siteCovs,
  numPrimary = max(dat_gdist$visit, na.rm = TRUE),
  yearlySiteCovs = yearlySiteCovs,
  dist.breaks = dist_breaks_km,
  survey = "point",
  unitsIn = "km")
summary(umf_gds)

# For gmultmix() removal sampling
umf_gmm <- unmarkedFrameGMM(
  y = yRemoval,
  siteCovs = siteCovs,
  obsCovs = NULL,
  yearlySiteCovs = yearlySiteCovs,
  numPrimary = max(dat_gdist$visit, na.rm = TRUE),
  type = "removal")
summary(umf_gmm)

# For gdistremoval()
umf_gdr <- unmarkedFrameGDR(
  yDistance = yDistance,
  yRemoval = yRemoval,
  numPrimary = max(dat_gdist$visit, na.rm = TRUE),
  siteCovs = siteCovs,
  obsCovs = NULL,
  yearlySiteCovs = yearlySiteCovs_gdistrem,
  dist.breaks = dist_breaks_km,
  unitsIn = "km")
summary(umf_gdr)

use_hab_cov <- unit %in% c("JELA", "SAAN", "VICK") # For these park units, there is more than one habitat type so this can be used as a model covariate

cat("STARTING GDISTSAMP MODEL---------------------")

### FIND BEST-FIT MODEL WITH GDISTSAMP()----
## First fit gdistsamp() key function----
gdist_null_hn <- gdistsamp(lambdaformula = ~yr_sc, phiformula = ~1, pformula = ~1, data = umf_gds, keyfun = "halfnorm", output = "density", unitsOut = "ha", mixture = "P", se = FALSE) 

gdist_null_exp <- gdistsamp(lambdaformula = ~yr_sc, phiformula = ~1, pformula = ~1, data = umf_gds, keyfun = "exp", output = "density", unitsOut = "ha", mixture = "P", se = FALSE)
  
gdist_null_hz <- gdistsamp(lambdaformula = ~yr_sc, phiformula = ~1, pformula = ~1, data = umf_gds, keyfun = "hazard", output = "density", unitsOut = "ha", mixture = "P", se = FALSE)

gdist_null_unif <- gdistsamp(lambdaformula = ~yr_sc, phiformula = ~1, pformula = ~1, data = umf_gds, keyfun = "uniform", output = "density", unitsOut = "ha", mixture = "P", se = FALSE)

## Model selection for gdistsamp() key function
gdist_null_list <- listN(gdist_null_hn, gdist_null_exp, gdist_null_hz, gdist_null_unif)
View(gdist_null_sel <- FuncFitList(gdist_null_list))
  
(gdist_null_best <- get(gdist_null_sel$model[1]))
  
## Add p covariates for gdistsamp()----
(gdist_p_wind <- update(gdist_null_best, pformula = ~weather_wind_comb))

if(use_hab_cov) {
  (gdist_p_hab <- update(gdist_null_best, pformula = ~hab_type))
  
  (gdist_p_habcomb <- update(gdist_null_best, pformula = ~hab_type_comb))
  
  (gdist_p_hab_wind <- update(gdist_null_best, pformula = ~hab_type + weather_wind_comb))
  
  (gdist_p_habcomb_wind <- update(gdist_null_best, pformula = ~hab_type_comb + weather_wind_comb))
  
  gdist_p_list <- listN(gdist_null_best, gdist_p_wind, gdist_p_hab, gdist_p_habcomb, gdist_p_hab_wind, gdist_p_habcomb_wind)
} else {
  gdist_p_list <- listN(gdist_null_best, gdist_p_wind)
}

## Model selection for gdistsamp() p covariates
View(gdist_p_sel <- FuncFitList(gdist_p_list))

(gdist_p_best_temp <- get(gdist_p_sel$model[1]))
gdist_p_best <- gdist_p_best_temp

# # >>>>>>>>>> RETRY BEST MODEL WITH any other key functions <2 dAIC from best-supported model. Compare AIC's
# gdist_p_best_temp@keyfun
# View(gdist_null_sel)
# 
# gdist_p_best_temp2 <- update(gdist_p_best_temp, keyfun = "exp")
# 
# temp_list <- listN(gdist_p_best_temp, gdist_p_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# 
# (gdist_p_best <- get(temp_sel$model[1])) # <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT


## Add phi covariates for gdistsamp()----
# If needed for convergence, add: ,control = list(maxit=1e3)
(gdist_phi_julian <- update(gdist_p_best, phiformula = ~julian_prop, starts = rep(0, length(gdist_p_best@opt$par)+1),control = list(maxit=1e3)))

(gdist_phi_julian_time <- update(gdist_p_best, phiformula =~julian_prop + start_frac_time, starts = rep(0, length(gdist_p_best@opt$par)+2),control = list(maxit=1e3)))
  
(gdist_phi_julian_time2 <- update(gdist_p_best, phiformula =~julian_prop + start_frac_time+ I(start_frac_time^2), starts = rep(0, length(gdist_p_best@opt$par)+3),control = list(maxit=1e3)))
  
(gdist_phi_time <- update(gdist_p_best, phiformula =~start_frac_time, starts = rep(0, length(gdist_p_best@opt$par)+1),control = list(maxit=1e3)))
  
(gdist_phi_time2 <- update(gdist_p_best, phiformula =~start_frac_time + I(start_frac_time^2), starts = rep(0, length(gdist_p_best@opt$par)+2),control = list(maxit=1e3)))
  
## Model selection for gdistsamp() phi covariates
gdist_phi_list <- listN(gdist_null_best, gdist_p_best_temp, gdist_p_best, gdist_phi_julian, gdist_phi_julian_time, gdist_phi_julian_time2, gdist_phi_time, gdist_phi_time2)
View(gdist_phi_sel <- FuncFitList(gdist_phi_list))

(gdist_phi_best_temp <- get(gdist_phi_sel$model[1]))
gdist_phi_best <- gdist_phi_best_temp  

# # >>>>>>>>>> RETRY BEST MODEL WITH any other key functions <2 dAIC from best-supported model. Compare AIC's
# gdist_phi_best_temp@keyfun
# View(gdist_null_sel)
# 
# gdist_phi_best_temp2 <- update(gdist_phi_best_temp, keyfun = "halfnorm")
# 
# temp_list <- listN(gdist_phi_best_temp, gdist_phi_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# (gdist_phi_best <- get(temp_sel$model[1])) # <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT
# 
# # >>>>>>>>>> RETRY BEST MODEL WITH any other p covariates <2 dAIC from best-supported model. Compare AIC's
# gdist_phi_best_temp@formlist$pformula
# View(gdist_p_sel)
# 
# gdist_phi_best_temp2 <- update(gdist_phi_best_temp, pformula = ~weather_wind_comb)
# 
# temp_list <- listN(gdist_phi_best_temp, gdist_phi_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# (gdist_phi_best <- get(temp_sel$model[1])) # <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT. IF CHANGED, CHECK WITH HALF NORMAL KEY FUNCTION AGAIN.

## Add lambda covariates for gdistsamp()----
gdist_phi_best@formlist$pformula
gdist_phi_best@formlist$phiformula

# >>>>>>>>> IF PFORMULA INCLUDES HABITAT, ALWAYS ADD MODELS WITH HABITAT AS A COVARIATE FOR ABUNDANCE BUT NOT FOR DETECTION

if(use_hab_cov) {
  (gdist_abund_yr_hab <- update(gdist_phi_best, lambdaformula = ~yr_sc + hab_type, starts = rep(0, length(gdist_phi_best@opt$par)+3),control = list(maxit=1e3)))
  
  (gdist_abund_yr_habcomb <- update(gdist_phi_best, lambdaformula = ~yr_sc + hab_type_comb, starts = rep(0, length(gdist_phi_best@opt$par)+3),control = list(maxit=1e3)))
  
  (gdist_abund_yrXhab <- update(gdist_phi_best, lambdaformula = ~yr_sc + hab_type + yr_sc:hab_type, starts = rep(0, length(gdist_phi_best@opt$par)+5),control = list(maxit=1e3)))
  
  (gdist_abund_yrXhabcomb <- update(gdist_phi_best, lambdaformula = ~yr_sc*hab_type_comb, starts = rep(0, length(gdist_phi_best@opt$par)+3),control = list(maxit=1e3)))
  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  # Add a full model
  (gdist_full <- update(gdist_null_best, lambdaformula = ~yr_sc + hab_type + yr_sc:hab_type, phiformula =~julian_prop + start_frac_time+ I(start_frac_time^2), pformula = ~hab_type + weather_wind_comb, starts = rep(0, 21),control = list(maxit=1e3)))
  
  gdist_abund_list <- listN(gdist_null_best, gdist_p_best_temp, gdist_p_best, gdist_phi_best_temp, gdist_phi_best, gdist_abund_yr_hab, gdist_abund_yr_habcomb, gdist_abund_yrXhab, gdist_abund_yrXhabcomb, gdist_full)

  } else {
    
    # Add a full model
    (gdist_full <- update(gdist_null_best, lambdaformula = ~yr_sc, phiformula =~julian_prop + start_frac_time+ I(start_frac_time^2), pformula = ~weather_wind_comb, starts = rep(0, 21),control = list(maxit=1e3)))
    
    gdist_abund_list <- listN(gdist_null_best, gdist_p_best_temp, gdist_p_best, gdist_phi_best_temp, gdist_phi_best, gdist_full)
  }

## Model selection for gdistsamp() lambda covariates
View(gdist_abund_sel <- FuncFitList(gdist_abund_list))
(gdist_abund_best_temp <- get(gdist_abund_sel$model[1]))
gdist_abund_best <- gdist_abund_best_temp

# # Here we can also use a LRT to determine if a more complex model is better supported. If p<0.05, the more complex (second) model is supported.
# LRT(gdist_abund_yr_hab, gdist_abund_yrXhab)

# # >>>>>>>>>> RETRY BEST MODEL WITH any other key functions <2 dAIC from best-supported model. Compare AIC's
# gdist_abund_best_temp@keyfun
# View(gdist_null_sel)
# 
# gdist_abund_best_temp2 <- update(gdist_abund_best_temp, keyfun = "halfnorm")
# 
# temp_list <- listN(gdist_abund_best_temp, gdist_abund_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# gdist_abund_best <- get(temp_sel$model[1])# <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT
# 
# # >>>>>>>>>> RETRY BEST MODEL WITH any other p covariates <2 dAIC from best-supported model. Compare AIC's
# gdist_abund_best_temp@formlist$pformula
# View(gdist_p_sel)
#      
# gdist_abund_best_temp2 <- update(gdist_abund_best_temp, pformula = ~1)
# # gdist_abund_best_temp3 <- update(gdist_abund_best_temp, pformula = ~hab_type_comb + weather_wind_comb)
# 
# temp_list <- listN(gdist_abund_best_temp, gdist_abund_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# (gdist_abund_best <- get(temp_sel$model[1])) # <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT. IF CHANGED, CHECK WITH HALF NORMAL KEY FUNCTION AGAIN.
# 
# # >>>>>>>>>> RETRY BEST MODEL WITH any other phi covariates <2 dAIC from best-supported model. Compare AIC's
# gdist_abund_best_temp@formlist$phiformula
# View(gdist_phi_sel)
# 
# gdist_abund_best_temp2 <- update(gdist_abund_best_temp, phiformula = ~julian_prop + start_frac_time)
# 
# temp_list <- listN(gdist_abund_best_temp, gdist_abund_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# (gdist_abund_best <- get(temp_sel$model[1])) # <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT. IF CHANGED, CHECK WITH HALF NORMAL KEY FUNCTION AGAIN.

## Final model selection list----
if(use_hab_cov) {
  gdist_final_list <- listN(gdist_null_hn, gdist_null_exp, gdist_null_hz, gdist_null_unif, gdist_p_wind, gdist_p_hab, gdist_p_habcomb, gdist_p_hab_wind, gdist_p_habcomb_wind, gdist_phi_julian, gdist_phi_julian_time, gdist_phi_julian_time2, gdist_phi_time, gdist_phi_time2, gdist_abund_yr_hab, gdist_abund_yr_habcomb, gdist_abund_yrXhab, gdist_abund_yrXhabcomb, gdist_full) #, gdist_null_best, gdist_p_best, gdist_phi_best, gdist_abund_best) 
} else {
  gdist_final_list <- listN(gdist_null_hn, gdist_null_exp, gdist_null_hz, gdist_null_unif, gdist_p_wind, gdist_phi_julian, gdist_phi_julian_time, gdist_phi_julian_time2, gdist_phi_time, gdist_phi_time2, gdist_full) #, gdist_null_best, gdist_p_best, gdist_phi_best, gdist_abund_best)
}

## All-models final model selection for gdistsamp()
# If the interim 'best' models all have duplicate formula in list, can rerun without them
View(gdist_final_sel <- FuncFitList(gdist_final_list))
(gdist_best_mod_temp <- get(gdist_final_sel$model[1]))

## Final tweaking of best model for gdistsamp()----
# Try best model with mixture "NB"
gdist_best_mod_temp_nb <- update(gdist_best_mod_temp, mixture = "NB", starts = rep(0, length(gdist_best_mod_temp@opt$par) + 1), se = FALSE, control = list(maxit = 1000)) # NOTE: gdistremoval can do "ZIP", but gdistsamp cannot

# gdist_best_mod_temp_nb <- gdistsamp(lambdaformula = ~yr_sc, phiformula = ~julian_prop, pformula = ~weather_wind_comb, data = umf_gds, keyfun = "hazard", output = "density", unitsOut = "ha", mixture = "NB", starts = rep(0, length(gdist_best_mod_temp@opt$par) + 1), se = FALSE, control = list(maxit = 1000))

View(gdist_mixture_sel <-FuncFitList(listN(gdist_best_mod_temp, gdist_best_mod_temp_nb)))
gdist_best_mod_temp2 <- get(gdist_mixture_sel$model[1])

# ## >>>>>>>>>>>> If NB is selected, then try simplifying the model
# gdist_alt1 <- update(gdist_best_mod_temp2, phiformula = ~julian_prop + start_frac_time)
# gdist_alt2 <- update(gdist_alt1, pformula = ~1)
# gdist_alt3 <- update(gdist_alt2, lambdaformula = ~yr_sc)
# gdist_alt4 <- update(gdist_alt3, phiformula = ~1)
# gdist_alt5 <- update(gdist_alt3, phiformula = ~start_frac_time)
# gdist_alt6 <- update(gdist_alt3, phiformula = ~julian_prop)
# 
# temp_list <- listN(gdist_best_mod_temp2, gdist_alt1, gdist_alt2, gdist_alt3, gdist_alt4, gdist_alt5, gdist_alt6)
# View(temp_sel <- FuncFitList(temp_list))
# 
# ## <<<<<<<<<<<<<<<<<<<<<<<<<<

# Determine if need larger K
(starting_k <- gdist_best_mod_temp2@K)

# k*2
gdist_best_mod_temp2_5k <- update(gdist_best_mod_temp2, K = starting_k*5)
data.frame(rbind(coef(gdist_best_mod_temp2), coef(gdist_best_mod_temp2_2k)))

# k*4
gdist_best_mod_temp2_10k <- update(gdist_best_mod_temp2, K = starting_k*10) #, control = list(maxit=300)) 
data.frame(rbind(coef(gdist_best_mod_temp2), coef(gdist_best_mod_temp2_5k), coef(gdist_best_mod_temp2_10k))) # ...do more K if needed...

# Add SE to final model
gdist_best_mod <- update(gdist_best_mod_temp2, starts = coef(gdist_best_mod_temp2), se = TRUE) # control = list(maxit=1e5), <<<<<<<<<<<< CHANGE THIS TO THE CORRECT K

### Examine final model----
summary(gdist_best_mod) # Remember that year is scaled, so the estimated trend is change per sd of year covariate, NOT change per year

## MOD NOTES

(GDIST_MOD_NOTES <- list(
  best_mod = paste0("Best gdistsamp model key function = ", gdist_best_mod@keyfun, "; formula = ", toString(gdist_best_mod@formula[-1])),
  gof = "GOOD! RANGED FROM 0.16 TO 0.42", 
  other = "NB BUT ESTIMATES STABILIZED!AIC DIFFERENCE FROM POISSON WAS <2AIC"))# Use this for any notes 

## Check model fit
gdist_best_mod_nm <- update(gdist_best_mod, method='Nelder-Mead')
(gdist_best_mod_pb <- parboot(gdist_best_mod_nm, fitstats2, nsim = 30, report = 1, ncores = 4))

# ## >>>>>>>>>>>>>> IF POOR MODEL FIT, REASSESS AND TRY OTHER MODELS!
# alt1 <- update(gdist_null_best, lambdaformula = ~yr_sc, phiformula =~julian_prop + start_frac_time, pformula = ~weather_wind_comb)
# alt2 <- update(gdist_null_best, lambdaformula = ~yr_sc, phiformula =~start_frac_time, pformula = ~weather_wind_comb)
# alt3 <- update(gdist_null_best, lambdaformula = ~yr_sc, phiformula =~julian_prop, pformula = ~weather_wind_comb)
# 
# View(gmult_final_sel <-FuncFitList(listN(gdist_best_mod, alt1, alt2, alt3)))


## PLOT MODEL OUTPUTS ----
## Create data frame of model predictions for specific covariate combinations
gdist_best_mod@formlist # <<<<< SEE WHAT COVARIATES TO INCLUDE

df_covs <- data.frame(
  julian_prop = mean(umf_gds@yearlySiteCovs$julian_prop, na.rm = TRUE),
  start_frac_time = mean(umf_gds@yearlySiteCovs$start_frac_time, na.rm = TRUE), 
  weather_wind_comb = "calm", 
  hab_type = unique(umf_gds@siteCovs$hab_type))
df_covs$hab_type_comb <- "not_forest"
df_covs$hab_type_comb[df_covs$hab_type == "forest"] <- "forest"

df_yr <- unique(umf_gds@siteCovs[c("yr_sc", "yr")]) %>% arrange(yr_sc)

newdata <- merge(df_covs, df_yr)
  
## PLOT MODEL PREDICTIONS for specific covariate combinations
predict_out <- FuncPlotPredict(func_name = "gdistsamp", mod = gdist_best_mod, newdata = newdata, incl_re = FALSE)
gdist_predict_dat <- predict_out$predict_dat
gdist_predict_plot <- predict_out$p
ggplotly(gdist_predict_plot)

# ## PLOT DISTANCE DETECTION FUNCTION. Get the predicted distance detection function parameters (this is NOT the probability of detection!)
# IF BEST MODEL HAS COVARIATES, THEN SHOULD PREDICT AT MEAN LEVELS
gdist_best_mod@keyfun

det_shape <- gdist_best_mod@estimates[3]@estimates %>% exp()
# If there is a single covariate on detection... MODIFY AS NECESSARY
if(length(det_shape) > 1) {
  det_shape <- gdist_best_mod@estimates[3]@estimates[1] %>% exp()
  det_shape2 <- (gdist_best_mod@estimates[3]@estimates[1] + gdist_best_mod@estimates[3]@estimates[2]) %>% exp()
  det_shape3 <- (gdist_best_mod@estimates[3]@estimates[1] + gdist_best_mod@estimates[3]@estimates[3]) %>% exp()
  det_shape4 <- (gdist_best_mod@estimates[3]@estimates[1] + gdist_best_mod@estimates[3]@estimates[4]) %>% exp()
  det_shape5 <- (gdist_best_mod@estimates[3]@estimates[1] + gdist_best_mod@estimates[3]@estimates[5]) %>% exp()
}

if(gdist_best_mod@keyfun == "hazard") {
  det_scale <- gdist_best_mod@estimates[4]@estimates %>% exp()
}

# Half normal
(gdist_detect_plot <- plot(function(x) unmarked::gxhn(x, sigma = det_shape), 0, max_bin_bound/1000, main = paste0("HALF-NORM Function: ", save_prefix), xlab="Distance (km)", ylab="Detection probability", cex.lab=0.7, cex.axis=0.7, las=1)) # col = "blue", add = TRUE, to add this line to another plot

# Half normal with COVARIATE (SAVES AS FUNCTION)
gdist_detect_plot <- function() {
  det_shape <-0.0377943
  det_shape2 <-0.02785479 
  plot(function(x) unmarked::gxhn(x, sigma = det_shape), 0, max_bin_bound/1000, main = paste0("HALF-NORM Function: ", save_prefix), xlab="Distance (km)", ylab="Detection probability", cex.lab=0.7, cex.axis=0.7, las=1)
  plot(function(x) unmarked::gxhn(x, sigma = det_shape2), 0, max_bin_bound/1000, col = "green", add = TRUE)
  legend('topright', c("Calm", "Windy"), col=c("black", "green"), lty=1, cex=1)
}
gdist_detect_plot()

# Hazard
(gdist_detect_plot <- plot(function(x) unmarked::gxhaz(x, shape = det_shape, scale = det_scale), 0, max_bin_bound/1000, main = paste0("HAZARD Function: ", save_prefix), xlab="Distance (km)", ylab="Detection probability", cex.lab=0.7, cex.axis=0.7, las=1))

# Hazard with COVARIATE (SAVES AS FUNCTION)
gdist_detect_plot <- function() {
  det_shape <- 0.0007122781
  det_shape2 <- 0.02675571 
  det_shape3 <- 0.001653192 
  det_shape4 <- 0.004203732 
  det_shape5 <- 0.005407489 
  det_scale <- 0.9352508  
  plot(function(x) unmarked::gxhaz(x, shape = det_shape, scale = det_scale), 0, max_bin_bound/1000, main = paste0("HAZARD Function: ", save_prefix), xlab="Distance (km)", ylab="Detection probability", cex.lab=0.7, cex.axis=0.7, las=1)
plot(function(x) unmarked::gxhaz(x, shape = det_shape2, scale = det_scale), 0, max_bin_bound/1000, col = "green", add = TRUE)
plot(function(x) unmarked::gxhaz(x, shape = det_shape3, scale = det_scale), 0, max_bin_bound/1000, col = "orange", add = TRUE)
plot(function(x) unmarked::gxhaz(x, shape = det_shape4, scale = det_scale), 0, max_bin_bound/1000, col = "brown", add = TRUE)
plot(function(x) unmarked::gxhaz(x, shape = det_shape5, scale = det_scale), 0, max_bin_bound/1000, col = "red", add = TRUE)

legend('topright', c("Developed", "Forest", "Mixed", "Open", "Shrub"), col=c("black", "green", "orange", "brown", "red"), lty=1, cex=1)
}
gdist_detect_plot()

# Exp
(gdist_detect_plot <- plot(function(x) unmarked::gxexp(x, rate = det_shape), 0, max_bin_bound/1000, main = paste0("NEG. EXP. Function: ", save_prefix), xlab="Distance (km)", ylab="Detection probability", cex.lab=0.7, cex.axis=0.7, las=1))

# Exp with COVARIATE (SAVES AS FUNCTION)
gdist_detect_plot <- function() {
  det_shape <-0.03170303
  det_shape2 <-0.01237553 
  plot(function(x) unmarked::gxexp(x, rate = det_shape), 0, max_bin_bound/1000, main = paste0("NEG. EXP. Function: ", save_prefix), xlab="Distance (km)", ylab="Detection probability", cex.lab=0.7, cex.axis=0.7, las=1)
  plot(function(x) unmarked::gxexp(x, rate = det_shape2), 0, max_bin_bound/1000, col = "green", add = TRUE)
  legend('topright', c("Calm", "Windy"), col=c("black", "green"), lty=1, cex=1)
}


## PLOT BUP DENSITIES per site-year 
# These are the estimated abundances per site-year. Approximately the same as multiplying the predicted density by survey area, but also incorporates information about the other covariates and observed counts?
# Ranef gives the empirical Bayes estimates of the superpopulation size. Gives similar mean to predict() but it seems to also consider the actual observed data. >>>> I think I can convert to density just by dividing by survey point area, but the documentation says it's the estimated superpopulation size.

(out <- FuncPlotBUP(func_name = "GDISTSAMP", mod = gdist_best_mod, dat = dat_gdist))
gdist_bup_plot <- out$p
gdist_bup_dat <- out$bup_dat

mean(out$bup_dat$bup_dens) # mean estimated density across all surveys

### SAVE ALL MODELS----
if(use_hab_cov) {
  save(spec, unit, max_bin_bound, combine_bins, max_time, pair_time_int, single_visit, dat_gmix, dat_gdist, umf_gds, gdist_null_hn, gdist_null_exp, gdist_null_hz, gdist_null_unif, gdist_p_wind, gdist_p_hab, gdist_p_habcomb, gdist_p_hab_wind, gdist_p_habcomb_wind, gdist_phi_julian, gdist_phi_julian_time, gdist_phi_julian_time2, gdist_phi_time, gdist_phi_time2, gdist_abund_yr_hab, gdist_abund_yr_habcomb, gdist_abund_yrXhab, gdist_abund_yrXhabcomb, gdist_full, gdist_final_sel, gdist_null_best, gdist_p_best, gdist_phi_best, gdist_abund_best, gdist_best_mod, gdist_predict_dat, gdist_predict_plot, gdist_bup_plot, gdist_bup_dat, gdist_detect_plot, GDIST_MOD_NOTES, file = here::here("Model_fits", "Gdistsamp", paste0("GDISTSAMP_", save_prefix, ".RData"))) 
  # } else {
    save(spec, unit, max_bin_bound, combine_bins, max_time, pair_time_int, single_visit, dat_gmix, dat_gdist, umf_gds, gdist_null_hn, gdist_null_exp, gdist_null_hz, gdist_null_unif, gdist_p_wind, gdist_phi_julian, gdist_phi_julian_time, gdist_phi_julian_time2, gdist_phi_time, gdist_phi_time2, gdist_full, gdist_final_sel, gdist_null_best, gdist_p_best, gdist_phi_best, gdist_abund_best, gdist_best_mod, gdist_best_mod_pb, gdist_predict_dat, gdist_predict_plot, gdist_bup_plot, gdist_bup_dat, gdist_detect_plot, GDIST_MOD_NOTES, file = here::here("Model_fits", "Gdistsamp", paste0("GDISTSAMP_", save_prefix, ".RData")))
  }



cat("STARTING GMULTMIX MODEL---------------------")

### FIND BEST-FIT MODEL WITH GMULTMIX()----
## Null model
gmult_null <- gmultmix(lambdaformula = ~yr_sc, phiformula = ~1, pformula = ~1, data = umf_gmm, se = FALSE) 

## Add p covariates for gmultmix()----
(gmult_p_wind <- gmultmix(lambdaformula = ~yr_sc, phiformula = ~1, pformula = ~weather_wind_comb, data = umf_gmm, se = FALSE))

if(use_hab_cov) {
  (gmult_p_hab <- gmultmix(lambdaformula = ~yr_sc, phiformula = ~1, pformula = ~hab_type, data = umf_gmm, se = FALSE))
  
  (gmult_p_habcomb <- gmultmix(lambdaformula = ~yr_sc, phiformula = ~1, pformula = ~hab_type_comb, data = umf_gmm, se = FALSE))
  
  (gmult_p_hab_wind <- gmultmix(lambdaformula = ~yr_sc, phiformula = ~1, pformula = ~hab_type + weather_wind_comb, data = umf_gmm, se = FALSE))
  
  (gmult_p_habcomb_wind <- gmultmix(lambdaformula = ~yr_sc, phiformula = ~1, pformula = ~hab_type_comb + weather_wind_comb, data = umf_gmm, se = FALSE))
  
  gmult_p_list <- listN(gmult_null, gmult_p_wind, gmult_p_hab, gmult_p_habcomb, gmult_p_hab_wind, gmult_p_habcomb_wind)
} else {
  gmult_p_list <- listN(gmult_null, gmult_p_wind)
}

## Model selection for gdistsamp() p covariates
View(gmult_p_sel <- FuncFitList(gmult_p_list))

(gmult_p_best <- get(gmult_p_sel$model[1]))

## Add phi covariates for gmultmix()----
# If needed for convergence, add: ,control = list(maxit=1e3)
(gmult_phi_julian <- update(gmult_p_best, phiformula = ~julian_prop))

(gmult_phi_julian_time <- update(gmult_p_best, phiformula =~julian_prop + start_frac_time)) #, starts = rep(0, 7),control = list(maxit=1e3)))

(gmult_phi_julian_time2 <- update(gmult_p_best, phiformula =~julian_prop + start_frac_time+ I(start_frac_time^2))) #, starts = rep(0, 8),control = list(maxit=1e3)))

(gmult_phi_time <- update(gmult_p_best, phiformula =~start_frac_time))

(gmult_phi_time2 <- update(gmult_p_best, phiformula =~start_frac_time + I(start_frac_time^2))) #, starts = rep(0, 7),control = list(maxit=1e3)))

## Model selection for gmultmix() phi covariates
gmult_phi_list <- listN(gmult_null, gmult_p_best, gmult_phi_julian, gmult_phi_julian_time, gmult_phi_julian_time2, gmult_phi_time, gmult_phi_time2)
View(gmult_phi_sel <- FuncFitList(gmult_phi_list))

(gmult_phi_best_temp <- get(gmult_phi_sel$model[1]))
gmult_phi_best <- gmult_phi_best_temp

# # >>>>>>>>>> RETRY BEST MODEL WITH any other p covariates <2 dAIC from best-supported model. Compare AIC's
# gmult_phi_best_temp@formlist$pformula
# View(gmult_p_sel)
# 
# gmult_phi_best_temp2 <- update(gmult_phi_best_temp, pformula = ~hab_type + weather_wind_comb)
# 
# temp_list <- listN(gmult_phi_best_temp, gmult_phi_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# 
# (gmult_phi_best <- get(temp_sel$model[1])) # <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT. IF CHANGED, CHECK WITH HALF NORMAL KEY FUNCTION AGAIN.


## Add lambda covariates for gmultmix()----
if(use_hab_cov) {
  (gmult_abund_yr_hab <- update(gmult_phi_best, lambdaformula = ~yr_sc + hab_type))
  
  (gmult_abund_yr_habcomb <- update(gmult_phi_best, lambdaformula = ~yr_sc + hab_type_comb))
  
  (gmult_abund_yrXhab <- update(gmult_phi_best, lambdaformula = ~yr_sc*hab_type))
  
  (gmult_abund_yrXhabcomb <- update(gmult_phi_best, lambdaformula = ~yr_sc*hab_type_comb))
  
  # Add a full model
  (gmult_full <- update(gmult_null, lambdaformula = ~yr_sc*hab_type, phiformula =~julian_prop + start_frac_time+ I(start_frac_time^2), pformula = ~hab_type + weather_wind_comb)) #, starts = rep(0, 20),control = list(maxit=1e3)))
  
#  # Add the model that was the best gdistsamp model (in case it's not included here)--NEED TO ADD IT IN THE FIT LISTS AND SAVE LIST THEN.
  # gdist_best_mod@formlist
  # (gmult_dist_mod <- gmultmix(lambdaformula = as.formula(gdist_best_mod@formlist$lambdaformula), phiformula = as.formula(gdist_best_mod@formlist$phiformula), pformula = as.formula(gdist_best_mod@formlist$pformula), data = umf_gmm, se = FALSE))
  
  gmult_abund_list <- listN(gmult_null, gmult_p_best, gmult_phi_best_temp, gmult_phi_best, gmult_abund_yr_hab, gmult_abund_yr_habcomb, gmult_abund_yrXhab, gmult_abund_yrXhabcomb, gmult_full)
  
} else {
  
  # Add a full model
  (gmult_full <- update(gmult_null, lambdaformula = ~yr_sc, phiformula =~julian_prop + start_frac_time+ I(start_frac_time^2), pformula = ~weather_wind_comb, starts = rep(0, 8),control = list(maxit=1e3)))
  
  # # Add the model that was the best gdistsamp model (in case it's not included here)--NEED TO ADD IT IN THE FIT LISTS AND SAVE LIST THEN.
  # gdist_best_mod@formlist
  # (gmult_dist_mod <- update(gmult_null, lambdaformula = gdist_best_mod@formlist$lambdaformula, phiformula =gdist_best_mod@formlist$phiformula, pformula = gdist_best_mod@formlist$pformula))
  
  gmult_abund_list <- listN(gmult_null, gmult_p_best, gmult_phi_best_temp, gmult_phi_best, gmult_full)
}

## Model selection for gmultmix() lambda covariates
View(gmult_abund_sel <- FuncFitList(gmult_abund_list))

# # Here we can also use a LRT to determine if a more complex model is better supported. If p<0.05, the more complex model is supported as being a statistically significant better fit
# LRT(gmult_full, gmult_null)

(gmult_abund_best_temp <- get(gmult_abund_sel$model[1]))
gmult_abund_best <- gmult_abund_best_temp

# # >>>>>>>>>> RETRY BEST MODEL WITH any other p covariates <2 dAIC from best-supported model. Compare AIC's
# gmult_abund_best_temp@formlist$pformula
# View(gmult_p_sel)
# 
# gmult_abund_best_temp2 <- update(gmult_abund_best_temp, pformula = ~1, starts = rep(0, 5),control = list(maxit=1e3))
# gmult_abund_best_temp3 <- update(gmult_abund_best_temp, pformula = ~hab_type_comb, starts = rep(0, 8),control = list(maxit=1e3))
# 
# temp_list <- listN(gmult_abund_best_temp, gmult_abund_best_temp2, gmult_abund_best_temp3)
# View(temp_sel <- FuncFitList(temp_list))
# 
# # >>>>>>>>>> RETRY BEST MODEL WITH any other phi covariates <2 dAIC from best-supported model. Compare AIC's
# gmult_abund_best_temp@formlist$phiformula
# View(gmult_phi_sel)
# 
# gmult_abund_best_temp2 <- update(gmult_abund_best_temp, phiformula = ~julian_prop + start_frac_time, starts = rep(0, 8),control = list(maxit=1e3))
# 
# temp_list <- listN(gmult_abund_best_temp, gmult_abund_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# 
# gmult_abund_best <- get(temp_sel$model[1]) # <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT. IF CHANGED, CHECK WITH HALF NORMAL KEY FUNCTION AGAIN.
# # gmult_abund_best <- gmult_abund_best_temp 
# gmult_abund_best

## Final model selection list----
if(use_hab_cov) {
  gmult_final_list <- listN(gmult_null, gmult_p_wind, gmult_p_hab, gmult_p_habcomb, gmult_p_hab_wind, gmult_p_habcomb_wind, gmult_phi_julian, gmult_phi_julian_time, gmult_phi_julian_time2, gmult_phi_time, gmult_phi_time2, gmult_abund_yr_hab, gmult_abund_yr_habcomb, gmult_abund_yrXhab, gmult_abund_yrXhabcomb, gmult_full) #, gmult_p_best, gmult_phi_best, gmult_abund_best)
} else {
  gmult_final_list <- listN(gmult_null, gmult_p_wind, gmult_phi_julian, gmult_phi_julian_time, gmult_phi_julian_time2, gmult_phi_time, gmult_phi_time2, gmult_full) #, gmult_p_best, gmult_phi_best, gmult_abund_best)
}

## All-models final model selection for gmultmix()
# If the interim 'best' models all have duplicate formula in list, can rerun without them
View(gmult_final_sel <- FuncFitList(gmult_final_list))
(gmult_best_mod_temp <- get(gmult_final_sel$model[1]))

## Final tweaking of best model for gmultmix()----
# Try best model with mixture "NB"
gmult_best_mod_temp_nb <- update(gmult_best_mod_temp, mixture = "NB", starts = rep(0, 8)) #,control = list(maxit=1e3)) # NOTE: gmultremoval can do "ZIP", but gmultmix cannot

View(gmult_mixture_sel <-FuncFitList(listN(gmult_best_mod_temp, gmult_best_mod_temp_nb)))

gmult_best_mod_temp2 <- get(gmult_mixture_sel$model[1])

# Determine if need larger K
(starting_k <- gmult_best_mod_temp2@K)

# k*2
gmult_best_mod_temp2_20k <- update(gmult_best_mod_temp2, K = starting_k*20)
data.frame(rbind(coef(gmult_best_mod_temp2), coef(gmult_best_mod_temp2_20k)))
           
# k*4
gmult_best_mod_temp2_30k <- update(gmult_best_mod_temp2, K = starting_k*30) 
data.frame(rbind(coef(gmult_best_mod_temp2), coef(gmult_best_mod_temp2_20k), coef(gmult_best_mod_temp2_30k))) # ...do more K if needed (remember to change the model name and also the K value)...

# Add SE to final model
gmult_best_mod <- update(gmult_best_mod_temp2, starts = coef(gmult_best_mod_temp2), se = TRUE) # , control = list(maxit=1e5),  <<<<<<<<<<<< CHANGE THIS TO THE CORRECT K AS NEEDED

### Examine final model----
summary(gmult_best_mod) # >>>>>>>>>>>>>>>> Remember that year is scaled, so the estimated trend is change per sd of year covariate, NOT change per year


## Check model fit
# gmult_best_mod_pb <- NULL # if it didn't work...

gmult_best_mod_nm <- update(gmult_best_mod, method='Nelder-Mead')
(gmult_best_mod_pb <- parboot(gmult_best_mod_nm, fitstats2, nsim = 30, report = 1, ncores = 4)) # switching to `Nelder-Mead` allows this to work better. This can be SLOW! For final analyses, use much larger nsim

# ## >>>>>>>>>>>>>> IF POOR MODEL FIT, REASSESS AND TRY OTHER MODELS!
# gmult_alt1 <- gmultmix(lambdaformula = ~yr_sc*hab_type, phiformula = ~julian_prop, pformula = ~hab_type + weather_wind_comb, data = umf_gmm, se = TRUE)
# gmult_alt2 <- gmultmix(lambdaformula = ~yr_sc + hab_type, phiformula = ~julian_prop, pformula = ~hab_type + weather_wind_comb, data = umf_gmm, se = TRUE)
# gmult_alt3 <- gmultmix(lambdaformula = ~yr_sc, phiformula = ~julian_prop, pformula = ~hab_type + weather_wind_comb, data = umf_gmm, se = TRUE)
# gmult_alt4 <- gmultmix(lambdaformula = ~yr_sc + hab_type, phiformula = ~julian_prop, pformula = ~weather_wind_comb, data = umf_gmm, se = TRUE)
# gmult_alt5 <- gmultmix(lambdaformula = ~yr_sc + hab_type, phiformula = ~julian_prop + start_frac_time, pformula = ~hab_type + weather_wind_comb, data = umf_gmm, se = TRUE)
# gmult_alt6 <- gmultmix(lambdaformula = ~yr_sc + hab_type, phiformula = ~julian_prop + start_frac_time + I(start_frac_time^2), pformula = ~hab_type + weather_wind_comb, data = umf_gmm, se = TRUE)
# 
# View(gmult_final_sel <-FuncFitList(listN(gmult_null, gmult_p_wind, gmult_p_hab, gmult_p_habcomb, gmult_p_hab_wind, gmult_p_habcomb_wind, gmult_phi_julian, gmult_phi_julian_time, gmult_phi_julian_time2, gmult_phi_time, gmult_phi_time2, gmult_abund_yr, gmult_abund_yr_hab, gmult_abund_yr_habcomb, gmult_abund_yrXhab, gmult_abund_yrXhabcomb, gmult_full, gmult_alt1, gmult_alt2, gmult_alt3, gmult_alt4, gmult_alt5, gmult_alt6))) 
# 
# (gmult_best_mod_temp <- get(gmult_final_sel$model[1]))
# 
# # <<<<<<<<<<<<<<< FROM HERE GO BACK AGAIN TO DO THE FINAL MODEL TWEAKS!

## NOTES
GMULT_MOD_NOTES <- list(
  best_mod = paste0("Best gmultmix model formula = ", toString(gmult_best_mod@formula[-1])), # <<<<<<<<< CHANGE GOF
  gof = "NO GOOD. RANGED 0 -0.48 (3 OF 6 WERE 0)",
  other = "BEST MODEL IS NB.")# Use this for any notes 

## PLOT MODEL OUTPUTS----
## Create data frame of model predictions for specific covariate combinations
gmult_best_mod@formlist # <<<<< SEE WHAT COVARIATES TO INCLUDE

df_covs <- data.frame(
  julian_prop = mean(umf_gmm@yearlySiteCovs$julian_prop, na.rm = TRUE),
  start_frac_time = mean(umf_gmm@yearlySiteCovs$start_frac_time, na.rm = TRUE), 
  weather_wind_comb = "calm", # this is the most common wind level
  hab_type = unique(umf_gmm@siteCovs$hab_type)
)
df_yr <- unique(umf_gmm@siteCovs[c("yr_sc", "yr")]) %>% arrange(yr_sc) # this is so we can make predictions by yr_sc, but see the corresponding yr

(newdata <- merge(df_covs, df_yr))

## Plot model predictions for specific covariate combinations
predict_out <- FuncPlotPredict(func_name = "gmultmix", mod = gmult_best_mod, newdata = newdata, incl_re = FALSE)
gmult_predict_dat <- predict_out$predict_dat
gmult_predict_plot <- predict_out$p
ggplotly(gmult_predict_plot)

## Plot bup densities per site-year 
# These are the estimated abundances per site-year. Differs from predict() in that it incorporates information about the other covariates and observed counts? Note that predict() gives output as predicted density, if model is set to output density
# Ranef gives the empirical Bayes estimates of the superpopulation size. Gives similar mean to predict() but it seems to also consider the actual observed data. >>>> I think I can convert to density just by dividing by survey point area, but the documentation says it's the estimated superpopulation size.

(out <- FuncPlotBUP(func_name = "GMULTMIX", mod = gmult_best_mod, dat = dat_gremoval))
gmult_bup_plot <- out$p
gmult_bup_dat <- out$bup_dat

### SAVE ALL MODELS----
if(use_hab_cov) {
  save(spec, unit, max_bin_bound, combine_bins, max_time, pair_time_int, single_visit, dat_gmix, dat_gremoval, umf_gmm, gmult_null, gmult_p_wind, gmult_p_hab, gmult_p_habcomb, gmult_p_hab_wind, gmult_p_habcomb_wind, gmult_phi_julian, gmult_phi_julian_time, gmult_phi_julian_time2, gmult_phi_time, gmult_phi_time2, gmult_abund_yr_hab, gmult_abund_yr_habcomb, gmult_abund_yrXhab, gmult_abund_yrXhabcomb, gmult_full, gmult_final_sel, gmult_p_best, gmult_phi_best, gmult_abund_best, gmult_best_mod, gmult_predict_dat, gmult_predict_plot, gmult_bup_plot, gmult_bup_dat, GMULT_MOD_NOTES, file = here::here("Model_fits", "gmultmix", paste0("gmultmix_", save_prefix, ".RData"))) # gmult_best_mod_pb, 
} else {
  save(spec, unit, max_bin_bound, combine_bins, max_time, pair_time_int, single_visit, dat_gmix, dat_gremoval, umf_gmm, gmult_null, gmult_p_wind, gmult_phi_julian, gmult_phi_julian_time, gmult_phi_julian_time2, gmult_phi_time, gmult_phi_time2, gmult_full, gmult_final_sel, gmult_p_best, gmult_phi_best, gmult_abund_best, gmult_best_mod, gmult_best_mod_pb, gmult_predict_dat, gmult_predict_plot, gmult_bup_plot, gmult_bup_dat, GMULT_MOD_NOTES, file = here::here("Model_fits", "gmultmix", paste0("gmultmix_", save_prefix, ".RData")))
}





cat("STARTING GDISTREMOVAL MODEL---------------------")

### FIND BEST-FIT MODEL WITH GDISTREMOVAL()----

# These are the best-fit models from distance and removal models. Make sure these are included in model selection here. Gdistremoval() SHOULD also pick a similar best model.
gdist_best_mod@keyfun # half norm
gdist_best_mod@formula # ~yr_sc ~ 1 ~ weather_wind_comb
gmult_best_mod@formula # ~yr_sc ~ 1 ~ weather_wind_comb

# REMEMBER, THIS MODEL ALLOWS RANDOM EFFECTS

## First fit detection function to null model ----
# All models will include yr_sc for trend
gdr_null_hn <- gdistremoval(lambdaformula = ~yr_sc, phiformula = ~1, removalformula = ~1, distanceformula = ~1, data = umf_gdr, keyfun = "halfnorm", output = "density", unitsOut = "ha", mixture = "P", se = FALSE) 

gdr_null_exp <- gdistremoval(lambdaformula = ~yr_sc, phiformula = ~1, removalformula = ~1, distanceformula = ~1, data = umf_gdr, keyfun = "exp", output = "density", unitsOut = "ha", mixture = "P", se = FALSE) 

gdr_null_hz <- gdistremoval(lambdaformula = ~yr_sc, phiformula = ~1, removalformula = ~1, distanceformula = ~1, data = umf_gdr, keyfun = "hazard", output = "density", unitsOut = "ha", mixture = "P", se = FALSE)

gdr_null_unif <- gdistremoval(lambdaformula = ~yr_sc, phiformula = ~1, removalformula = ~1, distanceformula = ~1, data = umf_gdr, keyfun = "uniform", output = "density", unitsOut = "ha", mixture = "P", se = FALSE)

## Model selection for gdistremoval() key function
gdr_keyfun_list <- listN(gdr_null_hn, gdr_null_exp, gdr_null_hz, gdr_null_unif)

View(gdr_keyfun_sel <- FuncFitList(gdr_keyfun_list))

(gdr_keyfun_best <- get(gdr_keyfun_sel$model[1]))

## Add distanceformula covariates for gdistremoval()----
(gdr_p_wind <- update(gdr_keyfun_best, distanceformula = ~weather_wind_comb))

if(use_hab_cov) {
  (gdr_p_hab <- update(gdr_keyfun_best, distanceformula = ~hab_type))
  
  (gdr_p_habcomb <- update(gdr_keyfun_best, distanceformula = ~hab_type_comb))
  
  (gdr_p_hab_wind <- update(gdr_keyfun_best, distanceformula = ~hab_type + weather_wind_comb))
  
  (gdr_p_habcomb_wind <- update(gdr_keyfun_best, distanceformula = ~hab_type_comb + weather_wind_comb))
  
  gdr_p_list <- listN(gdr_keyfun_best, gdr_p_wind, gdr_p_hab, gdr_p_habcomb, gdr_p_hab_wind, gdr_p_habcomb_wind)
} else {
  gdr_p_list <- listN(gdr_keyfun_best, gdr_p_wind)
}

## Model selection for gdistremoval() distanceformula covariates
View(gdr_p_sel <- FuncFitList(gdr_p_list))

(gdr_p_best <- get(gdr_p_sel$model[1]))

## Add removalformula covariates for gdistremoval()----
# I'm calling it phi here but it's actually the removal formula
# If needed for convergence, add: ,control = list(maxit=1e3)
(gdr_phi_julian <- update(gdr_p_best, removalformula = ~julian_prop))

(gdr_phi_julian_time <- update(gdr_p_best, removalformula =~julian_prop + start_frac_time))

(gdr_phi_julian_time2 <- update(gdr_p_best, removalformula =~julian_prop + start_frac_time+ I(start_frac_time^2))) #, starts = rep(0, 10))) #,control = list(maxit=1e3)))

(gdr_phi_time <- update(gdr_p_best, removalformula =~start_frac_time))

(gdr_phi_time2 <- update(gdr_p_best, removalformula =~start_frac_time + I(start_frac_time^2))) #, starts = rep(0, 12))) #,control = list(maxit=1e3)))

## Model selection for gdistremoval() phi covariates
gdr_phi_list <- listN(gdr_keyfun_best, gdr_p_best, gdr_phi_julian, gdr_phi_julian_time, gdr_phi_julian_time2, gdr_phi_time, gdr_phi_time2)
View(gdr_phi_sel <- FuncFitList(gdr_phi_list))

(gdr_phi_best_temp <- get(gdr_phi_sel$model[1]))
gdr_phi_best <- gdr_phi_best_temp

# # >>>>>>>>>> RETRY BEST MODEL WITH any other distanceformula covariates <2 dAIC from best-supported model. Compare AIC's
# gdr_phi_best_temp@formlist$distanceformula
# View(gdr_p_sel)
# 
# gdr_phi_best_temp2 <- update(gdr_phi_best_temp, distanceformula = ~weather_wind_comb)
# 
# temp_list <- listN(gdr_phi_best_temp, gdr_phi_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# 
# gdr_phi_best <- get(temp_sel$model[1]) # <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT. IF CHANGED, CHECK WITH HALF NORMAL KEY FUNCTION AGAIN.
# gdr_phi_best

## Add lambdaformula covariates for gdistremoval()----

if(use_hab_cov) {
  (gdr_abund_yr_hab <- update(gdr_phi_best, lambdaformula = ~yr_sc + hab_type))
  
  (gdr_abund_yr_habcomb <- update(gdr_phi_best, lambdaformula = ~yr_sc + hab_type_comb))
  
  (gdr_abund_yrXhab <- update(gdr_phi_best, lambdaformula = ~yr_sc*hab_type))
  
  (gdr_abund_yrXhabcomb <- update(gdr_phi_best, lambdaformula = ~yr_sc*hab_type_comb))
  
  # Add a full model
  (gdr_full <- gdistremoval(lambdaformula = ~yr_sc*hab_type, phiformula = ~1, removalformula =~julian_prop + start_frac_time+ I(start_frac_time^2), distanceformula = ~hab_type + weather_wind_comb, starts = rep(0, 15),control = list(maxit=1e3), data = umf_gdr, keyfun = "exp", output = "density", unitsOut = "ha", mixture = "P", se = FALSE)) 
  
  # # # IF HABITAT IS A DETECTION COV, ADD A MODEL THAT MOVES IT TO ABUNDANCE COV AND REMOVES FROM DETECTION COV. <<<<<<< Consider adding the model that was the best gdistsamp or gmultmix model (in case it's not included here)--NEED TO ADD IT IN THE FIT LISTS AND SAVE LIST THEN.
  # gdr_best_mod_alt <- NULL
  # (gdr_best_mod_alt <- gdistremoval(lambdaformula = ~yr_sc + hab_type, phiformula = ~1, removalformula =~start_frac_time, distanceformula = ~1, starts = rep(0, 8),control = list(maxit=1e3), data = umf_gdr, keyfun = "halfnorm", output = "density", unitsOut = "ha", mixture = "P", se = FALSE))
  # gdr_alt2 <- update(gdr_best_mod_alt, lambdaformula = ~yr_sc*hab_type, starts = rep(0, 10))
  
  # (gdr_alt2 <- gdistremoval(lambdaformula = ~yr_sc + hab_type_comb, phiformula = ~1, removalformula =~julian_prop, distanceformula = ~weather_wind_comb + hab_type, starts = rep(0, 10),control = list(maxit=1e3), data = umf_gdr, keyfun = "exp", output = "density", unitsOut = "ha", mixture = "P", se = FALSE)) 
  
  gdr_abund_list <- listN(gdr_keyfun_best, gdr_p_best, gdr_phi_best_temp, gdr_phi_best, gdr_abund_yr_hab, gdr_abund_yr_habcomb, gdr_abund_yrXhab, gdr_abund_yrXhabcomb, gdr_full) #, gdr_best_mod_alt, gdr_alt2)
  
} else {
  
  # Add a full model
  (gdr_full <- gdistremoval(lambdaformula = ~yr_sc, phiformula = ~1, removalformula =~julian_prop + start_frac_time+ I(start_frac_time^2), distanceformula = ~weather_wind_comb, starts = rep(0, 9),control = list(maxit=1e3), data = umf_gdr, keyfun = "exp", output = "density", unitsOut = "ha", mixture = "P", se = FALSE))
  
  gdr_abund_list <- listN(gdr_keyfun_best, gdr_p_best, gdr_phi_best_temp, gdr_phi_best, gdr_full)#, gdr_alt1, gdr_alt2)
}

## Model selection for gdistremoval() lambda covariates
View(gdr_abund_sel <- FuncFitList(gdr_abund_list))
(gdr_abund_best_temp <- get(gdr_abund_sel$model[1]))
gdr_abund_best <- gdr_abund_best_temp

# # >>>>>>>>>> RETRY BEST MODEL WITH any other distance covariates <2 dAIC from best-supported model. Compare AIC's
# gdr_abund_best_temp@formlist$distanceformula
# View(gdr_p_sel)
# 
# gdr_abund_best_temp2 <- update(gdr_abund_best_temp, distanceformula = ~weather_wind_comb, starts = NULL) #, starts = rep(0, 7),control = list(maxit=1e3))
# gdr_abund_best_temp3 <- update(gdr_abund_best_temp, distanceformula = ~hab_type_comb + weather_wind_comb, starts = rep(0, 8),control = list(maxit=1e3))
# 
# temp_list <- listN(gdr_abund_best_temp, gdr_abund_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# gdr_abund_best <- get(temp_sel$model[1]) # <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT
# 
# # >>>>>>>>>> RETRY BEST MODEL WITH any other phi covariates <2 dAIC from best-supported model. Compare AIC's
# gdr_abund_best_temp@formlist$removalformula
# View(gdr_phi_sel)
# 
# gdr_abund_best_temp2 <- update(gdr_abund_best_temp, phiformula = ~start_frac_time, starts = rep(0, 9),control = list(maxit=1e3))
# 
# temp_list <- listN(gdr_abund_best_temp, gdr_abund_best_temp2)
# View(temp_sel <- FuncFitList(temp_list))
# gdr_abund_best <- get(temp_sel$model[1]) # <<<<<<<<<<<<<< CHANGE TO 'NEW' BEST MODEL, IF DIFFERENT. IF CHANGED, CHECK WITH HALF NORMAL KEY FUNCTION AGAIN.
# # gdr_abund_best <- gdr_abund_best_temp # <<<<<< OTHERWISE

## Final model selection list----
if(use_hab_cov) {
  gdr_final_list <- listN(gdr_null_hn, gdr_null_exp, gdr_null_hz, gdr_null_unif, gdr_p_wind, gdr_p_hab, gdr_p_habcomb, gdr_p_hab_wind, gdr_p_habcomb_wind, gdr_phi_julian, gdr_phi_julian_time, gdr_phi_julian_time2, gdr_phi_time, gdr_phi_time2, gdr_abund_yr_hab, gdr_abund_yr_habcomb, gdr_abund_yrXhab, gdr_abund_yrXhabcomb, gdr_full, gdr_p_best, gdr_phi_best, gdr_abund_best) # >>>>>> gdr_best_mod_alt, 
} else {
  gdr_final_list <- listN(gdr_null_hn, gdr_null_exp, gdr_null_hz, gdr_null_unif, gdr_p_wind, gdr_phi_julian, gdr_phi_julian_time, gdr_phi_julian_time2, gdr_phi_time, gdr_phi_time2, gdr_full) #, gdr_p_best, gdr_phi_best, gdr_abund_best) # gdr_best_mod_alt, 
}

## All-models final model selection for gdistremoval()
# If the interim 'best' models all have duplicate formula in list, can rerun without them
View(gdr_final_sel <- FuncFitList(gdr_final_list))
(gdr_best_mod_temp <- get(gdr_final_sel$model[1]))

# # # Add the model that was the best gdistsamp or gmultmix model (in case it's not included here)--NEED TO ADD IT IN THE FIT LISTS AND SAVE LIST THEN. <<<<<<<<<<<
# gdist_best_mod@formula
# gmult_best_mod@formula
# 
# (gdr_alt1 <- gdistremoval(lambdaformula = ~yr_sc, phiformula = ~1, removalformula =~julian_prop, distanceformula = ~weather_wind_comb, data = umf_gdr, keyfun = "halfnorm", output = "density", unitsOut = "ha", mixture = "P", se = FALSE))
# 
# (gdr_alt2 <- gdistremoval(lambdaformula = ~yr_sc, phiformula = ~1, removalformula =~julian_prop + start_frac_time, distanceformula = ~weather_wind_comb, data = umf_gdr, keyfun = "halfnorm", output = "density", unitsOut = "ha", mixture = "P", se = FALSE))
# 
# # # Here we can also use a LRT to determine if a more complex model is better supported. If p<0.05, the more complex model is supported as being a statistically significant better fit
# LRT(gdr_p_best, gdr_abund_best)



## Final tweaking of best model for gdistremoval()----
# Try best model with mixture "NB"
gdr_best_mod_temp_nb <- update(gdr_best_mod_temp, mixture = "NB", starts = rep(0, 8)) # ,control = list(maxit=1e3)) # NOTE: gdistremoval can do "ZIP"

gdr_best_mod_temp_zip <- update(gdr_best_mod_temp, mixture = "ZIP", starts = rep(0,8)) # ,control = list(maxit=1e3)) 

View(gdr_mixture_sel <-FuncFitList(listN(gdr_best_mod_temp, gdr_best_mod_temp_nb, gdr_best_mod_temp_zip)))

gdr_best_mod_temp2 <- get(gdr_mixture_sel$model[1])

# Determine if need larger K
(starting_k <- gdr_best_mod_temp2@K)

# k*2
gdr_best_mod_temp2_5k <- update(gdr_best_mod_temp2, K = starting_k*5)
data.frame(rbind(coef(gdr_best_mod_temp2), coef(gdr_best_mod_temp2_20k)))

# k*4
gdr_best_mod_temp2_10k <- update(gdr_best_mod_temp2, K = starting_k*10) 

data.frame(rbind(coef(gdr_best_mod_temp2), coef(gdr_best_mod_temp2_5k), coef(gdr_best_mod_temp2_10k))) # ...do more K if needed (remember to change the model name and also the K value)...

# Add SE to final model----
# <<<<<<<<<<<< CHANGE THIS TO THE CORRECT K AS NEEDED
gdr_best_mod <- update(gdr_best_mod_temp2, starts = coef(gdr_best_mod_temp2), se = TRUE) # , control = list(maxit=1e5)

## MODEL NOTES
gdr_best_mod@K

(GDR_MOD_NOTES <- list(
  best_mod = paste0("Best gdistremoval model formula = ", toString(gdr_best_mod@formula[-1])), # <<<<<<<<< CHANGE GOF
  gof = "",
  other = paste0("Best model is ZIP. Estimates stabilized.")
))# Use this for any notes 

# # Also check K and add SE to gdr_best_mod_alt
# gdr_best_mod_alt_2k <- update(gdr_best_mod_alt, K = starting_k*2)
# gdr_best_mod_alt_4k <- update(gdr_best_mod_alt, K = starting_k*4) 
# 
# data.frame(rbind(coef(gdr_best_mod_alt), coef(gdr_best_mod_alt_2k), coef(gdr_best_mod_alt_4k))) # ...do more K if needed (remember to change the model name and also the K value)...
# 
# # Add SE to final model
# gdr_best_mod_alt <- update(gdr_best_mod_alt, K = starting_k, starts = coef(gdr_best_mod_alt), control = list(maxit=1e5), se = TRUE) # <<<<<<<<<<<< CHANGE THIS TO THE CORRECT K AS NEEDED

## Add random effects >>>>>> FIX LAMBDA FORMULA AS NEEDED
# gdr_best_mod_re <- NULL 
gdr_best_mod@formlist$lambdaformula

(gdr_best_mod_reLoc <- update(gdr_best_mod, lambdaformula = as.formula(paste0(deparse(gdr_best_mod@formlist$lambdaformula), "+ (1|location_name)")), starts = c(coef(gdr_best_mod), 0), control = list(maxit=1e3))) # , control = list(maxit=1e3), se = TRUE))

# (gdr_best_mod_reYr <- update(gdr_best_mod, lambdaformula = as.formula(paste0(deparse(gdr_best_mod@formlist$lambdaformula), "+ (1|yr_f)")), starts = c(coef(gdr_best_mod),8), control = list(maxit=1e3), se = TRUE))

# AIC, just out of curiosity. The one to use should include both year and location RE's because that is how the actual data are collected. BUT ranef cannot support models with >1 RE, so if one of the RE's has variance close to zero, use the other one
# gdr_re_list <- listN(gdr_best_mod, gdr_best_mod_reLoc, gdr_best_mod_reYr)
# View(gdr_re_sel <- FuncFitList(gdr_re_list))


# # Compare estimates
# re_coef_df <- data.frame(rbind(coef(gdr_best_mod), coef(gdr_best_mod_reLoc), coef(gdr_best_mod_reYr), coef(gdr_best_mod_reLocYr)))
# rownames(re_coef_df) <- c("gdr_best_mod", "gdr_best_mod_reLoc", "gdr_best_mod_reYr", "gdr_best_mod_reLocYr")
# View(re_coef_df)
# 
# # Compare CI's
# confint(gdr_best_mod, type = "lambda") # can also use method = "profile" (default is "normal")
# confint(gdr_best_mod_reLoc, type = "lambda")[1:2,]
# confint(gdr_best_mod_reYr, type = "lambda")[1:2,]
# confint(gdr_best_mod_reLocYr, type = "lambda")[1:2,]
# 
# confint(gdr_best_mod, type = "phi")
# confint(gdr_best_mod_reLoc, type = "phi")
# confint(gdr_best_mod_reYr, type = "phi")
# confint(gdr_best_mod_reLocYr, type = "phi")
# 
# confint(gdr_best_mod, type = "rem")
# confint(gdr_best_mod_reLoc, type = "rem")
# confint(gdr_best_mod_reYr, type = "rem")
# confint(gdr_best_mod_reLocYr, type = "rem")
# 
# confint(gdr_best_mod, type = "dist")
# confint(gdr_best_mod_reLoc, type = "dist")
# confint(gdr_best_mod_reYr, type = "dist")
# confint(gdr_best_mod_reLocYr, type = "dist")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## If the RE's worked, take the most appropriate RE model (typically should be gdr_best_mod_reLoc) and check if need higher K
(gdr_best_mod_re <- gdr_best_mod_reLoc)
# gdr_best_mod@AIC
# gdr_best_mod_re@AIC # this should be smaller

(starting_k <- gdr_best_mod_re@K)

# k * 2
gdr_best_mod_re <- update(gdr_best_mod_re, K = starting_k*5) # <<<<<<<<<<<<<<
data.frame(rbind(coef(gdr_best_mod_re), coef(gdr_best_mod_re_2k)))
#                  
# # k * 4
# gdr_best_mod_re_4k <- update(gdr_best_mod_re, K = starting_k*4) 
# data.frame(rbind(coef(gdr_best_mod_re), coef(gdr_best_mod_re_2k), coef(gdr_best_mod_re_4k))) # ...do more K if needed (remember to change the model name and also the K value)...


### Examine final model set----
# data.frame(rbind(coef(gdr_best_mod), coef(gdr_best_mod_re), coef(gdr_best_mod_alt))) # >>>>>>>>>>>>>>>> Remember that year is scaled, so the estimated trend is change per sd of year covariate, NOT change per year


# ## Check model fit  <<<<<<<<<<<< PICK UP FROM HERE
# gdr_best_mod_nm <- update(gdr_best_mod, method='Nelder-Mead')
# (gdr_best_mod_pb <- parboot(gdr_best_mod_nm, fitstats2, nsim = 30, report = 1, ncores = 4)) # switching to `Nelder-Mead` allows this to work better. This is SLOW! For final analyses, use much larger nsim
# 
# gdr_best_mod_re_nm <- update(gdr_best_mod_re, method='Nelder-Mead')
# (gdr_best_mod_re_pb <- parboot(gdr_best_mod2_nm, fitstats2, nsim = 30, report = 1, ncores = 4)) # switching to `Nelder-Mead` allows this to work better. This is SLOW! For final analyses, use much larger nsim


# <<<<<<<<<<<<<<< FROM HERE GO BACK AGAIN TO DO THE FINAL MODEL TWEAKS!

## PLOT MODEL OUTPUTS----
## Create data frame of model predictions for specific covariate combinations
gdr_best_mod@formlist # <<<<< SEE WHAT COVARIATES TO INCLUDE

df_covs <- data.frame(
  julian_prop = mean(umf_gdr@yearlySiteCovs$julian_prop, na.rm = TRUE),
  start_frac_time = mean(umf_gdr@yearlySiteCovs$start_frac_time, na.rm = TRUE), 
  weather_wind_comb = "calm")# this is the most common wind level

df_habs <- unique(umf_gdr@siteCovs[c("location_name", "hab_type", "hab_type_comb")])
 
df_yr <- unique(umf_gdr@siteCovs[c("yr_sc", "yr", "yr_f")]) %>% arrange(yr_sc) # this is so we can make predictions by yr_sc, but see the corresponding yr

(newdata <- merge(df_covs, df_yr))
(newdata <- merge(newdata, df_habs))

## Plot model predictions for specific covariate combinations
# Best mod
predict_out <- FuncPlotPredict(func_name = "gdistremoval", mod = gdr_best_mod, newdata = newdata, is_gdr = TRUE, incl_re = FALSE)
gdr_predict_dat <- predict_out$predict_dat
gdr_predict_plot <- predict_out$p
ggplotly(gdr_predict_plot)

# Best mod RE
predict_out_re <- FuncPlotPredict(func_name = "gdistremoval", mod = gdr_best_mod_re, newdata = newdata, is_gdr = TRUE, incl_re = TRUE)
gdr_predict_dat_re <- predict_out_re$predict_dat
gdr_predict_plot_re <- predict_out_re$p
ggplotly(gdr_predict_plot_re)
# gdr_predict_dat_re <- NULL
# gdr_predict_plot_re <- NULL

## Plot bup densities per site-year 
# These are the estimated abundances per site-year. Differs from predict() in that it incorporates information about the other covariates and observed counts? Note that predict() gives output as predicted density, if model is set to output density
# Ranef gives the empirical Bayes estimates of the superpopulation size. Gives similar mean to predict() but it seems to also consider the actual observed data. >>>> I think I can convert to density just by dividing by survey point area, but the documentation says it's the estimated superpopulation size.

out <- FuncPlotBUPGDR(func_name = "gdistremoval", mod = gdr_best_mod, dat = dat_gremoval)
(gdr_bup_plot <- out$p)
gdr_bup_dat <- out$bup_dat

out_re <- FuncPlotBUPGDR(func_name = "gdistremoval", mod = gdr_best_mod_re, dat = dat_gremoval)
(gdr_bup_plot_re <- out_re$p)
gdr_bup_dat_re <- out_re$bup_dat
# gdr_bup_dat_re <- NULL
# gdr_bup_plot_re <- NULL

out_alt <- FuncPlotBUPGDR(func_name = "gdistremoval", mod = gdr_best_mod_alt, dat = dat_gremoval)
(gdr_bup_plot_alt <- out_alt$p)
gdr_bup_dat_re <- out_re$bup_dat

### SAVE ALL MODELS----
# >>>>> ADD PB WHEN I CAN GET GOF TO WORK
if(use_hab_cov) {
  save(spec, unit, max_bin_bound, combine_bins, max_time, pair_time_int, single_visit, dat_gmix, dat_gdist, dat_gremoval, umf_gdr, umf_gds, umf_gmm, gdr_null_hn, gdr_null_exp, gdr_null_hz, gdr_null_unif, gdr_p_wind, gdr_p_hab, gdr_p_habcomb, gdr_p_hab_wind, gdr_p_habcomb_wind, gdr_phi_julian, gdr_phi_julian_time, gdr_phi_julian_time2, gdr_phi_time, gdr_phi_time2, gdr_abund_yr_hab, gdr_abund_yr_habcomb, gdr_abund_yrXhab, gdr_abund_yrXhabcomb, gdr_full, gdr_best_mod_alt, gdr_p_best, gdr_phi_best, gdr_abund_best, gdr_final_sel, gdr_best_mod, gdr_best_mod_re, gdr_predict_plot, gdr_predict_plot_re, gdr_predict_dat, gdr_predict_dat_re, gdr_bup_plot, gdr_bup_dat, gdr_bup_plot_re, gdr_bup_dat_re, GDR_MOD_NOTES, file = here::here("Model_fits", "gdistremoval", paste0("gdistremoval_", save_prefix, ".RData")))
} else {
  save(spec, unit, max_bin_bound, combine_bins, max_time, pair_time_int, single_visit, dat_gmix, dat_gdist, dat_gremoval, umf_gdr, umf_gds, umf_gmm, gdr_null_hn, gdr_null_exp, gdr_null_hz, gdr_null_unif, gdr_p_wind, gdr_phi_julian, gdr_phi_julian_time, gdr_phi_julian_time2, gdr_phi_time, gdr_phi_time2, gdr_full, gdr_p_best, gdr_phi_best, gdr_abund_best, gdr_final_sel, gdr_best_mod, gdr_best_mod_re, gdr_predict_plot, gdr_predict_plot_re, gdr_predict_dat, gdr_predict_dat_re, gdr_bup_plot, gdr_bup_dat, gdr_bup_plot_re, gdr_bup_dat_re, GDR_MOD_NOTES,file = here::here("Model_fits", "gdistremoval", paste0("gdistremoval_", save_prefix, ".RData"))) # gdr_best_mod_alt, 
}

### RESULTS NOTES----

## These DISTANCE models should not be used b/c failed GOF or estimates kept increasing w/K
# > AMCR (best model was NB and estimates kept changing with K)

## These REMOVAL models should not be used b/c failed GOF or estimates kept increasing w/K
# > AMCR (best model was NB and estimates kept changing with K)

## These DISTREMOVAL models should not be used b/c failed GOF or estimates kept increasing w/K
# > AMCR (best model was NB and estimates kept changing with K)

