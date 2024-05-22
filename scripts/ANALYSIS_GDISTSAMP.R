### N-MIXTURE GDISTSAMP in unmarked

## >>>>>>>>> SEE ?RANEF AND ?PARBOOT
## >>>>>>>> REMEMBER THAT FOR FINAL ANALYSES NEED TO USE BOOTSTRAP ON GDISTSAMP MODEL OUTPUTS B/C COULD NOT INCLUDE RANDOM EFFECTS!

### Unmarked gdistsamp()
# https://groups.google.com/g/unmarked/c/hcGgk6VxHaA
# https://darinjmcneil.weebly.com/hierarchical-distance-models.html
# https://groups.google.com/g/unmarked/c/OesoTLrawgE
# https://groups.google.com/g/unmarked/c/Wvsne3_jx1U/m/R6vNPLvNAAAJ?utm_medium=email&utm_source=footer

### Follow this one for post-model diagnostics and predictions
# https://cran.r-project.org/web/packages/unmarked/vignettes/distsamp.html

### Suggested
# Year-specific phi. Multiply to lambda (superpopulation) to get "availability" per year.
# To get SE: Nonparametric bootstrap- repeatedly (e.g., 1000 times) draw a sample (with replacement) of n sites out of your real data set containing n sites, fit the model and estimate phi and lambda, multiply them and save that product. The SD of the 1000 products will then be the SE of the desired quantity
# The 'g' in gdistsamp means calculating phi, because not all animals are available in all surveys
# Use of the hazard-rate key function ("hazard") typically requires a large sample size in order to get good parameter estimates. If you have a relatively small number of points/transects (<100), you should be cautious with the resulting models. Check your results against estimates from unmarked, which doesn't require as much data to get good estimates of the hazard-rate shape and scale parameters.

rm(list = ls())

### Load packages ----
library(tidyverse)
library(magrittr)
library(unmarked)
library(plotly)

### FUNCTIONS ----
FuncEstim <- function(mod) {
  if(!single_visit){print(backTransform(mod, type="phi"))}
  print(backTransform(mod, type="det"))
  
  # Empirical Bayes estimates of abundance at each site
  print(ranef(mod))
}

### READ DATA/SCRIPTS ----
# # These are the ones I'll fit models to
# common_spec <- read_csv("C:/Users/echeng/Desktop/GULN_birdtrends_working copy/Data_out/common_species_df.csv") %>% 
#   dplyr::filter(species_code %in% c("BLJA", "CARW", "GTGR", "NOCA", "WEVI")) %>%
#   dplyr::filter(unit_code != "GUIS-MS")
# common_spec$unit_code[common_spec$unit_code == "GUIS-FL"] <- "GUIS"


source(here::here("scripts", "ANALYSIS_STARTUP.R"))

dist_breaks_km <- c(0, 0.025, 0.05, 0.1)

### USER ENTRY ----

park <- "SAAN" # <<<<<<<<<< PICK UP FROM HERE
spec <- "NOCA"
yearmod <- TRUE
max_bin_bound <- 100
single_visit <- TRUE

tmb_trend<- readRDS("C:/Users/echeng/Desktop/GULN_birdtrends_working copy/tmb_mods/tmb_BITH_NOCA.RDS")
tmb_trend$call # This is the best-fit model for this park-species, from tmb
summary(tmb_trend)

(save_prefix <- paste0(park, "_", spec, ifelse(single_visit, "_1visit", "_2visit"), ifelse(yearmod, "_yearmod", "")))

### Format the data ----
dat <- FuncFormatUnmarkedDist(park = park, spec = spec, single_visit = single_visit)

table(dat$researcher, dat$yr, dat$subunit) # check how many observers
# BITH HAS 3 RESEARCHERS; PAAL HAS 2; JELA & VICK have 1

## Dataset with complete cases only
# PAAL: Dropped 3 records b/c missing weather_wind_comb
# SAAN: Dropped 1 record b/c missing hrs_since_rise
nrow(dat)
colSums(is.na(dat))

dat <- dat[complete.cases(dat),]
nrow(dat)


### N-MIXTURE ANALYSIS ----
table(dat$yr_visit) # make sure nothing funny. For BITH, there is one location sampled twice in 2018 but no reason to drop it.

dat_Nmix <- dat

### Format data for 'ubms' ----

(primary_surv <- sort(unique(na.omit(dat_Nmix$yr)))) # these are the years of survey

(secondary_surv <- 1:max(as.numeric(sub(".*_", "", dat_Nmix$yr_visit)), na.rm = TRUE)) # sequence from one to maximum number of repeat survey visits (over entire dataset)-- This should just be 2 visits.

template <- expand_grid(location_name = sort(unique(na.omit(dat_Nmix$location_name))), yr = primary_surv, visit = secondary_surv) %>%
  dplyr::mutate(
    loc_yr = paste(location_name, yr, sep = "_"),
    yr_visit = paste(yr, visit, sep = "_")) %>%
  dplyr::select(-yr)
head(template)
table(template$location_name, template$yr_visit) # Check for anything funny

expand_dat <- template %>%
  left_join(dat_Nmix, by = c("location_name", "yr_visit"))
glimpse(expand_dat)

table(expand_dat$yr_visit) # check that every year has the same number of repeat visits and the same number of locs. The additional records are NA, indicating missing data.

### Response variable 
temp_yDistance <- FuncMatrix_Nmix(dat = expand_dat, num_surv = max(secondary_surv), cov_vec = c("distbin_0", "distbin_1", "distbin_2"), response_mat = TRUE, stack = TRUE, incl_site_yr = TRUE, gdistsamp = TRUE) # convert it to stacked format

dim(temp_yDistance)
# View(temp_yDistance)
table(temp_yDistance$location_name)
table(temp_yDistance$loc_yr)

### Site covariates
### FOR GUIS AND SAAN INCL SUBUNIT!!! <<<<

site_covs <- dat_Nmix %>%
  dplyr::select(subunit, location_name, hab_type_100_comb, hab_type_200_comb, prop_understory_cover_50, prop_understory_cover_100, prop_understory_cover_200, understory_cover_sd_50, understory_cover_sd_100, understory_cover_sd_200) %>%
  dplyr::distinct() %>%
  dplyr::arrange(location_name)

table(site_covs$hab_type_100_comb)
table(site_covs$hab_type_200_comb) # For BITH there is only forest; VICK only has 5 non-forest; SAAN only has 4 non-forest
length(unique(site_covs$location_name)) == nrow(site_covs) # Make sure this is true, that no location names are duplicated

site_covs_st <- do.call("rbind", replicate(
  length(unique(na.omit(dat_Nmix$yr))), site_covs, simplify = FALSE)) # stack the data, one rep per year
site_covs_st$yr <- rep(sort(unique(na.omit(dat_Nmix$yr))), each = length(unique(na.omit(dat_Nmix$location_name)))) # for year trend
site_covs_st$yr_c <- scale(site_covs_st$yr, center = TRUE, scale = FALSE) %>% as.numeric() 
site_covs_st$yr_sc <- scale(site_covs_st$yr) %>% as.numeric()
site_covs_st$yr_f <- as.factor(site_covs_st$yr) # for year random effect
site_covs_st %<>% dplyr::mutate(loc_yr = paste0(location_name, "_", yr))


site_covs_st %<>% dplyr::arrange(location_name, loc_yr)

### IF ONLY ONE SURVEY/YR, INCLUDE 'JULIAN_PROP_C' AND 'HRS_SINCE_RISE_SC' B/C THIS IS GOING TO BE A COVARIATE ON ABUNDANCE, NOT DETECTION THEN
if(single_visit) {
  uniq_julian <- dat_Nmix %>%
    dplyr::select(location_name, yr, julian_prop_c) %>%
    dplyr::distinct()
  
  site_covs_st %<>%
    dplyr::left_join(uniq_julian, by = c("location_name", "yr"))
  
  uniq_hrs_rise <- dat_Nmix %>%
    dplyr::select(location_name, yr, hrs_since_rise_sc) %>%
    dplyr::distinct()
  
  site_covs_st %<>%
    dplyr::left_join(uniq_hrs_rise, by = c("location_name", "yr"))
  
}

dim(site_covs_st)
glimpse(site_covs_st)

all(temp_yDistance$loc_yr == site_covs_st$loc_yr) # Make sure all rows are in same loc_yr order as yDistance!!

# If order of records is correct, then convert to matrix with just the response
yDistance <- as.matrix(temp_yDistance %>% dplyr::select(-location_name, -loc_yr))

siteCovs <- site_covs_st

### `Observation` covariates (julian date, start time, and wind)--should have MXT rows, where T is the number of site visits 
glimpse(expand_dat) # check it 
unique(expand_dat$researcher) # Are there multiple researchers? <<<<<<<<< If yes, include as covariate

yearlySiteCovs <- expand_dat %>%
  dplyr::select(researcher, first_yr, location_name, loc_yr, yr_f, visit, julian_prop_c, hrs_since_rise_sc, weather_wind_comb) %>% # Add 'yr_f' if doing yearmod; otherwise add something for trend
  dplyr::arrange(loc_yr, visit) # site-yr major, visit minor order

dim(yearlySiteCovs)
head(yearlySiteCovs)

# For gdistsamp() distance sampling
umf_gds <- unmarkedFrameGDS(
  y = yDistance,
  siteCovs = siteCovs,
  numPrimary = max(yearlySiteCovs$visit, na.rm = TRUE),
  yearlySiteCovs = yearlySiteCovs,
  dist.breaks = dist_breaks_km,
  survey = "point",
  unitsIn = "km")
summary(umf_gds)
View(umf_gds)
# >>>>>>> PICK FUP FROM HERE 
# These park units have multiple researchers -- BITH, PAAL, SAAN
# These park units have only one habitat type -- BITH, PAAL
# JELA & VICK & SAAN each has only 5 or fewer non-forest habitats

cat("STARTING GDISTSAMP MODEL---------------------")

### IF JUST DOING THE 2VISIT, DO A SHORTCUT ----
in_list <- readRDS(here::here("gdistsamp_mods", paste0("gdistsamp_", park, "_", spec, "_1visit_yearmod.RDS"))) 

in_list$best_mod@call # >>>>> COPY IT BELOW, BUT MOVE JULIAN_PROP AND HRS_SINCE_RISE FROM LAMBDAFORMULA TO PHIFORMULA (AND LEAVE IN PFORMULA IF THERE)
mod1 <- gdistsamp(lambdaformula = ~yr_f + subunit + yr_f:subunit, phiformula = ~1, 
                  pformula = ~1, data = umf_gds, keyfun = "halfnorm", output = "abund", 
                  mixture = "P", se = TRUE)
summary(mod1)

mod2 <- gdistsamp(lambdaformula = ~yr_f + subunit + yr_f:subunit, phiformula = ~julian_prop_c + hrs_since_rise_sc + I(hrs_since_rise_sc^2), 
                  pformula = ~1, data = umf_gds, keyfun = "halfnorm", output = "abund", 
                  mixture = "P", se = TRUE)

mod3 <- gdistsamp(lambdaformula = ~yr_f + subunit + yr_f:subunit, phiformula = ~hrs_since_rise_sc + I(hrs_since_rise_sc^2), 
                  pformula = ~1, data = umf_gds, keyfun = "halfnorm", output = "abund", 
                  mixture = "P", se = TRUE)

mod4 <- gdistsamp(lambdaformula = ~yr_f + subunit + yr_f:subunit, phiformula = ~julian_prop_c, 
                  pformula = ~1, data = umf_gds, keyfun = "halfnorm", output = "abund", 
                  mixture = "P", se = TRUE)

# take out any unimportant phiformula covs

# If hrs_since_rise in pformula, take it out


## FINAL Model selection
final_mods_list <- listN(mod1, mod2, mod3, mod4)
(final_mods_table <- unmarked::fitList(fits=final_mods_list) %>% unmarked::modSel(.))

gdist_best_mod <- mod3 # If a more complex covariate model is within approx. 2AIC of best-fit model, use it.

# Then also add a more complex model just to see how different estimates are
gdist_complex_mod <- mod4

# >>>>>>> NOW CHECK MODEL FIT
out_list <- list(
  park = park,
  spec = spec,
  yearmod = yearmod,
  single_visit = single_visit,
  umf_gds = umf_gds,
  final_mods = final_mods_list,
  final_mods_table = final_mods_table,
  best_mod = gdist_best_mod,
  complex_mod = gdist_complex_mod,
  gdist_best_mod_pb = gdist_best_mod_pb,
  gdist_complex_mod_pb = gdist_complex_mod_pb,
  gdist_best_mod_pred = gdist_best_mod_pred,
  gdist_complex_mod_pred = gdist_complex_mod_pred,
  GDIST_MOD_NOTES = GDIST_MOD_NOTES
)
View(out_list)
saveRDS(out_list, here::here("gdistsamp_mods", paste0("gdistsamp_", save_prefix, ".RDS"))) 

### FIND BEST-FIT MODEL WITH GDISTSAMP()----
## First fit gdistsamp() key function----
# For test_mode, use yr_f; otherwise, use yr_sc
 gdist_null_list <- list()
for (distfun in c("halfnorm", "exp", "hazard")) {
  for (mix in c("P", "NB")) {
  gdist_null_list[[paste(distfun, mix, sep = "_")]] <- tryCatch(gdistsamp(lambdaformula = ~yr_f, phiformula = ~1, pformula = ~1, data = umf_gds, keyfun = distfun, output = "abund", mixture = mix, se = FALSE), error = function(e) {cat("error: ", distfun)})
  }
}

## Model selection for gdistsamp() key function
null_mods <- unmarked::fitList(fits=gdist_null_list)
(null_mods_table <- unmarked::modSel(null_mods))

gdist_null_best <- gdistsamp(lambdaformula = ~yr_f, phiformula = ~1, pformula = ~1, data = umf_gds, keyfun = "exp", output = "abund", mixture = "NB", se = FALSE) # <<<<<<<< ENTER THE BEST MODEL KEYFUN AND MIXTURE. I can't simply read the list element corresponding to the best_distshape b/c causes problems when I try to update from this model

### Add pformula covariates ----
# Add pformula = ~..., this is probability of detecting when a bird calls. Should include researcher, first_yr, wind, and habitat covariates. Also include hrs_since_rise. Possibly include julian_prop here IF only single survey (otherwise belongs in phiformula = ~...)
# May need to add this if model has problems: starts = rep(0, XXX),control = list(maxit=1e3)))

gdist_pfull <- update(gdist_null_best, pformula = ~weather_wind_comb + hab_type_200_comb + prop_understory_cover_200, se = TRUE, control = list(maxit=1e3)) # First p-model includes all likely p-covariates (anything with p> 0.5 from tmb trend model)
summary(gdist_pfull)

gdist_p1 <- update(gdist_null_best, pformula = ~ hab_type_200_comb + prop_understory_cover_200, se = TRUE, control = list(maxit=1e3)) # Then drop the less statistically significant ones (drop anything with p> 0.5)
summary(gdist_p1)

# # Have a model that drops julian_prop and hr_since_rise
gdist_p2 <- update(gdist_null_best, pformula = ~researcher + prop_understory_cover_50, se = TRUE, control = list(maxit=1e3)) # Then drop the less statistically significant ones (drop anything with p> 0.5)
summary(gdist_p2)

# # Also have a model with only researcher + first_yr, and another with only researcher
gdist_p3 <- update(gdist_null_best, pformula = ~researcher + first_yr, se = TRUE, control = list(maxit=1e3)) # Then drop the less statistically significant ones (drop anything with p> 0.5)
summary(gdist_p3)

gdist_p4 <- update(gdist_null_best, pformula = ~researcher, se = TRUE, control = list(maxit=1e3)) # Then drop the less statistically significant ones (drop anything with p> 0.5)
summary(gdist_p4)

## Model selection
(p_mods_table <- unmarked::fitList(fits=listN(gdist_null_best, gdist_pfull, gdist_p1)) %>% unmarked::modSel(.))

gdist_p_best <- gdist_p1

# If a more complex covariate model is within approx. 2AIC of best-fit model, use it.

### Plot distance function ----
# >>>> NOW FOR SINGLE VISIT, PLOT THE DETECTION FUNCTION USING THE BEST NON-NULL MODEL IF CLOSE TO NULL <<<
mod <- gdist_p_best # If important categorical covariate on detection, use that; otherwise, use the null model
summary(mod)

### If single visit, phibest = pbest ----
gdist_phi_best <- gdist_p_best

# ### Add phiformula covariates, if multiple visits per year ----
# gdist_phifull <- update(gdist_p_best, phiformula = ~julian_prop_c + hrs_since_rise_sc + I(hrs_since_rise_sc^2), se = TRUE)
# summary(gdist_phifull)
# 
# gdist_phi1 <- update(gdist_p_best, phiformula = ~julian_prop_c + hrs_since_rise_sc, se = TRUE)
# summary(gdist_phi1)
# 
# gdist_phi2 <- update(gdist_p_best, phiformula = ~julian_prop_c, se = TRUE)
# summary(gdist_phi2)
# 
# ## Model selection
# (phi_mods_table <- unmarked::fitList(fits=listN(gdist_null_best, gdist_p_best, gdist_phifull, gdist_phi1, gdist_phi2)) %>% unmarked::modSel(.))
# 
# gdist_phi_best <- gdist_phi1

### Add lambdaformula covariates ----
# This should be year-factor, julian_prop, hrs_since_rise, and habitat variables
gdist_lambdafull <- update(gdist_phi_best, lambdaformula = ~yr_f + hab_type_200_comb + julian_prop_c + prop_understory_cover_200, se = TRUE) # First model includes all likely abundance covariates (anything with p> 0.5 from tmb trend model)
summary(gdist_lambdafull)

gdist_lambda1 <- update(gdist_phi_best, lambdaformula = ~yr_f + julian_prop_c + prop_understory_cover_200, se = TRUE) # First model includes all likely abundance covariates (anything with p> 0.5 from tmb trend model)
summary(gdist_lambda1)

# Also add a full(ish) model ----
# just add the lambdafull to the pfull/phifull
gdist_full <- update(gdist_pfull, lambdaformula = ~yr_f + hab_type_200_comb + julian_prop_c + prop_understory_cover_200, se = TRUE)
summary(gdist_full)

## Model selection so far...
unmarked::fitList(fits=listN(gdist_null_best, gdist_p_best, gdist_lambdafull, gdist_lambda1, gdist_lambda2, gdist_lambda3, gdist_full)) %>% unmarked::modSel(.)

# Finally, IF NOT ALREADY IN ABOVE MODEL SET, add models putting julian_prop and hrs_since_rise in the "most appropriate" model component--detection or abundance

## FINAL Model selection
final_mods_list <- listN(gdist_null_best, gdist_p_best, gdist_lambdafull, gdist_lambda1, gdist_full)
(final_mods_table <- unmarked::fitList(fits=final_mods_list) %>% unmarked::modSel(.))

gdist_best_mod <- gdist_lambda1 # If a more complex covariate model is within approx. 2AIC of best-fit model, use it.

# Then also add a more complex model just to see how different estimates are
gdist_complex_mod <- gdist_lambdafull
  
# # Here we can also use a LRT to determine if a more complex model is better supported. If p<0.05, the more complex (second) model is supported.
LRT(gdist_null_best, gdist_p_best)
LRT(gdist_best_mod, gdist_complex_mod)

# Determine if need larger K
(starting_k <- gdist_best_mod@K)

# k*2
gdist_best_mod_2k <- update(gdist_best_mod, K = starting_k*2)
data.frame(rbind(coef(gdist_best_mod), coef(gdist_best_mod_2k))) # if results are the same, then starting K is sufficient

# Add SE to final model, if not already there. Can use 'starts = coef(gdist_best_mod)'


### Examine final model----
summary(gdist_best_mod) # Remember that IF year is scaled, the estimated trend is change per sd of year covariate, NOT change per year

### Check model fit ----
gdist_best_mod_nm <- update(gdist_best_mod, method='Nelder-Mead')
(gdist_best_mod_pb <- parboot(gdist_best_mod, fitstats2, nsim = 80, report = 1, ncores = 20))
plot(gdist_best_mod_pb) # Nice plots for GOF

# Also try with the complex model, if seems promising
gdist_complex_mod_nm <- update(gdist_complex_mod, method='Nelder-Mead')
(gdist_complex_mod_pb <- parboot(gdist_complex_mod_nm, fitstats2, nsim = 80, report = 1, ncores = 20))
plot(gdist_complex_mod_pb)

# ### Bootstrap abundance in surveyed points
# # Population size in sampled area
# mod <- gdist_best_mod
# mod <- gdist_complex_mod
# 
# Nhat <- function(mod) {
#   sum(bup(ranef(mod)))
# }
# 
# (pb_best_mod_total <- parboot(gdist_best_mod, Nhat, nsim=80, report=5))
# colSums(confint(ranef(gdist_best_mod))) # compare to empirical Bayes confidence intervals
# 
# (pb_complex_mod_total <- parboot(gdist_complex_mod, Nhat, nsim=80, report=5))
# colSums(confint(ranef(gdist_complex_mod))) # compare to empirical Bayes confidence intervals

### MOD NOTES ----
(GDIST_MOD_NOTES <- list(
  best_mod = paste0("Best gdistsamp model key function = ", gdist_best_mod@keyfun, "; formula = ", toString(gdist_best_mod@formula[-1])),
  gof = "Good", 
  other = "Distance function fit is good for forest habitat (which is most of the data), not good for non-forest. Best model doesn't include habitat type anyway"))# Use this for any notes 

### PLOT MODEL OUTPUTS  ----
## Year effects plots
# We can get simple covariate effects plots (i.e., marginal effects plots) using the `plotEffects` function.
# This function only allows one covariate at a time to be varying, while the other covariates are held at their median values or reference levels.
gdist_best_mod_pred <- plotEffectsData(gdist_best_mod, 'lambda', covariate='yr_f')
plotEffects(gdist_best_mod, 'lambda', covariate='yr_f')

gdist_complex_mod_pred <- plotEffectsData(gdist_complex_mod, 'lambda', covariate='yr_f')
plotEffects(gdist_complex_mod, 'lambda', covariate='yr_f')

# Alternatively, if we want to specify ourselves, or for interactions in plots...

mod = gdist_best_mod
mod = gdist_complex_mod
nd = merge(
  data.frame(yr_f = sort(unique(mod@data@siteCovs$yr_f))),
  data.frame(subunit = c("GUIS-FL", "GUIS-MS")))
nd$hrs_since_rise_sc = 0
                
nd = data.frame(yr_sc = sort(unique(mod@data@siteCovs$yr_sc)),
                yr_f = sort(unique(mod@data@siteCovs$yr_f)),
                first_yr = 0,
                researcher = "Swanson, Erin", # <<<<<<<<<< BITH
                # researcher = "Pruitt,  Kenneth", # <<<< PAAL
                weather_wind_comb = "calm",
                hab_type_100_comb = "forest", # <<<<<<<<<
                hab_type_200_comb = "forest", # <<<<<<<<<
                julian_prop_c = 0,
                hrs_since_rise_sc = 0,
                understory_cover_sd_50 = mean(mod@data@siteCovs$understory_cover_sd_50),
                understory_cover_sd_200 = mean(mod@data@siteCovs$understory_cover_sd_200),
                understory_cover_sd_100 = mean(mod@data@siteCovs$understory_cover_sd_100),
                prop_understory_cover_200 = mean(mod@data@siteCovs$prop_understory_cover_200),
                prop_understory_cover_100 = mean(mod@data@siteCovs$prop_understory_cover_100),
                prop_understory_cover_50 = mean(mod@data@siteCovs$prop_understory_cover_50))
nd

# Get yearmod predictions
gdist_best_mod_pred <- unmarked::predict(mod, type = "lambda", newdata = nd, re.form = NA) %>%
  janitor::clean_names()

gdist_complex_mod_pred <- unmarked::predict(mod, type = "lambda", newdata = nd, re.form = NA) %>%
  janitor::clean_names()


### SAVE RESULTS ----

out_list <- list(
  park = park,
  spec = spec,
  yearmod = yearmod,
  single_visit = single_visit,
  umf_gds = umf_gds,
  final_mods = final_mods_list,
  final_mods_table = final_mods_table,
  best_mod = gdist_best_mod,
  complex_mod = gdist_complex_mod,
  gdist_best_mod_pb = gdist_best_mod_pb,
  gdist_complex_mod_pb = gdist_complex_mod_pb,
  gdist_best_mod_pred = gdist_best_mod_pred,
  gdist_complex_mod_pred = gdist_complex_mod_pred,
  GDIST_MOD_NOTES = GDIST_MOD_NOTES
)
View(out_list)
saveRDS(out_list, here::here("gdistsamp_mods", paste0("gdistsamp_", save_prefix, ".RDS"))) 






# ### Inference from top model
# # Get model parameter estimates
# est <- coef(mod_habxyear)
# # Effect of hardwood relative to river levee/bottomland [e^beta]
# hardwood_effect <- exp(est[2])
# # Convert the effect size to a percent increase
# pct_higher <- round((hardwood_effect-1) * 100)
# # Print result
# names(pct_higher) <- "Percent hardwood higher than bottomland"
# pct_higher
# 
# 
# ### Percent yearly decline in river levee/bottomland hardwoods:
# # Since bottomland is the reference level for habitat, this is just
# # The effect size of the parameter for "year"
# year_effect <- exp(est[3])
# # Convert to a percent decline
# pct_decline_bottom <- round((1 - year_effect)*100, 1)
# names(pct_decline_bottom) <- "Percent yearly decline bottomland"
# stopifnot(pct_decline_bottom == 2.7)
# pct_decline_bottom
# 
# ### Percent yearly decline in hardwood plantations:
# # For yearly decline in hardwoods, the effect size is 
# # year + interaction of year and hardwood
# hardwood_year_effect <- exp(est[3] + est[4])
# # Convert to percent decline
# pct_decline_plant <- round((1 - hardwood_year_effect)*100, 1)
# names(pct_decline_plant) <- "Percent yearly decline hardwood plantation"
# stopifnot(pct_decline_plant == 7.0)
# pct_decline_plant
# 
# ### Effects plots
# # We can get simple covariate effects plots (i.e., marginal effects plots) using the `plotEffects` function.
# # This function only allows one covariate at a time to be varying, while the other covariates are held at their median values or reference levels.
# # For example, comparing predicted habitat:
# plotEffects(mod_habxyear, 'lambda', covariate='Habitat')
# 
# # Or changes by year:
# plotEffects(mod_habxyear, 'lambda', covariate='Year')
# 
# # Note that since the habitat covariate is held at its reference level (river levee/bottomland), the plot above represents changes in abundance over time for river levee/bottomland habitat only.
# # You can get the data used to make this plot (e.g. to customize the plot yourself) using `plotEffectsData`:
# plotEffects(gdist_best_mod, 'lambda', covariate='yr_f')
# plotEffectsData(gdist_best_mod, 'lambda', covariate='yr_f')
# 
# ## For this model we have an interaction between habitat type and year, which `plotEffects` cannot currently handle.
# # To show this interaction, we'll need to build the plot manually using `predict` and `ggplot`.
# nd <- data.frame(Habitat=factor(levels(umf@siteCovs$Habitat),
#                                 levels=levels(umf@siteCovs$Habitat)),
#                 Year=rep(0:17, each=2))
# pr <- predict(mod_habxyear, 'lambda', newdata=nd, re.form=NA)
# pr <- cbind(pr, nd)
# 
# ### Build plot:
# ggplot(data=pr, aes(x=Year+2005, y=Predicted)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper, fill=Habitat), alpha=0.2) +
#   geom_line(aes(col=Habitat)) +
#   labs(y="ACFL abundance and 95% CI", x="Year") +
#   theme_bw(base_size=14) +
#   theme(panel.grid=element_blank(), legend.pos=c(0.8,0.8),
#         strip.background=element_rect("white"))
# 
# ### Save plot to file:
# tiff("Figure_2.tiff", height=5, width=7, units='in', res=300, compression='lzw')
# ggplot(data=pr, aes(x=Year+2005, y=Predicted)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper, fill=Habitat), alpha=0.2) +
#   geom_line(aes(col=Habitat)) +
#   labs(y="ACFL abundance and 95% CI", x="Year") +
#   theme_bw(base_size=14) +
#   theme(panel.grid=element_blank(), legend.pos=c(0.8,0.8),
#         strip.background=element_rect("white"))
# dev.off()
# 
# 
# ## Latent Abundance estimates for each site-year
# 
# # The `ranef` function generates posterior distributions of latent abundance or occurrence (abundance in this example) for each site-year.
# # The `bup` function converts these posterior distributions into a point estimate of abundance for each site-year. See `?ranef` and `?bup` for more information.
# r <- ranef(mod_habxyear)
# head(bup(r))
# 
# ### You can apply the `predict` function to the output from `ranef` in order to generate derived posterior distributions - for example, posteriors for the mean latent abundance in each habitat type.
# # See `?predict` for more information on this functionality.
# 
# # Function to calculate the mean across sites of each habitat type
# hab_comp <- function(x){
#   c(river = mean(x[umf@siteCovs$Habitat == "River Levee"]),
#     hardwood = mean(x[umf@siteCovs$Habitat == "Hardwood Plantation"]))
# }
# 
# set.seed(123)
# ### Take samples and calculate mean for each sample, then calculate stats
# out_mat <- predict(r, func=hab_comp)
# hab_mean <- t(apply(out_mat, 1, function(x){
#                       c(mean=mean(x), quantile(x, c(0.025,0.975)))
#                     }
#             ))
# stopifnot(all(round(hab_mean[,1], 2) == c(5.24, 6.06)))
# hab_mean
