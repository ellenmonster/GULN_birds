### GULN BIRDS - glmmTMB

### READINGS ----
# > https://stats.stackexchange.com/questions/86269/what-is-the-effect-of-having-correlated-predictors-in-a-multiple-regression-mode
# > https://stats.stackexchange.com/questions/78828/is-there-a-difference-between-controlling-for-and-ignoring-other-variables-i/78830#78830

### FINAL MODELS ----
# > Should use % forest w/in 100m as covariate, instead of the two-level habitat type. But remember that different parks may define "forest" differently. Also, sometimes only part of the 100-m point buffer lies within park boundary and we only have habitat data for what lies within park boundary, so those habitat estimates (e.g., % forest) may no misleading
# > brms cross-validate for "model selection" (instead of AICc)
# > bootstrap the CI's
# > ORDINATION to identify points that are very "different" from expectations based on habitat? Remember there were a couple sites at one park that I was going to suggest dropping b/c numbers often inexplicably very different there
# > With year-factor models, glmmTMB RE variances seem way too small, and get those strangely huge 95%CI's for years with 0 detections (for Poisson model; then for ZIP, the estimates are NaN for 95% CI's)
# > With year-factor models, rank-deficiency when researcher changes by year (e.g., GUIS had Mark Woodrey 2013-2016, then Adam Rohnke 2017-2022)... why isn't it parsing out that researcher effect?
# > QUESTION... For zero-inflated models, can/should the zero inflation be a function of location (which is a RE)?
# > For each park, write out the pre-analysis data-cleaning rules such as dropping certain researcher, any "odd" locations that should be excluded (perhaps for a subset of sites)--get QE thoughts about when it's okay to just drop a researcher who only surveyed for a year or two of the early or late years
# > For BITH, TC29 is surveyed twice in 2018, so there is just one record that is 2018_2. In final, decide if this should be excluded from analyses?
# > Some glmmTMB RE variance estimates are near zero--these may be wrong! May need to try with a different package

### NOTES ----
# > The habitat variables are treated as if these sites are static, which is certainly not true. TO be appropriate covariates, we really need habitat data to be updated over time.
# > VICK-- only use data from Daniel and only starting in 2012. I dropped 3 VICK surveys from 2010 to keep a consistent max 2 surveys per year for N-mixture. Do the same for GLMM just for ability to compare estimates between approaches. <<<<<
# > PAAL: Dropped 3 records b/c missing weather_wind_comb. Also drop hab_type_100 b/c too many missing. Has multiple researchers, so add as predictor--BUT REALLY should consider dropping the first two years so can just have one researcher [for the low pop 100m models I did drop the first two years--so not comparable to the higher pop models]. PAAL only has non-forest habitat.
# > BITH: Has multiple researchers, so add as predictor AND REMOVE FIRST THREE YEARS B/C RESEARCHER CHANGES CORRELATED WITH YEARS. BITH only has forest habitat.
# > SAAN: Dropped 1 record b/c missing hrs_since_rise. Add subunits: SAAN1-6 ARE 'RA', AND REST ARE 'MI' SUBUNIT and yr*subunit interaction. SAAN-RA has no surveys in 2013??
# > GUIS: Add subunits and yr*subunit interaction. GUIS =-MS only has forest habitat, so rank deficient in a model with subunit and forest. Something weird with GUIS-FL b/c from 2013-2017, alternate panels and the understory sd values are very different for panel 1 (GUISFL16 & 27 lowest) vs panel 2 (GUISFL 28 & 29 very high)
# > JELA: Only has one researcher

### GENERAL STEPS ----
# 0. SUBSET the data as decided, e.g., drop additional (beyond 2 per year) visits, limit to certain observers, drop certain years
# 1. EDA plots (incl. party tree and single cov. plots and interaction plots) to identify likely covariates to include in model
#   - for parks with subunits, if it looksl like covariates have very different effects by subunit, then analyze subunits separately
# 2. Run model with all combinations of the likely covariates
#   - Start by identifying best distribution--Poisson for starts, but use NB or Conway or zero-inflated...
#   - If model doesn't fit, try adjusting control values
#   - Use 'performance' to determine if best-fit is poisson, zero-inflation, NB or other distributions, etc.
#   - Using best-fit distribution, run models with all combinations of the likely covariates (incl. interactions)
#   - Use 'performance' to determine best-fit
#   - Compare coefficient estimates and CI's for top models to see if there are models with pretty different estimates, massive CI's, etc.
# 3. Run all sorts of model diagnostics (packages 'performance' and 'dharma')
# 4. Based on the model diagnostics see if I should try more models with additional covariates, quadratic forms, or with covariates on dispersion or ZI, etc. Repeat Steps 2&3 as many times as necessary to get an acceptable model or to determine that there is not a good-fit model

### FINAL ANALYSES NOTES ----
## VICK
# > My list of top 20 species is similar to what GULN said, but my total number of detections is pretty different even when I calculate straight from the data they sent me--check to see if I'm doing something funny
# > VF42 was visited 3 times in 2023. Removed the 2023_3 visit so would not screw up N-mixture

rm(list = ls())

### Load packages ----

pkgs <- c("tidyverse", "partykit", "magrittr", "DHARMa", "emmeans", "ggeffects", "tidybayes", "broom", "broom.mixed", "brms", "glmmTMB", "ggstats", "performance", "AICcmodavg", "here", "modelsummary", "cv")
installed_pkgs <- pkgs %in% installed.packages()
if (length(pkgs[!installed_pkgs]) > 0) install.packages(pkgs[!installed_pkgs], repos = "https://cloud.r-project.org" , dep=TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

### START ANALYSIS ----
source(here::here("scripts", "ANALYSIS_STARTUP.R")) # Read in formatted data
park = "VICK"
# num_species <- 20
# 
# # FIND MOST ABUNDANT SPECIES ----
# # These calculations are WITHOUT dropping detections beyond 100m and some other filters
# top_species <- df_full_obs %>%
#   dplyr::filter(unit_code == park, landbird == "yes", park_breeding == "yes") %>%
#   dplyr::group_by(species_code, common_name) %>%
#   dplyr::summarize(
#     overall_total = sum(sum_indiv, na.rm = TRUE),
#     num_per_pointsurv = mean(sum_indiv, na.rm = TRUE)) %>%
#   ungroup() %>%
#   dplyr::slice_max(n = num_species, order_by = num_per_pointsurv)
# # View(top_species)
# top_species$common_name
# 
# ### ANALYSIS FOR SELECTED SPECIES ----
# # SUBSET AND SCALE COVARIATES ----
# top_species$species_code

spec = "BLJA" # go through each of the top species...

dat <- FuncFormatFullObs(park = park, spec = spec, limit_100m = TRUE) 
# dat %<>% dplyr::filter(subunit == "GUIS-FL"); dat$subunit <- as.character(dat$subunit)
# dat %<>% dplyr::filter(subunit == "MI"); dat$subunit <- as.character(dat$subunit)

table(dat$yr_visit) # make sure nothing funny
table(dat$location_name, dat$yr)

mean(dat$sum_indiv, na.rm = TRUE)
range(dat$sum_indiv, na.rm = TRUE)
sum(dat$sum_indiv)
sum(is.na(dat$sum_indiv))

table(dat$researcher, dat$yr, dat$subunit) # check how many observers
summary(dat)

# Dataset with complete cases only
# PAAL: Dropped 3 records b/c missing weather_wind_comb
# SAAN: Dropped 1 record b/c missing hrs_since_rise
nrow(dat)
colSums(is.na(dat))

dat <- dat[complete.cases(dat),] # <<<<<<<<<< SHOULD DO THIS AFTER DETERMINING COVARIATES?
nrow(dat)

# SCALE PREDICTORS----
dat$prop_understory_cover_50_sc <- scale(dat$prop_understory_cover_50) %>% as.vector()
dat$prop_understory_cover_100_sc <- scale(dat$prop_understory_cover_100) %>% as.vector()
dat$prop_understory_cover_200_sc <- scale(dat$prop_understory_cover_200) %>% as.vector()
dat$understory_cover_sd_50_sc <- scale(dat$understory_cover_sd_50) %>% as.vector()
dat$understory_cover_sd_100_sc <- scale(dat$understory_cover_sd_100) %>% as.vector()
dat$understory_cover_sd_200_sc <- scale(dat$understory_cover_sd_200) %>% as.vector()
dat$perc_forest_50_sc <- scale(dat$perc_forest_50) %>% as.vector()
dat$perc_forest_100_sc <- scale(dat$perc_forest_100) %>% as.vector()
dat$perc_forest_200_sc <- scale(dat$perc_forest_200) %>% as.vector()
dat$perc_opendev_50_sc <- scale(dat$perc_opendev_50) %>% as.vector()
dat$perc_opendev_100_sc <- scale(dat$perc_opendev_100) %>% as.vector()
dat$perc_opendev_200_sc <- scale(dat$perc_opendev_200) %>% as.vector()

dat %>% dplyr::group_by(subunit) %>% summarize(tot = sum(sum_indiv))


### EDA ----
# PARTYKIT: Classification tree. Think of it as piecewise linear. Gives us different information from linear regression--may be useful for identifying interactions and polynomial predictors

## SAAN & GUIS--add subunit

# dat %<>% dplyr::filter(subunit=="MI") # <<<<< SAAN ONLY, AS NEEDED

## as needed... researcher + first_yr (BITH, PAAL, GUIS?)
names(dat)

dat_tree <- partykit::ctree(sum_indiv ~ subunit, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~ researcher, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~ hab_type_200 +hab_type_100+hab_type_50, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~ hab_type_200_comb +hab_type_100_comb + hab_type_50_comb, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~ hab_type_200 +hab_type_100 + hab_type_50 + hab_type_200_comb +hab_type_100_comb + hab_type_50_comb, data = dat) # <<<<<<<<<<<<<<<
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~ perc_opendev_50 + perc_opendev_100 + perc_opendev_200, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~ perc_forest_50 + perc_forest_100 + perc_forest_200, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~ hab_type_200_comb +hab_type_100_comb + hab_type_50_comb + perc_forest_50 + perc_forest_100 + perc_forest_200 + perc_opendev_50 + perc_opendev_100 + perc_opendev_200, data = dat) # <<<<<<<<<
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~prop_understory_cover_50 + prop_understory_cover_100 + prop_understory_cover_200, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~understory_cover_sd_50 + understory_cover_sd_100+understory_cover_sd_200, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~prop_understory_cover_50 + prop_understory_cover_100 + prop_understory_cover_200 + understory_cover_sd_50 + understory_cover_sd_100 + understory_cover_sd_200, data = dat) # <<<<<<<<
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~weather_wind_comb + weather_temperature_cs + julian_prop + hrs_since_rise, data = dat) # <<<<<<<<<<<<
plot(dat_tree)

# pull the main ones together now...
dat_tree <- partykit::ctree(sum_indiv ~prop_understory_cover_50 + hrs_since_rise, data = dat)
plot(dat_tree)


### MODEL NOTES ----
## NOTES ON BITH MODEL
# > exclude habitat at 200m scale bc only one level
# > incl researcher & first_yr b/c many--researcher is confounded with year. So this one is a bit funny b/c a lot of the plots seem to suggest much higher counts in the early years, but the model is attributing those to researcher differences, not to actual high numbers... so BBS models, e.g., say NOCA is declining everywhere but BITH NOCA says nonsign. increase. Will be interesting to see NOCA trend at parks with less researcher turnover
## PAAL--exclude hab b/c only non-forest. Has two researchers.
## JELA--only one researcher

### YEAR MODELS ----
### 2-VISIT, 100M LIMIT MODELS
# Note that we are using up a lot of degrees of freedom with calculating each year as a factor level, so that reduces the additional covariates we can use.

# Below, use the covariates identified initially from EDA as being potentially important predictors. May need to add other models if model diagnostics suggest we missed an important covariate.

# Options for compois model problems: diagnose(mod); update(XXXXX,control=glmmTMBControl(optimizer=optim, parallel = 8, optArgs=list(method="BFGS"), optCtrl=list(iter.max=1e3,eval.max=1e3)))

# Start with Poisson null model
mod0pois <-glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), data = dat, family = poisson)
performance::check_overdispersion(mod0pois) # Check for overdispersion. dispersion ratio >1 is overdispersion. Note that the dispersion ratio is calculated by... sum((residuals(mod0pois, type = "pearson"))^2)/df.residual(mod0pois)
performance::check_distribution(mod0pois) # See if negative binomial might be better, or if need zero-inflation

# Try NB null model IF it seems appropriate
mod0nb <- glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), data = dat, family = nbinom2)

# Try conway (slower, but use parallel processing)
mod0comp <- glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), data = dat, family="compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))

# ....year heterogeneity
mod0comp_disp_yr <- glmmTMB(formula = sum_indiv ~ yr_c_f + (1 | location_name), data = dat, family="compois", ziformula = ~0, dispformula = ~yr_c_f, control = glmmTMBControl(optCtrl = list(trace = 1), parallel = 8))

# Zero-inflated Poisson
mod0pois_zi <-glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), data = dat, family = poisson, ziformula = ~1)

# ....ZI by location
mod0pois_zi_loc <-glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), data = dat, family = poisson, ziformula = ~location_name)

# AICc to compare
mod_list <- tibble::lst(mod0pois, mod0comp, mod0comp_disp_yr, mod0pois_zi, mod0pois_zi_loc, mod0nb)
performance::compare_performance(mod_list, rank = TRUE) # CANNOT INCLUDE MODELS THAT DROP DATA HERE, THE DATASET HAS TO BE THE SAME FOR ALL COMPARED MODELS!
# aictab(mod_list, second.ord = TRUE)
ggstats::ggcoef_compare(mod_list) 

# Assign best one to mod0 (if underdispersed, it may be comp instead of pois regardless of AICc?)
mod0 <- mod0pois

### Try partykit-suggested interactions (use centered.scaled)
# Compare models then do diagnostics and adjust as necessary
# Adjustments include replacing correlated values in interactions with quadratics
mod_partybest <- update(mod0, formula = ~. + prop_understory_cover_50_sc*hrs_since_rise_c)
summary(mod_partybest)


# AICc to compare
temp_mod_list <- c(mod_list, tibble::lst(mod_party11, mod_party12))
mod_list <- temp_mod_list

mod_list <- tibble::lst(mod0, mod_partybest, mod_party1, mod_party2, mod_party3, mod_party4, mod_party5, mod_party6, mod_party7, mod_party8, mod_party9) 
(p <- performance::compare_performance(mod_list, rank = TRUE))
plot(p)
# BIC usually penalizes # of parameters more than AIC does
# ICC is more for 
ggstats::ggcoef_compare(tibble::lst(mod_list$mod_party5, mod_list$mod_party6, mod_list$mod_party8, mod_list$mod_party10))

# Check diagnostics on best-fit party model and, if needed, try additional models below

# Add covariates

## ONE BY ONE... ----
# For each group, identify the most significant scale (if any)--do not change these, they cover all available covariates
# If only linear is signif, re-run with only linear
mod1 <- update(mod0, formula = ~. + hab_type_50_comb)
summary(mod1) # 

mod2 <- update(mod0, formula = ~. + hab_type_100_comb)
summary(mod2) # 

mod3 <- update(mod0, formula = ~. + hab_type_200_comb)
summary(mod3) 

mod4 <- update(mod0, formula = ~. + perc_forest_50_sc + I(perc_forest_50_sc^2))
summary(mod4) # 

mod5 <- update(mod0, formula = ~. + perc_forest_100_sc + I(perc_forest_100_sc^2))
summary(mod5) # <<< .52, 0.046

mod6 <- update(mod0, formula = ~. + perc_forest_200_sc + I(perc_forest_200_sc^2))
summary(mod6) 

mod7 <- update(mod0, formula = ~. + perc_opendev_50_sc + I(perc_opendev_50_sc^2))
summary(mod7) #

mod8 <- update(mod0, formula = ~. + perc_opendev_100_sc + I(perc_opendev_100_sc^2))
summary(mod8) # 

mod9 <- update(mod0, formula = ~. + perc_opendev_200_sc + I(perc_opendev_200_sc^2))
summary(mod9) 

mod10 <- update(mod0, formula = ~. + understory_cover_sd_50_sc + I(understory_cover_sd_50_sc^2))
summary(mod10) # 

mod11 <- update(mod0, formula = ~. + understory_cover_sd_100_sc + I(understory_cover_sd_100_sc^2))
summary(mod11) # 

mod12 <- update(mod0, formula = ~. + understory_cover_sd_200_sc)#+ I(understory_cover_sd_200_sc^2))
summary(mod12) # 0.026

mod13 <- update(mod0, formula = ~.+ prop_understory_cover_50_sc) # + I(prop_understory_cover_50_sc^2))
summary(mod13) # <<<<<<<<<<<<<<<

mod14 <- update(mod0, formula = ~.+ prop_understory_cover_100_sc + I(prop_understory_cover_100_sc^2))
summary(mod14) 

mod15 <- update(mod0, formula = ~.+ prop_understory_cover_200_sc + I(prop_understory_cover_200_sc^2))
summary(mod15) # 

mod16 <- update(mod0, formula = ~. + julian_prop_c + I(julian_prop_c^2))
summary(mod16) # 

mod17 <- update(mod0, formula = ~. + hrs_since_rise_c)# + I(hrs_since_rise_c^2))
summary(mod17) # 0.007

mod18 <- update(mod0, formula = ~. + weather_wind_comb)
summary(mod18) # 

mod19 <- update(mod0, formula = ~. + weather_temperature_cs)
summary(mod19) # 

# Compare models
mod_list <- tibble::lst(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15, mod16, mod17, mod18, mod19) # see how these rankings compare to EDA
performance::compare_performance(mod_list, rank = TRUE)

## FULL MODEL... ----
# Include all variables that were statistically signif (or close) in single models OR in partykit, possibly ordered by compare_performance. For each, use the scale that is "most" signif.
# Include polynomials if signif in regression tree

mod_full <- update(mod0, formula = ~. + prop_understory_cover_50_sc + I(prop_understory_cover_50_sc^2) + understory_cover_sd_200_sc + hrs_since_rise_c + I(hrs_since_rise_c^2) + perc_forest_100_sc + I(perc_forest_100_sc^2))
summary(mod_full)

# >>> NOTES:

## NOW SUBSETS REMOVING THE LEAST IMPORTANT ONE BY ONE, THEN IN GROUPS... 
# remove one at a time, check the model, remove next least important
# swap out some correlated ones also

mod_sub1 <- update(mod_full, formula = ~. -I(prop_understory_cover_50_sc^2))
summary(mod_sub1)

mod_sub2 <- update(mod_full, formula = ~.-I(prop_understory_cover_50_sc^2) - I(hrs_since_rise_c^2))
summary(mod_sub2)

mod_sub3 <- update(mod_full, formula = ~. -I(prop_understory_cover_50_sc^2) - I(hrs_since_rise_c^2) - perc_forest_100_sc - I(perc_forest_100_sc^2))
summary(mod_sub3)

mod_sub4 <- update(mod_full, formula = ~. -I(prop_understory_cover_50_sc^2) - I(hrs_since_rise_c^2) -understory_cover_sd_200_sc)
summary(mod_sub4)

mod_sub5 <- update(mod_full, formula = ~. -I(prop_understory_cover_50_sc^2) - I(hrs_since_rise_c^2) - perc_forest_100_sc - I(perc_forest_100_sc^2)-understory_cover_sd_200_sc)
summary(mod_sub5)



# Compare models
# If it's not the slimmest model that is chosen, go back and drop things in different orders from the best-fit-so-far
# Include mod0 and party models
mod_list <- tibble::lst(mod0, mod_partybest, mod_full, mod_sub1, mod_sub2, mod_sub3, mod_sub4, mod_sub5) #, mod_sub3, mod_sub4, mod_sub5, mod_sub6, mod_sub7, mod_sub8)
(p <- performance::compare_performance(mod_list, rank = TRUE))
plot(p)
ggstats::ggcoef_compare(tibble::lst(mod13, mod_sub5, mod_sub3)) # do this for subset of models

# CROSS-VALIDATE subset of models ----
#Include the BIC-best-supported. Default metric is MSE. Smaller MSE is better. Does this tend to go with the BIC-best-fit?
mod_dat <- dat
(cv_yearmod_list <- cv::cv(cv::models(mod_best_comp = mod_best_comp, mod_sub3 = mod_sub3, mod_sub5 = mod_sub5), data = mod_dat, clusterVariables = "location_name", ncores = 8, k = 20))
# SAVE THE CV LIST
saveRDS(cv_yearmod_list, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_cv_yearmod_list.RDS")))

# # Check also with LOO--this often hangs up
# (LOO_year_mod_list <- cv::cv(cv::models(mod13 = mod13, mod_sub3 = mod_sub3, mod_sub5 = mod_sub5), data = mod_dat, ncores = 8, k = "loo", start = TRUE))


## THIS IS THE MODEL TO EVALUATE...

year_mod <- mod_sub5
mod_dat <- dat 

### SAVE THE FULL MOD_LIST
# Full model list <<<<<<< BRING TOGETHER ALL MOD_LISTS
mod_list <- tibble::lst(mod0pois, mod0comp,mod0pois_zi, mod0pois_zi_loc, mod0nb, mod_partybest, mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15, mod16, mod17, mod18, mod19, mod_full, mod_sub1, mod_sub2, mod_sub3, mod_sub4, mod_sub5) 
saveRDS(mod_list, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_yearmod_modlist.RDS")))

### SAVE THE FINAL 100M 2-VISIT YEAR MODEL ----
saveRDS(year_mod, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_yearmod.RDS")))
saveRDS(mod_dat, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_dat.RDS")))


## [AS NEEDED] TRY ADJUSTED MODELS ----
# DO THIS ONLY IF DIAGNOSTICS ARE NOT GOOD!! COME BACK AND MAKE MODEL ADJUSTMENTS BASED ON DIAGNOSTICS... DROP UNUSUAL SITES, ADD INTERACTIONS, OTHER CHANGES AS NEEDED... ADD THEM TO THE FULL LIST AND REPEAT
# try dropping unusual sites
# NOTE: Cannot do model selection when data differ
mod_adj1 <- update(mod_sub3, dispfo)

##COMPARE BEST SET OF MODELS 
mod_list = tibble::lst(mod0pois, mod0comp, mod0pois_zi, mod0nb, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod_full, mod_sub1, mod_sub2, mod_sub3, mod_sub4, mod_sub5, mod_sub6, mod_adj1, mod_adj2, mod_adj3, mod_adj4, mod_adj5, mod_adj6, mod_adj7, mod_adj8)

mod_list = tibble::lst(year_mod, mod_adj1, mod_adj2, mod_adj3, mod_adj4, mod_adj5, mod_adj6, mod_adj7, mod_adj8)
(p <- performance::compare_performance(mod_list, rank = TRUE))
plot(p) # Check which models perform best with R^2 vs BIC vs AIC

ggstats::ggcoef_compare(tibble::lst(mod_full, mod_sub3, mod_sub5)) # do this for subset of models

## [AS NEEDED] Go back and check different model distributions on best model ----
year_mod$call

mod_best_pois <-glmmTMB(formula = sum_indiv ~ yr_c_f + (1 | location_name) + 
                          understory_cover_sd_50_sc + I(understory_cover_sd_50_sc^2) + 
                          prop_understory_cover_100_sc + I(prop_understory_cover_100_sc^2) + 
                          hrs_since_rise_c, data = dat, family = poisson)

mod_best_zi <-glmmTMB(formula = sum_indiv ~ yr_c_f + (1 | location_name) + 
                        understory_cover_sd_50_sc + I(understory_cover_sd_50_sc^2) + 
                        prop_understory_cover_100_sc + I(prop_understory_cover_100_sc^2) + 
                        hrs_since_rise_c, data = dat, family = poisson, ziformula = ~1)

mod_best_nb <- glmmTMB(formula = sum_indiv ~ yr_c_f + (1 | location_name) + 
                         understory_cover_sd_50_sc + I(understory_cover_sd_50_sc^2) + 
                         prop_understory_cover_100_sc + I(prop_understory_cover_100_sc^2) + 
                         hrs_since_rise_c, data = dat, family = nbinom2)

mod_best_comp <- glmmTMB(formula = sum_indiv ~ yr_c_f + (1 | location_name) + 
                           prop_understory_cover_50_sc + hrs_since_rise_c, data = dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))

# Comparison of different distributions with best-fit model... redo submodels with 'best' distribution if it changes here
mod_list_distrib = tibble::lst(year_mod, mod_best_comp) #mod_best_pois, mod_best_zi, mod_best_nb)
performance::compare_performance(mod_list_distrib, rank = TRUE)

## [AS NEEDED] DROP SITES OR OTHER ADJUSTMENTS FROM BEST-FIT MODEL ----
# Do this if DHARMA results suggest one or a few sites are driving odd diagnostics
# mod_drop3 <- update(year_mod, data = dat %>% dplyr::filter(!location_name %in% c("VF36", "VM07")) %>% droplevels(.)) # mod_sub4, but drop site with the highest single survey count (N=4)
ggstats::ggcoef_compare(tibble::lst(year_mod, mod_best_comp))

new_mod_list <- c(mod_list, tibble::lst(mod_dropVF31))
saveRDS(new_mod_list, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_yearmod_modlist.RDS")))


### DIAGNOSTICS ON BEST-FIT MODEL ----
mod <- year_mod
# mod_dat <- dat %>% dplyr::filter(!location_name %in% c("VF36", "VM07")) #dat # dat %>% dplyr::filter(location_name != "VF31") %>% droplevels(.)

# mod_list <- readRDS(here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_yearmod_modlist.RDS")))

summary(mod) # make sure no absurd standard errors
ggstats::ggcoef_table(mod, intercept = FALSE, table_witdhs = c(2, 1)) + labs(title = "XXX")

## pp-checks
performance::check_model(mod)
performance::check_predictions(mod, check_range = TRUE)
performance::check_heterogeneity_bias(mod)
performance::check_collinearity(mod) # If interaction terms are included in a model, high VIF values are expected. 
performance::check_autocorrelation(mod) 
performance::check_overdispersion(mod)
performance::check_zeroinflation(mod)

## CHECK FOR HIGH LEVERAGE LOCATIONS ----
## >>> CANNOT YET CHECK FOR INFLUENTIAL POINTS IN GLMMTMB MODELS, BUT THEY ARE WORKING ON IT... https://github.com/glmmTMB/glmmTMB/pull/1065
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
# For influential POINTS...
mod_leverage_case <- influence_mixed(mod, groups=".case")

# For influential locations...
mod_leverage <- influence_mixed(mod, groups="location_name")

# Convert to data frame, then heat map
x <- as.data.frame(mod_leverage$`fixed.effects[-location_name]`)
x_long <- x %>% rownames_to_column(var = "location_name") %>% tidyr::pivot_longer(data = ., cols = -location_name, names_to = "FE", values_to = "coef")

x_long_sub <- x_long %>% dplyr::filter(!FE %like% "yr_c_f")

p <- ggplot(x_long_sub, aes(x = FE, y = location_name, fill = coef)) + geom_tile() + geom_text(aes(label = round(coef, 2)), size = 4, color = "white") # this heat map shows for each FE, what the estimated coefficient is if you exclude that location from the an
ggplotly(p)

### Dharma checks for glmmTMB ----
simres <- DHARMa::simulateResiduals(fittedModel = mod, n = 1000, refit = FALSE, plot = TRUE) # re-simulate all levels, to test model structure as a whole. The default is resimulating on all levels (unconditional sims) instead of simulating conditional on the fitted random effects.
simres <- DHARMa::simulateResiduals(fittedModel = mod, n = 1000, refit = FALSE, re.form = NULL, plot = TRUE) #re-simulate conditional on the fitted random effects, i.e., the RE's are not re-simulated. More sensitive for detecting problems, but it doesn't test all possible options.

DHARMa::testZeroInflation(simres)

DHARMa::testDispersion(simres) # test for over/under dispersion
DHARMa::testOutliers(simres, type = "bootstrap")



# Plot residuals vs each independent factor ----
# Identify the ones with problems--were any of these identified as important earlier?
# Is weirdness due to just one or two odd sites?
DHARMa::plotResiduals(simres, mod_dat$yr_c)
DHARMa::plotResiduals(simres, as.factor(mod_dat$yr_c_f))
DHARMa::plotResiduals(simres, mod_dat$hab_type_50_comb, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$hab_type_100_comb, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$hab_type_200_comb, rank  = FALSE)

DHARMa::plotResiduals(simres, mod_dat$perc_forest_50, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$perc_forest_100, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$perc_forest_200, rank  = FALSE) # <

DHARMa::plotResiduals(simres, mod_dat$perc_opendev_50, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$perc_opendev_100, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$perc_opendev_200, rank  = FALSE)

DHARMa::plotResiduals(simres, mod_dat$prop_understory_cover_50, rank  = FALSE) #
DHARMa::plotResiduals(simres, mod_dat$prop_understory_cover_100, rank  = FALSE) # 
DHARMa::plotResiduals(simres, mod_dat$prop_understory_cover_200, rank  = FALSE) # <<


DHARMa::plotResiduals(simres, mod_dat$understory_cover_sd_50, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$understory_cover_sd_100, rank  = FALSE) # <<
DHARMa::plotResiduals(simres, mod_dat$understory_cover_sd_200, rank  = FALSE)

DHARMa::plotResiduals(simres, mod_dat$weather_wind_comb, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$julian_prop_c, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$hrs_since_rise_c, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$researcher, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$first_yr, rank  = FALSE)

# ### Append Dharma residuals to model data ----
dat_append <- cbind(mod_dat, predicted = DHARMa::getFitted(mod), resid = DHARMa::getResiduals(mod))
dat_append$scaled_resid = scales::rescale(dat_append$resid, to = c(0,1), from = range(dat_append$resid, na.rm = TRUE))
dat_append$rank_pred = datawizard::ranktransform(dat_append$understory_cover_sd_100)
plotly::ggplotly(ggplot(dat_append, aes(x=understory_cover_sd_100, y = scaled_resid, group = 1, text = paste0(location_name, "<br>Observed: ", sum_indiv, "<br>PredictedResponse: ", round(predicted, 2), "<br>Cov: ", understory_cover_sd_100, "<br>Resid: ", round(resid,2)))) + geom_point() + geom_quantile(), tooltip = "text")

### SAVE THE FINAL 100M 2-VISIT YEAR MODEL ----
saveRDS(year_mod, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_yearmod.RDS")))
saveRDS(mod_dat, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_dat.RDS")))

### CHECK THE PREDICTION PLOT ----
# NOTE THESE ARE WALD CI'S!!!

## ggemmeans to get year predictions <<<<<<<<<<<<<<<<<<<<<<< CHANGE THIS IN FINAL SUMMARIES
# ggeffects::predict_response(mod, terms = c("yr_c_f", "subunit"), margin = "mean_reference", type = "zero_inflated")) # > Note that to get emmeans, 'margin = "marginalmeans"' but then can't do zero-inflation with that
# (temp_tmb_pred_yearmod <- ggemmeans(mod_tmb_yearmod, terms=c("yr_c_f"), type = "fixed") %>% # # type = ifelse(tmb_zi==TRUE, "fe.zi", "fixed"))) #### ZERO-INFLATION PART DOESN' WORK!!!

tmb_zi <- ifelse(mod$modelInfo$allForm$ziformula=="~1", TRUE, FALSE)

tmb_pred_yearmod <- ggpredict(mod, terms=c("yr_c_f"), type = ifelse(tmb_zi==TRUE, "zero_inflated", "fixed")) # This is same as below

(tmb_pred_yearmod <- ggeffects::predict_response(mod, terms=c("yr_c_f"), margin = "mean_mode", type = ifelse(tmb_zi==TRUE, "zero_inflated", "fixed")))

(tmb_pred_yearmod %<>% 
   as.data.frame() %>%
   janitor::clean_names() %>%
    dplyr::mutate(yr_c = as.numeric(as.character(x))) %>%
    dplyr::select(
      group,
      yr_c,
      yr_predicted = predicted,
      yr_conf_low = conf_low,
      yr_conf_high = conf_high))
tmb_pred_yearmod$yr <- tmb_pred_yearmod$yr_c+(range(dat$yr) %>% median())

ggplot(data = tmb_pred_yearmod) + 
  geom_errorbar(
    mapping = aes(x = yr, y = yr_predicted, ymin = yr_conf_low, ymax = yr_conf_high, color = group), width = 0) +
  geom_point(
    mapping = aes(x = yr, y = yr_predicted, color = group), size = 2.5) +
  facet_grid(rows = vars(group))


### TRY BRMS ----
# Figure out how to do best weakly informative priors
# If it doesn't work, try on personal laptop with cmdstanr
mod_brms <- brms::brm(
  formula = sum_indiv ~ yr_c_f + (1 | location_name) + 
    prop_understory_cover_50_sc + hrs_since_rise_c,
  data = mod_dat,
  family = brmsfamily("poisson"),
  # family = brmsfamily("com_poisson"),
  # sample_prior = "only", # <<<<<<<<<<< CHECK PRIORS
  iter = 4000, warmup = 1000, chains = 4, cores = 4,
  file = here::here("brms_mods", paste0("brms_", park, "_", spec, "_100m_yearmod.RDS"))
)
pp_check(mod_brms, ndraws = 200)
mod_brms <- add_criterion(mod_brms, criterion = c("loo", "waic"), force_save = TRUE)
loo(mod_brms)
saveRDS(mod_dat, here::here("brms_mods", paste0("brms_", park, "_", spec, "_100m_yearmod_dat.RDS")))


### TREND MODELS ----
# For random slopes...+ (1 + yr_c | location_name)). You can actually omit the '1+' b/c it is assumed--to manually remove intercept RE, would do '0+'. In all models, include random intercept AND slope. For some models, it may also make sense to have a fixed effect interaction of habitat type with year. The base GLMM model is therefore:  Overall trend (yr_c) plus year to year fluctuations concordant across locations (1|yr_f) plus random variation among locations in location-specific trends ( 1 + yr_c | location_name). Also, note   https://academic.oup.com/esr/article-abstract/35/2/258/5306121?redirectedFrom=fulltext
# >>>>> CHECK IF I SHOULD HAVE ADDED DISPFORMULA =~1, OR IF THAT IS DEFAULT
# If removed any locations for best-fit mod, use that same data set for trend model

# Start with the best-fit year model, then try models with hab*yr interaction if that seems appropriate and use diagnostics, etc. to see if should drop/add other covariates
# Try dropping covariates that in year model were only driven by one or two sites


# mod_dat <- dat %>% dplyr::filter(!location_name %in% c("VO08")) %>% droplevels()
# 
# p <- ggplot(dat, aes(x = yr, y = sum_indiv, color = hab_type_100, fill = hab_type_100, group_by = location_name)) +
#          stat_smooth(method = "lm", se = TRUE)
# plotly::ggplotly(p)
# 
# raw_trend <- lm(formula = sum_indiv ~ yr_c*location_name, data = mod_dat)
# trend_coef <- broom::tidy(raw_trend)
# trend_df <- trend_coef[grep("yr_c:location_name", trend_coef$term), ]
# trend_df$term <- gsub(pattern = "yr_c:location_name", replacement = "", trend_df$term) # then color by estimate



mod$call
### BASE TREND MODEL ----
# BASE model is trend with year trend, year RE and location intercept&slope RE 
trend_mod_base <- glmmTMB(sum_indiv ~ yr_c + (1 | yr_c_f) + (1 + yr_c | location_name), data = mod_dat, family = "poisson")
                        # "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8)) 
summary(trend_mod_base)

# IF RE variance is practically zero, consider dropping it
trend_mod_base2 <- glmmTMB(sum_indiv ~ yr_c + (1 + yr_c | location_name), data = mod_dat, family = "poisson")
summary(trend_mod_base2)

# IF the RE's are highly correlated, drop the slope
trend_mod_base3 <- glmmTMB(sum_indiv ~ yr_c + (1|location_name), data = mod_dat, family = "poisson")
summary(trend_mod_base3)

# Compare base models
# Exponentiate to get % of population in next year. exp(trend_coef)^(1:10) to see % over next 10 years. Need to bootstrap CI to get correct CI on trend.
trend_mod_base_list <- tibble::lst(trend_mod_base, trend_mod_base2, trend_mod_base3)
performance::compare_performance(trend_mod_base_list, rank = TRUE)
ggstats::ggcoef_compare(trend_mod_base_list) 
(cv_trend_base_mod <- cv::cv(models(trend_mod0=trend_mod0, trend_mod0_Xhab=trend_mod0_Xhab), data = mod_dat, reps = 5))
 
trend_mod_base_RE <- trend_mod_base3 # <<<<<<<<<<<<< ENTER HERE

## HABITAT INTERACTION WITH YEAR-SLOPE----
trend_mod_base_hab50 <- update(trend_mod_base_RE, formula = ~.+ yr_c*hab_type_50_comb)
summary(trend_mod_base_hab50)

trend_mod_base_hab100 <- update(trend_mod_base_RE, formula = ~.+ yr_c*hab_type_100_comb)
summary(trend_mod_base_hab100) # p = 0.043

trend_mod_base_hab200 <- update(trend_mod_base_RE, formula = ~.+ yr_c*hab_type_200_comb)
summary(trend_mod_base_hab200)


# Check other distributions
trend_mod0_compois <- update(trend_mod0, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))
summary(trend_mod0_compois)

performance::compare_performance(tibble::lst(trend_mod0, trend_mod0_compois), rank = TRUE)

trend_mod <- trend_mod3 # ASSIGN THE BEST-FIT MODEL

## [AS NEEDED] ADJUST MODEL AFTER CHECKING DIAGNOSTICS ----
trend_mod1 <- update(trend_mod, formula = ~. + prop_understory_cover_50_sc)
summary(trend_mod1)

trend_mod2 <- update(trend_mod, formula = ~. + hrs_since_rise_c)
summary(trend_mod2)

trend_mod3 <- update(trend_mod, formula = ~. +prop_understory_cover_50_sc + hrs_since_rise_c)
summary(trend_mod3)

trend_mod4 <- update(trend_mod, formula = ~. +prop_understory_cover_50_sc + hrs_since_rise_c + understory_cover_sd_200_sc)
summary(trend_mod4)

test <- update(trend_mod, formula = sum_indiv ~ yr_c*hab_type_100_comb + (1 | location_name) + prop_understory_cover_50_sc + hrs_since_rise_c)
summary(test)


trend_mod_list <- tibble::lst(trend_mod_base, trend_mod_base2, trend_mod_base3, trend_mod1, trend_mod2, trend_mod3, trend_mod4)
performance::compare_performance(trend_mod_list, rank = TRUE)
ggstats::ggcoef_compare(tibble::lst(trend_mod1, trend_mod3, trend_mod4)) 

### EXPLORE RANDOM EFFECTS ----
re <-ranef(trend_mod0)[1]$cond$location_name %>%
  janitor::clean_names()
re
ggplot(re, aes(x = intercept, y = yr_c)) + geom_point()

### ASSIGN TREND_MOD
trend_mod <- trend_mod3
mod <- trend_mod

### SAVE TRENDMOD RESULTS ----
saveRDS(trend_mod_list, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_trendmod_modlist.RDS")))
saveRDS(trend_mod, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_trendmod.RDS")))
saveRDS(mod_dat, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_trendmod_dat.RDS")))

### TRY BRMS ----
# Figure out how to do best weakly informative priors
# If it doesn't work, try on personal laptop with cmdstanr
trendmod_brms <- brms::brm(
  formula = sum_indiv ~ yr_c * hab_type_100_comb + (1 |yr_c_f) + (1 | location_name) + prop_understory_cover_50_sc + I(prop_understory_cover_50_sc^2) + understory_cover_sd_100_sc + 
    I(understory_cover_sd_100_sc^2) + hrs_since_rise_c + I(hrs_since_rise_c^2) + 
    perc_opendev_200_sc,
  data = mod_dat,
  family = brmsfamily("com_poisson"),
  # sample_prior = "only", # <<<<<<<<<<< CHECK PRIORS
  iter = 4000, warmup = 1000, chains = 4, cores = 4,
  file = here::here("brms_mods", paste0("brms_", park, "_", spec, "_100m_trendmod.RDS"))
)
pp_check(trendmod_brms, ndraws = 500)



### PLOT glmmTMB PREDICTIONS OF FINAL TRENDS MODEL ----
# Check the estimated trend
em=emmeans::emtrends(trend_mod, specs = c("yr_c", "hab_type_100_comb"), var = "yr_c") # trends at mean year. Ignores a 'type' argument, it seems
summary(em, infer=c(TRUE,TRUE),null=0) %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  dplyr::mutate(signif = p_value < 0.05)

(temp_tmb_pred <- ggemmeans(trend_mod, terms=c("yr_c", "hab_type_100_comb"), type = "fixed"))


### PLOT TREND WITH YEAR ESTIMATES ----
tmb_zi_trendmod <- ifelse(trend_mod$modelInfo$allForm$ziformula=="~1", TRUE, FALSE)

(temp_trendmod_mmeans <- ggeffects::predict_response(trend_mod, terms=c("yr_c", "hab_type_100_comb"), margin = "marginalmeans", type = ifelse(tmb_zi_trendmod==TRUE, "zero_inflated", "fixed"))) # this is the same as emmeans. For factors, it weights them equally regardless of how many are actually in each factor level
    
(temp_trendmod_mmode <- ggeffects::predict_response(trend_mod, terms=c("yr_c", "hab_type_100_comb"), margin = "mean_mode", type = ifelse(tmb_zi_trendmod==TRUE, "zero_inflated", "fixed")))    

temp_trendmod_mmeans %<>% 
  as.data.frame 

FuncFormatTrendDat <- function(pred, dat) {
  pred_out <- pred %>% 
    as.data.frame()
  
  pred_out$yr <- pred_out$x+(range(dat$yr) %>% median())
  

  pred_out %<>%  
    janitor::clean_names() %>%
    dplyr::mutate(yr_c = as.numeric(x))
  
  if("group" %in% names(pred_out)) {
    
    pred_out %<>%
      dplyr::select(
        group,
        yr,
        yr_c,
        trend_predicted = predicted,
        trend_conf_low = conf_low,
        trend_conf_high = conf_high)
  } else {
    pred_out %<>%
      dplyr::select(
        yr,
        yr_c,
        trend_predicted = predicted,
        trend_conf_low = conf_low,
        trend_conf_high = conf_high)
  }
  
  return(pred_out)
}

tmb_pred_trendmod_mmeans <- FuncFormatTrendDat(pred = temp_trendmod_mmeans, dat)
# tmb_pred_trendmod_mmode <- FuncFormatTrendDat(pred = temp_trendmod_mmode, dat) # <<<<<<<<<<<<< FIX--THIS ISN'T QUITE RIGHT AS A COMPARISON TO YEARMOD

(p <- ggplot() + 
  geom_errorbar(data = tmb_pred_yearmod, aes(x = yr, y = yr_predicted, ymin = yr_conf_low, ymax = yr_conf_high), width = 0, color = "blue") +
  geom_point(data = tmb_pred_yearmod, aes(x = yr, y = yr_predicted), size = 2.5, color = "blue") +
  geom_ribbon( # glmmTMB TREND RESULTS AS RIBBON
    data = tmb_pred_trendmod_mmeans, aes(x = yr, y = trend_predicted, ymin = trend_conf_low, ymax = trend_conf_high),
    alpha = 0.1, fill = "orange") +
  geom_line(data =tmb_pred_trendmod_mmeans, aes(x = yr, y= trend_predicted), color = "orange", linetype = "dashed") +
  # geom_ribbon(
  #   data = tmb_pred_trendmod_mmode, aes(x = yr, y = trend_predicted, ymin = trend_conf_low, ymax = trend_conf_high),
  #   alpha = 0.1, fill = "red") +
  # geom_line(data =tmb_pred_trendmod_mmode, aes(x = yr, y= trend_predicted), color = "red", linetype = "dashed") +
  scale_x_continuous(breaks = unique(sort(dat$yr))) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title = paste0("Trends in relative abundance of ", spec, " at ", park), x = "Year", y = "Index of abundance (# of birds/point count)") +
  theme_bw() +
  facet_grid(rows = vars(group)))

saveRDS(p, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_trendplot.RDS")))
