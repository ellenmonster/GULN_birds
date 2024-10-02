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

## PAAL
# > Drop PAAL21 due to no habitat info

rm(list = ls())

### Load packages ----

pkgs <- c("tidyverse", "partykit", "magrittr", "DHARMa", "emmeans", "ggeffects", "tidybayes", "broom", "broom.mixed", "brms", "glmmTMB", "ggstats", "performance", "AICcmodavg", "here", "modelsummary", "cv")
installed_pkgs <- pkgs %in% installed.packages()
if (length(pkgs[!installed_pkgs]) > 0) install.packages(pkgs[!installed_pkgs], repos = "https://cloud.r-project.org" , dep=TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# options(mc.cores = 4, brms.backend = "cmdstanr")
### START ANALYSIS ----
source(here::here("scripts", "ANALYSIS_STARTUP.R")) # Read in formatted data
park = "PAAL"
spec = "EAME" # go through each of the top species...

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

dat <- FuncFormatFullObs(park = park, spec = spec, limit_100m = TRUE) # <<<<<<<<<<<<< DROP_SITES OR NOT
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

### TEMPORARY UNTIL CORRECTED HABITAT DATA AVAILABLE----
# DROP PAAL21 B/C NO HABITAT INFO (IT'S OUTSIDE OF PARK)
dat %<>% dplyr::filter(location_name != "PAAL21") # <<<<<<<<<<<<<<<<<<<<<<<

### Check data missingness----
# PAAL: 3 records missing weather_wind_comb; 2 records missing weather_temperature
# SAAN: 1 record missing hrs_since_rise
colSums(is.na(dat)) %>% sort() # which covariates have some missing data? Keep this in mind -- if it's an important covariate in the final model, should drop the records and rerun the model selection
sum(dat$sum_indiv) # total detections in final dataset
full_dat <- dat # hold this aside in case need to put records back b/c missingness wont' affect final model

nrow(dat)
dat <- dat[complete.cases(dat),] 
colSums(is.na(dat))
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

# >>>>> Check for complete separation
dat %>% dplyr::group_by(weather_wind_comb) %>% summarize(mean(sum_indiv))
dat %>% dplyr::group_by(yr) %>% summarize(mean(sum_indiv))

# # >>>>> Check for covariate patterns that may explain unusually low or high years
# table(dat$yr_c, dat$weather_wind_comb)

### EDA ----

## >>>>> FIRST, LOOK ON DASHBOARD FOR YEAR TRENDS, ODD YEARS (OR YRS WITH NO DETECTIONS), UNUSUAL HIGH NUMBERS (E.G., LIKELY INFLUENTIAL POINTS FOR TREND MODEL), NUMBER OF LOCATIONS WITH NO DETECTIONS EVER, PREDICTOR RELATIONSHIPS. ALSO LOOK UP HABITAT FOR THE SPECIES TO SET EXPECTATIONS ABOUT PREDICTORS. LOOK FOR EVIDENCE OF COMPLETE SEPARATION.

# >>>> CHECK WHICH BIRDS IF ANY SHOW DIFFERENT PATTERN FOR HABITAT COVER @ DIFFERENT SCALES

# PARTYKIT: Classification tree. Think of it as piecewise linear. Gives us different information from linear regression--may be useful for identifying interactions and polynomial predictors. Sometimes it identifies a fork, but the numbers don't differ much on either side of the fork or the sample size is really small for one fork or the cutoff for the fork is near a boundary. Consider those with a grain of salt.

## SAAN & GUIS--add subunit

# dat %<>% dplyr::filter(subunit=="MI") # <<<<< SAAN ONLY, AS NEEDED

## as needed... researcher + first_yr (BITH, PAAL, GUIS?)

dat_tree <- partykit::ctree(sum_indiv ~prop_understory_cover_50_sc + prop_understory_cover_100_sc +  perc_forest_50_sc + perc_forest_100_sc + perc_opendev_50_sc + perc_opendev_100_sc +  understory_cover_sd_100_sc + weather_wind_comb + weather_temperature_cs + julian_prop_sc + hrs_since_rise_sc , data = dat) # 
plot(dat_tree)

### MODEL NOTES ----
## NOTES ON BITH MODEL
# > exclude habitat at 200m scale bc only one level
# > incl researcher & first_yr b/c many--researcher is confounded with year. So this one is a bit funny b/c a lot of the plots seem to suggest much higher counts in the early years, but the model is attributing those to researcher differences, not to actual high numbers... so BBS models, e.g., say NOCA is declining everywhere but BITH NOCA says nonsign. increase. Will be interesting to see NOCA trend at parks with less researcher turnover
## PAAL--exclude hab_type b/c only non-forest. Has multiple researchers.
## JELA--only one researcher

### YEAR MODELS ----
### 2-VISIT, 100M LIMIT MODELS
# Note that we are using up a lot of degrees of freedom with calculating each year as a factor level, so that reduces the additional covariates we can use.

# Below, use the covariates identified initially from EDA as being potentially important predictors. May need to add other models if model diagnostics suggest we missed an important covariate.

# Options for compois model problems: diagnose(mod); update(XXXXX,control=glmmTMBControl(optimizer=optim, parallel = 8, optArgs=list(method="BFGS"), optCtrl=list(iter.max=1e3,eval.max=1e3)))

# Start with Poisson null model
mod0pois <-glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), data = dat, family = poisson)
summary(mod0pois)
performance::check_overdispersion(mod0pois, alternative = "two.sided") # Check for overdispersion. dispersion ratio >1 is overdispersion. Note that the dispersion ratio is calculated by... sum((residuals(mod0pois, type = "pearson"))^2)/df.residual(mod0pois). Statistical signif means to reject equidispersion (i.e., indicates it's overdispersed or underdispersed).
# performance::check_distribution(mod0pois) # See if negative binomial might be better, or if need zero-inflation

# # # Try NB null model IF it seems appropriate
# mod0nb <- glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), data = dat, family = nbinom2)
# summary(mod0nb)

# Try conway (slower, but use parallel processing)
mod0comp <- glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), data = dat, family="compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))
summary(mod0comp)

# Zero-inflated Poisson
mod0pois_zi <-glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), data = dat, family = poisson, ziformula = ~1)
# summary(mod0pois_zi)
# plogis(-21.62)

# AICc to compare
mod_list <- tibble::lst(mod0pois, mod0comp, mod0pois_zi)# mod0nb) #, mod0pois_zi_loc)#, mod0nb)
performance::compare_performance(mod_list, rank = TRUE, metrics = "common") # CANNOT INCLUDE MODELS THAT DROP DATA HERE, THE DATASET HAS TO BE THE SAME FOR ALL COMPARED MODELS!

# Assign best one to mod0 (prioritize BIC)

mod0 <- mod0pois


## ONE BY ONE... ----
# TRy quadratic if EDA suggests that is more appropriate
# >>> CHECK FOR VERY STRANGE ESTIMATES!
# For each group, identify the most significant scale (if any)--do not change these, they cover all available covariates
# TRy quadratic if EDA suggests that is more appropriate

mod1 <- update(mod0, formula = ~. + perc_forest_100_sc)# + I(perc_forest_100_sc^2))
summary(mod1)  #

mod2 <- update(mod0, formula = ~. + perc_opendev_100_sc)# + I(perc_opendev_200_sc^2))
summary(mod2) # 

mod3 <- update(mod0, formula = ~. + understory_cover_sd_100_sc)#+ I(understory_cover_sd_100_sc^2))
summary(mod3) # 

mod4 <- update(mod0, formula = ~.+ prop_understory_cover_100_sc ) #+ I(prop_understory_cover_100_sc^2))
summary(mod4) # 

mod5 <- update(mod0, formula = ~. + julian_prop_sc)# + I(julian_prop_sc^2))
summary(mod5) # 

mod6 <- update(mod0, formula = ~. + hrs_since_rise_sc)# + I(hrs_since_rise_sc^2))
summary(mod6) # 

mod7 <- update(mod0, formula = ~. + weather_wind_comb)
summary(mod7) #

mod8 <- update(mod0, formula = ~. + weather_temperature_cs)
summary(mod8) #

mod9 <- update(mod0, formula = ~. + prop_understory_cover_50_sc)
summary(mod9) # 

mod10 <- update(mod0, formula = ~. + perc_forest_50_sc)# + I(perc_forest_50_sc^2))
summary(mod10)  #

mod11 <- update(mod0, formula = ~. + perc_opendev_50_sc)# + I(perc_opendev_200_sc^2))
summary(mod11) # 

# Compare models
mod_list <- tibble::lst(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11) 
(perf <- performance::compare_performance(mod_list, rank = TRUE, metrics = "common")) # How much variability does the model actually explain?


## FULL MODEL... ----
# Include all variables that were statistically signif (or close) in single models OR in partykit, possibly ordered by compare_performance. For each, use the scale that is "most" signif.
# Include polynomials if signif in regression tree

mod_full <- update(mod0, formula = ~. + prop_understory_cover_100_sc + weather_wind_comb + understory_cover_sd_100_sc + hrs_since_rise_sc)
summary(mod_full)

## NOW SUBSETS REMOVING THE LEAST IMPORTANT ONE BY ONE, THEN IN GROUPS - IF THERE ARE >2 COVARIATES... 
# remove one at a time, check the model, remove next least important
# swap out some correlated ones also
mod_sub1 <- update(mod_full, formula = ~. - understory_cover_sd_100_sc)
summary(mod_sub1)

mod_sub2 <- update(mod_full, formula = ~. - understory_cover_sd_100_sc - weather_wind_comb)
summary(mod_sub2)


### ASSIGN BEST-FIT MODEL (CHECK FOR RED FLAGS!), DO DIAGNOSTICS, AND SAVE ----
# Include mod0 and party models and best-fit single predictor models
# If it's not the BIC-supported model that is chosen, go back and drop things in different orders from the best-fit-so-far
# !!! Check the standard errors and coefficients--if very large SE or extreme coefficients, there is a problem. Perhaps even caused by a single location.
# For zero-inflated models (or any distribution NOT Poisson, check the additional distribution parameter for nonsensicality. For example, if best-fit is zero-inflated but the zero-inflation estimate is near zero, the problem is probably one data point, site, etc.)
# Check if top models all have similar estimates and CI's
# How much variability do the models actually explain? (maybe they all suck)
# Prioritize best-fit BIC, but check if coef similar with best-fit AIC and highest R2 models. If diagnostics really not good...
#- Okay if just a few diagnostics are wonky.
#- is it because a single site? (can try leverage)
#- cross-validate to see which model is best CV
#- try brms?
#- try increasingly complex models that are best-fit AIC until diagnostics are good--but only if it makes sense ecologically and there is not much difference across the estimates & CI's compared to slimmer model, the added predictors are significant or almost, and no issues with huge SE's or weird betas. 
# If dashboard shows especially low counts in 2021 for VICK, it may be necessary to include wind covariate (for better model diag.) even if not signif.


# Full model list <<<<<<< BRING TOGETHER ALL MOD_LISTS
mod_list <- tibble::lst(mod0pois, mod0comp, mod0pois_zi, mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod_full, mod_sub1, mod_sub2, mod_adj1, mod_adj2) 
(p <- performance::compare_performance(mod_list, rank = TRUE))

# Assign best-fit model, so can check diagnostics
year_mod <- mod_sub2
mod_dat <- dat 

### EVAL. MODEL FIT, THEN RETURN MISSING DATA RECORDS, IF POSSIBLE ----
# If covariates with missing data are not in the final model, return the relevant data records to the dataset, then rerun models and model comparisons again
# For PAAL, missing data covariates are temp and wind


### SAVE THE FULL MOD_LIST, MOD, & DAT
saveRDS(mod_list, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_yearmod_modlist.RDS")))
saveRDS(year_mod, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_yearmod.RDS")))
saveRDS(mod_dat, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_dat.RDS")))


## [IF NEEDED] CROSS-VALIDATE subset of models ----
#Include the BIC-best-supported. Default metric is MSE. Smaller MSE is better. Does this tend to go with the BIC-best-fit? Include others with more AIC support and statistically signif. predictors that help with diagnostics
mod_dat <- dat
(cv_yearmod_list <- cv::cv(cv::models(mod_sub2 = mod_sub2, mod_full = mod_full), data = mod_dat, ncores = 8, k = 10))
(cv_yearmod_loc_list <- cv::cv(cv::models(mod_sub2 = mod_sub2, mod_full = mod_full), data = mod_dat, clusterVariables = "location_name", ncores = 8, k = 10))
# SAVE THE CV LIST
saveRDS(cv_yearmod_list, here::here("tmb_mods_SIMPLE", paste0("tmb_", park, "_", spec, "_cv_yearmod_list.RDS")))
saveRDS(cv_yearmod_loc_list, here::here("tmb_mods_SIMPLE", paste0("tmb_", park, "_", spec, "_cv_yearmod_loc_list.RDS")))

# # Check also with LOO--this often hangs up
# (LOO_year_mod_list <- cv::cv(cv::models(mod13 = mod13, mod_sub3 = mod_sub3, mod_sub5 = mod_sub5), data = mod_dat, ncores = 8, k = "loo", start = TRUE))

## [IF NEEDED] CHECK FOR HIGH LEVERAGE LOCATIONS ----
## >>> CANNOT YET CHECK FOR INFLUENTIAL POINTS IN GLMMTMB MODELS, BUT THEY ARE WORKING ON IT... https://github.com/glmmTMB/glmmTMB/pull/1065
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
# [IF NEEDED] For influential POINTS ----
# mod_leverage_case <- influence_mixed(mod, groups=".case")

# For influential locations... DO THIS IF DIAGNOSTIC PROBLEMS
mod_leverage <- influence_mixed(mod, groups="location_name", ncores = 8)


# Convert to data frame, then heat map
x <- as.data.frame(mod_leverage$`fixed.effects[-location_name]`)
x_long <- x %>% rownames_to_column(var = "location_name") %>% tidyr::pivot_longer(data = ., cols = -location_name, names_to = "FE", values_to = "coef")

x_long_sub <- x_long
# x_long_sub <- x_long %>% dplyr::filter(!FE %like% "yr_c_f")

p <- ggplot(x_long_sub, aes(x = FE, y = location_name, fill = coef)) + geom_tile() + geom_text(aes(label = round(coef, 2)), size = 4, color = "white") # this heat map shows for each FE, what the estimated coefficient is if you exclude that location from the an
ggplotly(p)

# saveRDS(p, here::here("tmb_mods_SIMPLE", paste0("tmb_", park, "_", spec, "_100m_yearmod_levtable.RDS")))
saveRDS(p, here::here("tmb_mods_SIMPLE", paste0("tmb_", park, "_", spec, "_100m_trendmod_levtable.RDS")))

# [AS NEEDED] TRY ADJUSTED MODELS ----
# DO THIS ONLY IF DIAGNOSTICS ARE NOT GOOD!! COME BACK AND MAKE MODEL ADJUSTMENTS BASED ON DIAGNOSTICS... DROP UNUSUAL SITES, ADD INTERACTIONS, OTHER CHANGES AS NEEDED... ADD THEM TO THE FULL LIST AND REPEAT
# try dropping unusual sites
# NOTE: Cannot do model selection when data differ

mod_adj1 <- update(year_mod, family="compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))
summary(mod_adj1)

mod_adj2 <- update(mod_sub1, family="compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))
summary(mod_adj2)


# Compare models...
(p <- performance::compare_performance(tibble::lst(mod_full, mod_sub1, mod_sub2, mod_adj1, mod_adj2), rank = TRUE))

### DIAGNOSTICS ON BEST-FIT MODEL ----
# NOTE: https://github.com/glmmTMB/glmmTMB/issues/888 talks about glmmTMB w/DHARMA--doesn't have capability to simulate conditional on RE's, and may lead to apparent issues of autocorrelated residuals
# year_mod <- update(mod_adj1, data = dat %>% dplyr::filter(!location_name %in% c("VO08")))

mod <- mod_adj2 #year_mod
# mod <- trend_mod
mod_dat <- dat

# mod <- gldat# mod <- glmmTMB(formula = sum_indiv ~ yr_sc + (1 + yr_sc | location_name) + perc_opendev_200_sc, data = mod_dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))
# mod_dat <- dat %>% dplyr::filter(!location_name %in% c("VF36", "VO08")) %>% droplevels(.) #dat # dat %>% dplyr::filter(location_name != "VF31") %>% droplevels(.)

# mod_list <- readRDS(here::here("tmb_mods_SIMPLE", paste0("tmb_", park, "_", spec, "_100m_yearmod_modlist.RDS")))

summary(mod) # make sure no absurd standard errors
ggstats::ggcoef_table(mod, intercept = FALSE, table_witdhs = c(2, 1)) + labs(title = "XXX") # THESE ARE WALD CI'S!


### Performance checks ----
performance::check_model(mod, check = c("pp_check", "vif", "qq","linearity", "reqq"))
performance::check_predictions(mod, check_range = TRUE)
performance::check_collinearity(mod) # If interaction terms are included in a model, high VIF values are expected. Can run check on the model without interaction to make sure not a problem. High VIF doesn't bias the estimate but makes the 95%CI larger.
performance::check_autocorrelation(mod) # if autocorrelation detected, there may be unmodeled clustering of the data
performance::check_overdispersion(mod)
performance::check_zeroinflation(mod)

# If high VIF issue, check pair plots
GGally::ggpairs(mod$frame %>% dplyr::select(-yr_sc, -location_name))
# 
# # ## check residuals on own
# resid <- simulate_residuals(mod, quantile_function = qnorm, iterations = 500)
# df <- cbind(mod$frame, sim_resid = residuals(resid), fit = fitted(mod))
# ggplot(df, aes(x = fit, y = sim_resid, color = location_name)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + facet_wrap(~location_name) # plot residuals against fitted

 # >>>>>>>> PICK UP FROM HERE. EAME DHARMA RESULTS NOT GOOD. MOD_SUB2 IS "BEST-FIT" BUT THERE ARE ISSUES W/HETEROSCEDASTICITY? THIS IS WITH USING COMPLETE CASES ONLY.


### DHARMA CHECKS ----
# NOTE: IF the diagnostics are "wonky" (no obvious pattern, but jumping around somewhat randomly) then try running it in brms and checking diagnostics there (are there just some odd sites, etc.?)

simres <- DHARMa::simulateResiduals(fittedModel = mod, n = 1000, refit = FALSE, plot = TRUE) # re-simulate all levels, to test model structure as a whole. The default is resimulating on all levels (unconditional sims) instead of simulating conditional on the fitted random effects.
# simres <- DHARMa::simulateResiduals(fittedModel = mod, n = 1000, refit = FALSE, re.form = NULL, plot = TRUE) #re-simulate conditional on the fitted random effects, i.e., the RE's are not re-simulated. More sensitive for detecting problems, but it doesn't test all possible options.

DHARMa::testZeroInflation(simres)

DHARMa::testDispersion(simres) # test for over/under dispersion
DHARMa::testOutliers(simres, type = "bootstrap")

# Plot residuals vs each independent factor ----
# Identify the ones with problems--were any of these identified as important earlier?
# Is weirdness due to just one or two odd sites?
DHARMa::plotResiduals(simres, mod_dat$yr_sc)
DHARMa::plotResiduals(simres, as.factor(mod_dat$yr_c_f))

DHARMa::plotResiduals(simres, mod_dat$perc_forest_50, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$perc_forest_100, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$perc_opendev_50, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$perc_opendev_100, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$prop_understory_cover_50, rank  = FALSE) #
DHARMa::plotResiduals(simres, mod_dat$prop_understory_cover_100, rank  = FALSE) # 
DHARMa::plotResiduals(simres, mod_dat$understory_cover_sd_100, rank  = FALSE) # <<

DHARMa::plotResiduals(simres, mod_dat$weather_wind_comb, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$weather_temperature_cs, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$julian_prop_sc, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$hrs_since_rise_sc, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$researcher, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$first_yr, rank  = FALSE)

# # If there are diagnostic problems below, not a big deal--unless a lot of problems with clear trends
# DHARMa::plotResiduals(simres, mod_dat$perc_forest_50, rank  = FALSE)
# DHARMa::plotResiduals(simres, mod_dat$perc_forest_200, rank  = FALSE) # 
# 
# DHARMa::plotResiduals(simres, mod_dat$perc_opendev_50, rank  = FALSE)
# DHARMa::plotResiduals(simres, mod_dat$perc_opendev_200, rank  = FALSE)
# 
# DHARMa::plotResiduals(simres, mod_dat$prop_understory_cover_200, rank  = FALSE) # <<
# 
# DHARMa::plotResiduals(simres, mod_dat$understory_cover_sd_50, rank  = FALSE)
# DHARMa::plotResiduals(simres, mod_dat$understory_cover_sd_200, rank  = FALSE)
# 
# 
# 
# # ### Append Dharma residuals to model data ----
# dat_append <- cbind(mod_dat, predicted = DHARMa::getFitted(mod), resid = DHARMa::getResiduals(mod))
# dat_append$scaled_resid = scales::rescale(dat_append$resid, to = c(0,1), from = range(dat_append$resid, na.rm = TRUE))
# dat_append$rank_pred = datawizard::ranktransform(dat_append$understory_cover_sd_100)
# plotly::ggplotly(ggplot(dat_append, aes(x=understory_cover_sd_100, y = scaled_resid, group = 1, text = paste0(location_name, "<br>Observed: ", sum_indiv, "<br>PredictedResponse: ", round(predicted, 2), "<br>Cov: ", understory_cover_sd_100, "<br>Resid: ", round(resid,2)))) + geom_point() + geom_quantile(), tooltip = "text")

### SAVE THE FINAL 100M 2-VISIT YEAR MODEL ----
saveRDS(year_mod, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_yearmod.RDS")))
saveRDS(mod_dat, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_dat.RDS")))

### CHECK THE PREDICTION PLOT ----
# NOTE THESE ARE WALD CI'S!!!

## ggemmeans to get year predictions <<<<<<<<<<<<<<<<<<<<<<< CHANGE THIS IN FINAL SUMMARIES
# ggeffects::predict_response(mod, terms = c("yr_c_f", "subunit"), margin = "mean_reference", type = "zero_inflated")) # > Note that to get emmeans, 'margin = "marginalmeans"' but then can't do zero-inflation with that
# (temp_tmb_pred_yearmod <- ggemmeans(mod_tmb_yearmod, terms=c("yr_c_f"), type = "fixed") %>% # # type = ifelse(tmb_zi==TRUE, "fe.zi", "fixed"))) #### ZERO-INFLATION PART DOESN' WORK!!!

year_mod_zi <- ifelse(year_mod$modelInfo$allForm$ziformula=="~1", TRUE, FALSE)

# tmb_pred_yearmod <- ggpredict(year_mod, terms=c("yr_c_f"), type = ifelse(year_mod_zi==TRUE, "zero_inflated", "fixed")) # This is same as below

# # CHeck interaction effects!! [IF APPLICABLE]
# tmb_pred_yearmod_interact <- ggeffects::predict_response(year_mod, terms=c("perc_opendev_50_sc", "odd3"), margin = "mean_mode", type = ifelse(year_mod_zi==TRUE, "zero_inflated", "fixed")) %>% plot()
# saveRDS(tmb_pred_yearmod_interact, here::here("tmb_mods_SIMPLE", paste0("tmb_", park, "_", spec, "_100m_yearmod_interaction_plot.RDS")))

# Check year estimates
(tmb_pred_yearmod <- ggeffects::predict_response(year_mod, terms=c("yr_c_f [all]"), margin = "mean_mode", type = ifelse(year_mod_zi==TRUE, "zero_inflated", "fixed")))

(tmb_pred_yearmod %<>% 
    as.data.frame() %>%
    janitor::clean_names() %>%
    dplyr::mutate(yr_f = as.numeric(as.character(x))) %>%
    dplyr::select(
      group,
      yr_f,
      yr_predicted = predicted,
      yr_conf_low = conf_low,
      yr_conf_high = conf_high))
tmb_pred_yearmod$yr <- tmb_pred_yearmod$yr_f+(range(dat$yr) %>% median())

ggplot(data = tmb_pred_yearmod) + 
  geom_errorbar(
    mapping = aes(x = yr, y = yr_predicted, ymin = yr_conf_low, ymax = yr_conf_high, color = group), width = 0) +
  geom_point(
    mapping = aes(x = yr, y = yr_predicted, color = group), size = 2.5) +
  facet_grid(rows = vars(group))


### TRY BRMS ----
# Since we have standardized data, Normal(0, 1) is reasonable for reg coefs
# If it doesn't work, try on personal laptop with cmdstanr
mod_brms <- brms::brm(
  formula = sum_indiv ~ yr_c_f + (1 | location_name) + perc_forest_100_sc +  
    prop_understory_cover_50_sc + weather_temperature_cs + I(perc_forest_100_sc^2),
  data = mod_dat,
  # family = brmsfamily("poisson"),
  family = brmsfamily("com_poisson"),
  prior = c(prior(normal(0,1), class = b)),
  # sample_prior = "only", # <<<<<<<<<<< CHECK PRIORS
  iter = 4000, warmup = 1000, chains = 1, cores = 8)

pp_check(mod_brms, ndraws = 200)
mod_brms <- add_criterion(mod_brms, criterion = c("loo"))
loo(mod_brms)
saveRDS(mod_brms, here::here("brms_mods", paste0("brms_", park, "_", spec, "_100m_yearmod.RDS")))
saveRDS(mod_dat, here::here("brms_mods", paste0("brms_", park, "_", spec, "_100m_yearmod_dat.RDS")))


### TREND MODELS ----
# For random slopes...+ (1 + yr_sc | location_name)). You can actually omit the '1+' b/c it is assumed--to manually remove intercept RE, would do '0+'. In all models, include random intercept AND slope. For some models, it may also make sense to have a fixed effect interaction of habitat type with year. The base GLMM model is therefore:  Overall trend (yr_sc) plus year to year fluctuations concordant across locations (1|yr_f) plus random variation among locations in location-specific trends ( 1 + yr_sc | location_name). Also, note   https://academic.oup.com/esr/article-abstract/35/2/258/5306121?redirectedFrom=fulltext
# >>>>> CHECK IF I SHOULD HAVE ADDED DISPFORMULA =~1, OR IF THAT IS DEFAULT
# Start with the best-fit year model, then try models with hab*yr interaction if that seems appropriate and use diagnostics, etc. to see if should drop/add other covariates
# Try dropping covariates that in year model were only driven by one or two sites


# mod_dat <- dat %>% dplyr::filter(!location_name %in% c("VO08")) %>% droplevels()
# 
# p <- ggplot(dat, aes(x = yr, y = sum_indiv, color = hab_type_100, fill = hab_type_100, group_by = location_name)) +
#          stat_smooth(method = "lm", se = TRUE)
# plotly::ggplotly(p)
# 
# raw_trend <- lm(formula = sum_indiv ~ yr_sc*location_name, data = mod_dat)
# trend_coef <- broom::tidy(raw_trend)
# trend_df <- trend_coef[grep("yr_sc:location_name", trend_coef$term), ]
# trend_df$term <- gsub(pattern = "yr_sc:location_name", replacement = "", trend_df$term) # then color by estimate



mod$call
### BASE TREND MODEL & CHECK RE VARIANCES----
# BASE model is trend with year trend, year RE and location intercept&slope RE 
trend_mod_base_compois <- glmmTMB(sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc | location_name), data = mod_dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8)) 
summary(trend_mod_base_compois)

trend_mod_base_poisson <- glmmTMB(sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc | location_name), data = mod_dat, family = "poisson")
summary(trend_mod_base_poisson)
# 
# trend_mod_base_pois_zi <- glmmTMB(sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc | location_name), data = mod_dat, family = "poisson", ziformula = ~1)
# summary(trend_mod_base_pois_zi)

# ### SIMPLIFY IF MODEL CONVERGENCE PROBLEMS ----
# ## [ONLY IF MODELS WON'T CONVERGE] TRY DROPPING SOME RANDOM EFFECTS ----# IF RE variance is practically zero, consider dropping it
# trend_mod_base2 <- glmmTMB(sum_indiv ~ yr_sc + (1 + yr_sc | location_name), data = mod_dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8)) 
# summary(trend_mod_base2)
# 
# # IF the RE's are highly correlated AND if location slope RE is nearly zero, drop the slope
# # plot the relationship btwn intercept and slope RE
# 
# trend_mod_base3 <- glmmTMB(sum_indiv ~ yr_sc + (1|yr_c_f) + (1|location_name), data = mod_dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8)) 
# summary(trend_mod_base3)
# 
# # If the RE's are highly correlated BUT LOCATION SLOPE VARIANCE IS NOT NEARLY ZERO drop the correlation
# trend_mod_base4 <- glmmTMB(formula = sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc ||location_name), data = mod_dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))
# summary(trend_mod_base4)
# 
# trend_mod_base_pois_zi_wind <- glmmTMB(sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc | location_name), data = mod_dat, family = "poisson", ziformula = ~weather_wind_comb)
# summary(trend_mod_base_pois_zi_wind)
# 
# anova(trend_mod_base4, trend_mod_base3)

# Compare base models ----
# Exponentiate to get % of population in next year. exp(trend_coef)^(1:10) to see % over next 10 years. Need to bootstrap CI to get correct CI on trend.
trend_mod_base_list <- tibble::lst(trend_mod_base_poisson, trend_mod_base_compois)
performance::compare_performance(trend_mod_base_list, rank = TRUE)
# ggstats::ggcoef_compare(tibble::lst(trend_mod0, trend_mod_base4))

# # Check for trend-habitat interaction
# (cv_trend_base_mod <- cv::cv(models(trend_mod0=trend_mod0, trend_mod0_Xhab=trend_mod0_Xhab), data = mod_dat, reps = 5))

trend_mod_base_RE <- trend_mod_base_poisson # <<<<<<<<<<<<< ENTER HERE, BIC-SUPPORTED MODEL

## ADD COVARIATES FROM BEST-FIT YEARMOD, START WILL FULL MODEL THEN SUBSET ----
trend_mod_full <- update(trend_mod_base_RE, formula = ~. + weather_temperature_cs)
summary(trend_mod_full)

trend_mod_adj1 <- update(trend_mod_full, formula = ~. + weather_temperature_cs - hrs_since_rise_sc)
summary(trend_mod_adj1)

trend_mod_sub1 <- update(trend_mod_full, formula = ~. - weather_temperature_cs)
# summary(trend_mod_sub1)
# 
# 
# anova(trend_mod_sub1, trend_mod_full)
trend_mod_list <- tibble::lst(trend_mod_base_poisson, trend_mod_full)
performance::compare_performance(trend_mod_list, rank = TRUE)
# ggstats::ggcoef_compare(tibble::lst(trend_mod_full, trend_mod_sub1))
# 
# 
# ggstats::ggcoef_table(trend_mod_full, intercept = FALSE, table_witdhs = c(2, 1)) + labs(title = "XXX") # make sure no crazy estimates or CI's

## [AS NEEDED] SET PRIORS ----
trend_mod2_prior <- update(trend_mod2, priors = c(prior(normal(0,1), class = "fixef")))

trend_mod_prior <- glmmTMB(sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc | location_name) + hrs_since_rise_sc + understory_cover_sd_200_sc + weather_wind_comb, data = mod_dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8), priors = c(prior(normal(0,1), class = "fixef"))) 



## [AS NEEDED] HABITAT INTERACTION WITH YEAR-SLOPE----
trend_mod_base_hab50 <- update(trend_mod0, formula = ~.+ yr_sc*hab_type_50_comb)
summary(trend_mod_base_hab50)

trend_mod_base_hab100 <- update(trend_mod_sub1, formula = ~.+ yr_sc*hab_type_100_comb)
summary(trend_mod_base_hab100) 

trend_mod_base_hab200 <- update(trend_mod_base_RE, formula = ~.+ yr_sc*hab_type_200_comb)
summary(trend_mod_base_hab200)


### [IF NEEDED] EXPLORE RANDOM EFFECTS ----
# Plot correlation btwn location intercept RE & slope RE
trend_mod <- trend_mod_adj1

re <-ranef(trend_mod)[1]$cond$location_name %>%
  janitor::clean_names()
re
ggplot(re, aes(x = intercept, y = yr_sc)) + geom_point()
ggplot(re, aes(intercept)) + geom_histogram()
ggplot(re, aes(yr_sc)) + geom_histogram()

### ASSIGN TREND_MOD ----
trend_mod <- trend_mod_full
mod <- trend_mod

### SAVE TRENDMOD RESULTS ----
saveRDS(trend_mod_list, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_trendmod_modlist.RDS")))
saveRDS(trend_mod, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_trendmod.RDS")))
saveRDS(mod_dat, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_trendmod_dat.RDS")))

### TRY BRMS ----
# Brms can be useful for identifying highly influential points and seeing if those points are responsible for model diagnostic issues.
mod_dat <- dat

trendmod_brms <- brms::brm(
  formula =  sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc|location_name) + hrs_since_rise_sc + perc_forest_100_sc + weather_wind_comb + prop_understory_cover_50_sc,
  # formula =  sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc ||location_name) + hrs_since_rise_sc + understory_cover_sd_200_sc,
  # formula =  sum_indiv ~ yr_sc + (1 | yr_c_f) + (1|location_name) + hrs_since_rise_sc + understory_cover_sd_200_sc + weather_wind_comb,
  # family = brmsfamily("poisson"),
  family = brmsfamily("com_poisson"),
  prior = c(prior(normal(0,1), class = b)),
  # sample_prior = "only", # <<<<<<<<<<< CHECK PRIORS
  data = mod_dat, 
  iter = 4000, warmup = 1000, chains = 1, cores = 8)
pp_check(trendmod_brms, ndraws = 500)
trendmod_brms <- add_criterion(trendmod_brms, criterion = c("loo"))

saveRDS(trendmod_brms, here::here("brms_mods_SIMPLE", paste0("brms_", park, "_", spec, "_100m_trendmod.RDS"))) # <<<<<

# # If there are INFLUENTIAL POINTS, they are usually the ones that are particularly high or low in early or late years
(trendmod_loo <- loo(trendmod_brms)) %>% plot(., label_points = TRUE)

saveRDS(trendmod_loo, here::here("brms_mods_SIMPLE", paste0("brms_", park, "_", spec, "_100m_trendmod_loo.RDS"))) # <<<<
# dat_loo <- cbind(dat[, c("location_name", "event_date", "sum_indiv")], mod_loo$diagnostics$pareto_k)
# View(dat_loo)




## Evaluate BRMS models (TRENDMOD_BRMS OR MOD_BRMS) ----
#### USER INPUT ON MODEL CHECKS BELOW... <<<<<<<<<<<<<<<<<<<<
mod_brms = trendmod_brms # <<<<<<<<<<<<< SELECT BEST MODEL
summary(mod_brms) # slightly shorter summary of model results

rstan::check_hmc_diagnostics(mod_brms$fit)
bayes_R2(mod_brms)

# Coefficient plot
mcmc_plot(mod_brms) # quick coefficient plot


# Model convergence check
plot(mod_brms) # plots showing posterior distribution of parameters, and trace plots


# Posterior predictive checks
# pp_check(mod_brms, type = "dens_overlay", ndraws = 200)
pp_check(mod_brms, type = "stat_2d", ndraws = 200)
# pp_check(mod_brms, type = "loo_pit_qq", ndraws = 200)
pp_check(mod_brms, type = "rootogram", ndraws = 200)
pp_check(mod_brms, type = "bars", ndraws = 200)


pp_check(mod_brms, type = "scatter_avg_grouped", group = "yr_c_f", ndraws = 200) + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
pp_check(mod_brms, type = "bars_grouped", group = "yr_c_f", ndraws = 200)
pp_check(mod_brms, type = "bars_grouped", group = "location_name", ndraws = 200)

# If categorical predictors...
pp_check(mod_brms, type = "bars_grouped", group = "weather_wind_comb", ndraws = 200) # <<<<<<<< CHANGE, DEPENDING ON BEST MODEL


## Plot conditional effects
# plot(conditional_effects(mod_brms, effects ="understory_cover_sd_200_sc:prop_understory_cover_50_sc", points=TRUE, point_args = list(width = 0.3, height = 0.1, alpha = 0.4)))

# Save conditional effects plots as list, for use with ggplot functions
(p <- plot(conditional_effects(mod_brms), points=TRUE, point_args = list(width = 0.3, height = 0.1, alpha = 0.4), plot = FALSE))

# saveRDS(p, here::here("brms_mods_SIMPLE", paste0("brms_", park, "_", spec, "_100m_condeffects.RDS"))) # <<<
saveRDS(p, here::here("brms_mods_SIMPLE", paste0("brms_", park, "_", spec, "_100m_trendmod_condeffects.RDS"))) # <<<

## DHARMA
simres <- DHARMa.helpers::dh_check_brms(mod_brms)
## NOW GO TO DHARMA CODE



### PLOT glmmTMB PREDICTIONS OF FINAL TRENDS MODEL ----
# Check the estimated trend
em=emmeans::emtrends(trend_mod, specs = c("yr_sc"), var = c("yr_sc"))
# em=emmeans::emtrends(trend_mod, specs = c("yr_sc", "hab_type_100_comb"), var = "yr_sc") # trends at mean year. Ignores a 'type' argument, it seems
summary(em, infer=c(TRUE,TRUE),null=0) %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  dplyr::mutate(signif = p_value < 0.05)

(temp_tmb_pred <- ggemmeans(trend_mod, terms=c("yr_sc[all]"), type = "fixed"))


### PLOT TREND WITH YEAR ESTIMATES ----
tmb_zi_trendmod <- ifelse(trend_mod$modelInfo$allForm$ziformula=="~1", TRUE, FALSE)

# (temp_trendmod_mmeans <- ggeffects::predict_response(trend_mod, terms=c("yr_sc", "hab_type_100_comb"), margin = "marginalmeans", type = ifelse(tmb_zi_trendmod==TRUE, "zero_inflated", "fixed"))) # this is the same as emmeans. For factors, it weights them equally regardless of how many are actually in each factor level
# (temp_trendmod_mmeans <- ggeffects::predict_response(trend_mod, terms=c("yr_sc"), margin = "marginalmeans", type = ifelse(tmb_zi_trendmod==TRUE, "zero_inflated", "fixed")))

# (temp_trendmod_mmode <- ggeffects::predict_response(trend_mod, terms=c("yr_sc", "hab_type_100_comb"), margin = "mean_mode", type = ifelse(tmb_zi_trendmod==TRUE, "zero_inflated", "fixed")))    
(temp_trendmod_mmode <- ggeffects::predict_response(trend_mod, terms=c("yr_sc [all]"), margin = "mean_mode", type = ifelse(tmb_zi_trendmod==TRUE, "zero_inflated", "fixed")))  

temp_trendmod_mmode %<>% 
  as.data.frame 

FuncFormatTrendDat <- function(pred, dat) {
  pred_out <- pred %>% 
    as.data.frame()
  
  
  pred_out %<>%  
    janitor::clean_names() %>%
    dplyr::mutate(yr_sc = as.numeric(x))
  
  if("group" %in% names(pred_out)) {
    
    pred_out %<>%
      dplyr::select(
        group,
        yr_sc,
        trend_predicted = predicted,
        trend_conf_low = conf_low,
        trend_conf_high = conf_high)
  } else {
    pred_out %<>%
      dplyr::select(
        yr_sc,
        trend_predicted = predicted,
        trend_conf_low = conf_low,
        trend_conf_high = conf_high)
  }
  
  return(pred_out)
}

# tmb_pred_trendmod_mmeans <- FuncFormatTrendDat(pred = temp_trendmod_mmeans, dat)
tmb_pred_trendmod_mmode <- FuncFormatTrendDat(pred = temp_trendmod_mmode, dat) # <<<<<<<<<<<<< FIX--THIS ISN'T QUITE RIGHT AS A COMPARISON TO YEARMOD
tmb_pred_trendmod_mmode %<>% dplyr::left_join(dat %>% dplyr::select(yr, yr_sc) %>% dplyr::mutate(yr_sc = round(yr_sc, 2)) %>% distinct(), by = "yr_sc")

(p <- ggplot() + 
    geom_errorbar(data = tmb_pred_yearmod, aes(x = yr, y = yr_predicted, ymin = yr_conf_low, ymax = yr_conf_high), width = 0, color = "blue") +
    geom_point(data = tmb_pred_yearmod, aes(x = yr, y = yr_predicted), size = 2.5, color = "blue") +
    # geom_ribbon( # glmmTMB TREND RESULTS AS RIBBON
    #   data = tmb_pred_trendmod_mmeans, aes(x = yr, y = trend_predicted, ymin = trend_conf_low, ymax = trend_conf_high),
    #   alpha = 0.1, fill = "orange") +
    # geom_line(data =tmb_pred_trendmod_mmeans, aes(x = yr, y= trend_predicted), color = "orange", linetype = "dashed") +
    geom_ribbon(
      data = tmb_pred_trendmod_mmode, aes(x = yr, y = trend_predicted, ymin = trend_conf_low, ymax = trend_conf_high),
      alpha = 0.1, fill = "orange") +
    geom_line(data =tmb_pred_trendmod_mmode, aes(x = yr, y= trend_predicted), color = "orange", linetype = "dashed") +
    scale_x_continuous(breaks = unique(sort(dat$yr))) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title = paste0("Trends in relative abundance of ", spec, " at ", park), x = "Year", y = "Index of abundance (# of birds/point count)") +
    theme_bw() +
    facet_grid(rows = vars(group)))

saveRDS(p, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_100m_trendplot.RDS")))