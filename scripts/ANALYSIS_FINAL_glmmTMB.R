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
spec = "REVI" # go through each of the top species...

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

dat <- FuncFormatFullObs(park = park, spec = spec, limit_100m = TRUE, drop_sites = FALSE) # <<<<<<<<<<<<< DROP_SITES OR NOT
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
sum(dat$sum_indiv) # total detections in final dataset

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

if(park == "VICK"){
  dat %<>% dplyr::mutate(odd3 = factor(ifelse(prop_understory_cover_50<0.4, "odd3_locs", "all_other_locs"))) # This is useful for some bird species... the 'odd3' locations
}


dat %>% dplyr::group_by(subunit) %>% summarize(tot = sum(sum_indiv))

# >>>>> Check for complete separation
dat %>% dplyr::group_by(weather_wind_comb) %>% summarize(mean(sum_indiv))
dat %>% dplyr::group_by(odd3) %>% summarize(mean(sum_indiv))
dat %>% dplyr::group_by(yr_c_f) %>% summarize(mean(sum_indiv))

# >>>>> Check for covariate patterns that may explain unusually low or high years
table(dat$yr_c_f, dat$weather_wind_comb)
table(dat$yr_c_f, dat$odd3)

### EDA ----

## >>>>> FIRST, LOOK ON DASHBOARD FOR YEAR TRENDS, ODD YEARS (OR YRS WITH NO DETECTIONS), UNUSUAL HIGH NUMBERS (E.G., LIKELY INFLUENTIAL POINTS FOR TREND MODEL), NUMBER OF LOCATIONS WITH NO DETECTIONS EVER, PREDICTOR RELATIONSHIPS. ALSO LOOK UP HABITAT FOR THE SPECIES TO SET EXPECTATIONS ABOUT PREDICTORS. LOOK FOR EVIDENCE OF COMPLETE SEPARATION.

# >>>> CHECK WHICH BIRDS IF ANY SHOW DIFFERENT PATTERN FOR HABITAT COVER @ DIFFERENT SCALES

# PARTYKIT: Classification tree. Think of it as piecewise linear. Gives us different information from linear regression--may be useful for identifying interactions and polynomial predictors. Sometimes it identifies a fork, but the numbers don't differ much on either side of the fork or the sample size is really small for one fork or the cutoff for the fork is near a boundary. Consider those with a grain of salt.

## SAAN & GUIS--add subunit

# dat %<>% dplyr::filter(subunit=="MI") # <<<<< SAAN ONLY, AS NEEDED

## as needed... researcher + first_yr (BITH, PAAL, GUIS?)


dat_tree <- partykit::ctree(sum_indiv ~ prop_understory_cover_50_sc + prop_understory_cover_100_sc, data = dat) # 
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~odd3 + prop_understory_cover_50_sc + prop_understory_cover_100_sc, data = dat) # 
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~prop_understory_cover_50_sc + prop_understory_cover_100_sc +  perc_forest_100_sc +perc_opendev_100_sc +  understory_cover_sd_100_sc + weather_wind_comb + weather_temperature_cs + julian_prop_sc + hrs_since_rise_sc , data = dat) # 
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
summary(mod0pois_zi)
plogis(-21.62)

# AICc to compare
mod_list <- tibble::lst(mod0pois, mod0comp, mod0pois_zi)# mod0nb) #, mod0pois_zi_loc)#, mod0nb)
performance::compare_performance(mod_list, rank = TRUE, metrics = "common") # CANNOT INCLUDE MODELS THAT DROP DATA HERE, THE DATASET HAS TO BE THE SAME FOR ALL COMPARED MODELS!

# Assign best one to mod0 (prioritize BIC)

mod0 <- mod0comp

# Add covariates

## ONE BY ONE... ----
# TRy quadratic if EDA suggests that is more appropriate
# >>> CHECK FOR VERY STRANGE ESTIMATES!


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

mod10 <- update(mod0, formula = ~. + odd3)
summary(mod10) # 

# Compare models
mod_list <- tibble::lst(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9) #, mod10) # >>>>> if odd3 is overwhelming, retry it without
(perf <- performance::compare_performance(mod_list, rank = TRUE, metrics = "common")) # How much variability does the model actually explain?

## FULL MODEL... ----
# Include all variables that were statistically signif (or close) in single models OR in partykit, possibly ordered by compare_performance. For each, use the scale that is "most" signif.
# Include polynomials if signif in regression tree

mod_full <- update(mod0, formula = ~. + odd3 + perc_forest_100_sc + prop_understory_cover_50_sc)

## NOW SUBSETS REMOVING THE LEAST IMPORTANT ONE BY ONE, THEN IN GROUPS... 
# remove one at a time, check the model, remove next least important
# swap out some correlated ones also
# >> IF ODD3 WAS IMPORTANT, TRY MODELS WITH UNDERSTORY50 QUADRATIC AS SUBSTITUTE

mod_sub1 <- update(mod_full, formula = ~. - weather_temperature_cs)
summary(mod_sub1)

mod_sub2 <- update(mod_full, formula = ~. - hrs_since_rise_sc - weather_temperature_cs)
summary(mod_sub2)

mod_sub3 <- update(mod_full, formula = ~. - hrs_since_rise_sc)
summary(mod_sub3)

mod_sub4 <- update(mod_full, formula = ~. - understory_cover_sd_100_sc - weather_temperature_cs)
summary(mod_sub4)

# >>>> TRY UNDERSTORY PREDICTOR AS ODD3 SUBSTITUTE <<<<
mod_adj1 <- update(mod_sub1, formula = ~. + prop_understory_cover_50_sc)
summary(mod_adj1)

mod_adj2 <- update(mod_sub1, formula = ~. + prop_understory_cover_50_sc + I(prop_understory_cover_50_sc^2))
summary(mod_adj2)

mod_adj3 <- update(mod_sub1, formula = ~. + odd3)
summary(mod_adj3)

mod_adj4 <- update(mod_sub1, formula = ~. + odd3 - understory_cover_sd_100_sc)
summary(mod_adj4)

# >>>>>>>>>>>> TESTER! <<<<<<<<<<<<<<<<<<<<
mod_list <- tibble::lst(mod2, mod1, mod_full, mod_sub1, mod_sub2, mod_sub3, mod_sub4, mod_adj1, mod_adj2, mod_adj3, mod_adj4)
(p <- performance::compare_performance(mod_list, rank = TRUE, metrics = "common"))
plot(p)
ggstats::ggcoef_compare(tibble::lst(mod_sub1, mod_sub5)) # do this for subset of models


### ASSIGN BEST-FIT MODEL (CHECK FOR RED FLAGS!) AND SAVE ----
# Use BIC as primary criterion, but if BIC is similar use the model with better AIC and/or R^2
# !!! Check the standard errors and coefficients--if very large SE or extreme coefficients, there is a problem. Perhaps even caused by a single location.
# For zero-inflated models (or any distribution NOT Poisson, check the additional distribution parameter for nonsensicality. For example, if best-fit is zero-inflated but the zero-inflation estimate is near zero, the problem is probably one data point, site, etc.)
# For VICK, check EDA for an "odd3" effect that should be modeled--especially if the best-fit (w/o odd3 predictor) is a very complicated model> BUT>>>> better to use a quadratic understory_50 if can b/c it is more meaningful
# Check if top models all have similar estimates and CI's
# How much variability do the models actually explain? (maybe they all suck)
# Prioritize best-fit BIC, but check if coef similar with best-fit AIC and highest R2 models. If diagnostics really not good...
  #- is it because a single site? (can try leverage)
  #- cross-validate to see which model is best CV
  #- try brms?
  #- try increasingly complex models that are best-fit AIC until diagnostics are good--but only if it makes sense ecologically and there is not much difference across the estimates & CI's compared to slimmer model, the added predictors are significant or almost, and no issues with huge SE's or weird betas. 


### SAVE THE FULL MOD_LIST----
# Full model list <<<<<<< BRING TOGETHER ALL MOD_LISTS
mod_list <- tibble::lst(mod0pois, mod0comp, mod0pois_zi, mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod_full, mod_sub1, mod_sub2, mod_sub3, mod_sub4, mod_adj1, mod_adj2, mod_adj3, mod_adj4) 
# (p <- performance::compare_performance(mod_list, rank = TRUE))
saveRDS(mod_list, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_100m_yearmod_modlist.RDS")))

### SAVE THE FINAL 100M 2-VISIT YEAR MODEL----
year_mod <- mod_adj4 # <<<<<<<<<
mod_dat <- dat 

saveRDS(year_mod, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_100m_yearmod.RDS")))
saveRDS(mod_dat, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_100m_dat.RDS")))


## [IF NEEDED] CROSS-VALIDATE subset of models ----
#Include the BIC-best-supported. Default metric is MSE. Smaller MSE is better. Does this tend to go with the BIC-best-fit? Include others with more AIC support and statistically signif. predictors that help with diagnostics
mod_dat <- dat
(cv_yearmod_list <- cv::cv(cv::models(mod_sub2 = mod_sub2, mod_full = mod_full), data = mod_dat, ncores = 8, k = 10))
(cv_yearmod_loc_list <- cv::cv(cv::models(mod_sub2 = mod_sub2, mod_full = mod_full), data = mod_dat, clusterVariables = "location_name", ncores = 8, k = 10))
# SAVE THE CV LIST
saveRDS(cv_yearmod_list, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_cv_yearmod_list.RDS")))
saveRDS(cv_yearmod_loc_list, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_cv_yearmod_loc_list.RDS")))

# # Check also with LOO--this often hangs up
# (LOO_year_mod_list <- cv::cv(cv::models(mod13 = mod13, mod_sub3 = mod_sub3, mod_sub5 = mod_sub5), data = mod_dat, ncores = 8, k = "loo", start = TRUE))

## [IF NEEDED] CHECK FOR HIGH LEVERAGE LOCATIONS ----
## >>> CANNOT YET CHECK FOR INFLUENTIAL POINTS IN GLMMTMB MODELS, BUT THEY ARE WORKING ON IT... https://github.com/glmmTMB/glmmTMB/pull/1065
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
# For influential POINTS...
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

saveRDS(p, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_100m_yearmod_levtable.RDS")))
# saveRDS(p, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_100m_trendmod_levtable.RDS")))

## [AS NEEDED] TRY ADJUSTED MODELS ----
# DO THIS ONLY IF DIAGNOSTICS ARE NOT GOOD!! COME BACK AND MAKE MODEL ADJUSTMENTS BASED ON DIAGNOSTICS... DROP UNUSUAL SITES, ADD INTERACTIONS, OTHER CHANGES AS NEEDED... ADD THEM TO THE FULL LIST AND REPEAT
# try dropping unusual sites
# NOTE: Cannot do model selection when data differ

mod_adj3 <- update(mod_sub1, family = poisson, ziformula = ~1)
summary(mod_adj3)


## [AS NEEDED] DROP SITES OR OTHER ADJUSTMENTS FROM BEST-FIT MODEL ----
# Do this if DHARMA results suggest one or a few sites are driving odd diagnostics
# mod_drop1 <- update(year_mod, data = dat %>% dplyr::filter(!location_name %in% c("VF36", "VM07")) %>% droplevels(.)) # mod_sub4, but drop site with the highest single survey count (N=4)


### GLMMTMB DIAGNOSTICS ON BEST-FIT MODEL ----
# NOTE: https://github.com/glmmTMB/glmmTMB/issues/888 talks about glmmTMB w/DHARMA--doesn't have capability to simulate conditional on RE's, and may lead to apparent issues of autocorrelated residuals
# year_mod <- update(mod_adj1, data = dat %>% dplyr::filter(!location_name %in% c("VO08")))
# BIC favors excluding interaction but R2 and AIC and DHARMA diagnostics much better with interaction; coefficients similar for both, but interaction is signif
mod <- year_mod
mod <- trend_mod
mod_dat <- dat


summary(mod) # make sure no absurd standard errors
ggstats::ggcoef_table(mod, intercept = FALSE, table_witdhs = c(2, 1)) + labs(title = "XXX") # THESE ARE WALD CI'S!


### GLMMTMB PERFORMANCE CHECKS ----
performance::check_model(mod, check = c("pp_check", "vif", "qq","linearity", "reqq"))
performance::check_predictions(mod, check_range = TRUE)
performance::check_collinearity(mod) # If interaction terms are included in a model, high VIF values are expected. Can run check on the model without interaction to make sure not a problem. High VIF doesn't bias the estimate but makes the 95%CI larger.
performance::check_autocorrelation(mod) # if autocorrelation detected, there may be unmodeled clustering of the data
performance::check_overdispersion(mod)
performance::check_zeroinflation(mod)

# ## check residuals on own
# resid <- simulate_residuals(mod, quantile_function = qnorm, iterations = 500)
# df <- cbind(mod$frame, sim_resid = residuals(resid), fit = fitted(mod))
# ggplot(df, aes(x = fit, y = sim_resid, color = location_name)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + facet_wrap(~location_name) # plot residuals against fitted


### DHARMA CHECKS ----
# NOTE: IF the diagnostics are "wonky" (no obvious pattern, but jumping around somewhat randomly) then try running it in brms and checking diagnostics there (are there just some odd sites, etc.?)

# simres <- DHARMa::simulateResiduals(fittedModel = mod, n = 1000, refit = FALSE, plot = TRUE) # re-simulate all levels, to test model structure as a whole. The default is resimulating on all levels (unconditional sims) instead of simulating conditional on the fitted random effects.
# simres <- DHARMa::simulateResiduals(fittedModel = mod, n = 1000, refit = FALSE, re.form = NULL, plot = TRUE) #re-simulate conditional on the fitted random effects, i.e., the RE's are not re-simulated. More sensitive for detecting problems, but it doesn't test all possible options.

DHARMa::testZeroInflation(simres)

DHARMa::testDispersion(simres) # test for over/under dispersion
DHARMa::testOutliers(simres, type = "bootstrap")

# Plot residuals vs each independent factor
# Identify the ones with problems--were any of these identified as important earlier?
# Is weirdness due to just one or two odd sites?
DHARMa::plotResiduals(simres, mod_dat$yr_sc)
DHARMa::plotResiduals(simres, as.factor(mod_dat$yr_c_f))
DHARMa::plotResiduals(simres, as.factor(mod_dat$odd3))

DHARMa::plotResiduals(simres, mod_dat$perc_forest_100, rank  = FALSE)
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

# If there are diagnostic problems below, not a big deal--unless a lot of problems with clear trends
DHARMa::plotResiduals(simres, mod_dat$perc_forest_50, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$perc_forest_200, rank  = FALSE) # 

DHARMa::plotResiduals(simres, mod_dat$perc_opendev_50, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$perc_opendev_200, rank  = FALSE)

DHARMa::plotResiduals(simres, mod_dat$prop_understory_cover_200, rank  = FALSE) # <<

DHARMa::plotResiduals(simres, mod_dat$understory_cover_sd_50, rank  = FALSE)
DHARMa::plotResiduals(simres, mod_dat$understory_cover_sd_200, rank  = FALSE)

### [IF NEEDED] Append Dharma residuals to model data ----
dat_append <- cbind(mod_dat, predicted = DHARMa::getFitted(mod), resid = DHARMa::getResiduals(mod))
dat_append$scaled_resid = scales::rescale(dat_append$resid, to = c(0,1), from = range(dat_append$resid, na.rm = TRUE))
dat_append$rank_pred = datawizard::ranktransform(dat_append$understory_cover_sd_100)
plotly::ggplotly(ggplot(dat_append, aes(x=understory_cover_sd_100, y = scaled_resid, group = 1, text = paste0(location_name, "<br>Observed: ", sum_indiv, "<br>PredictedResponse: ", round(predicted, 2), "<br>Cov: ", understory_cover_sd_100, "<br>Resid: ", round(resid,2)))) + geom_point() + geom_quantile(), tooltip = "text")


### CHECK THE PREDICTION PLOT ----

## ggemmeans to get year predictions <<<<<<<<<<<<<<<<<<<<<<< CHANGE THIS IN FINAL SUMMARIES
# ggeffects::predict_response(mod, terms = c("yr_c_f", "subunit"), margin = "mean_reference", type = "zero_inflated")) # > Note that to get emmeans, 'margin = "marginalmeans"' but then can't do zero-inflation with that
# (temp_tmb_pred_yearmod <- ggemmeans(mod_tmb_yearmod, terms=c("yr_c_f"), type = "fixed") %>% # # type = ifelse(tmb_zi==TRUE, "fe.zi", "fixed"))) #### ZERO-INFLATION PART DOESN' WORK!!!

year_mod_zi <- ifelse(year_mod$modelInfo$allForm$ziformula=="~1", TRUE, FALSE)

# tmb_pred_yearmod <- ggpredict(year_mod, terms=c("yr_c_f"), type = ifelse(year_mod_zi==TRUE, "zero_inflated", "fixed")) # This is same as below

# # CHeck interaction effects!! [IF APPLICABLE]
# tmb_pred_yearmod_interact <- ggeffects::predict_response(year_mod, terms=c("perc_opendev_50_sc", "odd3"), margin = "mean_mode", type = ifelse(year_mod_zi==TRUE, "zero_inflated", "fixed")) %>% plot()
# saveRDS(tmb_pred_yearmod_interact, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_100m_yearmod_interaction_plot.RDS")))

# Check year estimates
(brms_pred_yearmod <- ggeffects::predict_response(yearmod_brms, terms=c("yr_c_f [all]", "odd3"), margin = "mean_mode", type = ifelse(year_mod_zi==TRUE, "zero_inflated", "fixed")))

(brms_pred_yearmod %<>% 
   as.data.frame() %>%
   janitor::clean_names() %>%
    dplyr::mutate(yr_c = as.numeric(as.character(x))) %>%
    dplyr::select(
      group,
      yr_c,
      yr_predicted = predicted,
      yr_conf_low = conf_low,
      yr_conf_high = conf_high))
brms_pred_yearmod$yr <- brms_pred_yearmod$yr_c+(range(dat$yr) %>% median())

ggplot(data = brms_pred_yearmod) + 
  geom_errorbar(
    mapping = aes(x = yr, y = yr_predicted, ymin = yr_conf_low, ymax = yr_conf_high, color = group), width = 0) +
  geom_point(
    mapping = aes(x = yr, y = yr_predicted, color = group), size = 2.5) +
  facet_grid(rows = vars(group))


### BRMS YEAR MODEL----

mod_dat <- dat

reg_priors <- c(set_prior("normal(0, 1)", class = "b"),
                set_prior("normal(-1.5, 1.5)", class = "b", coef = "odd3odd3_locs")) # this sets odd3 counts to very close to zero and is NOT appropriate when odd3 counts are HIGHER than other locs

yearmod_brms <- brms::brm(
  formula = sum_indiv ~ yr_c_f + (1 | location_name) + hrs_since_rise_sc + perc_forest_100_sc + odd3,
  data = mod_dat,
  # family = brmsfamily("poisson"),
  family = brmsfamily("com_poisson"),
  prior = reg_priors,
  # prior = c(prior(normal(0, 1), class = b)),
  # sample_prior = "only", # <<<<<<<<<<< CHECK PRIORS
  iter = 4000, warmup = 1000, chains = 4, cores = 8)

summary(yearmod_brms)

pp_check(yearmod_brms, ndraws = 500)
pp_check(yearmod_brms, type = "bars_grouped", group = "odd3", ndraws = 200)


yearmod_brms <- add_criterion(yearmod_brms, criterion = c("loo"))

saveRDS(yearmod_brms, here::here("brms_mods_FINAL", paste0("brms_", park, "_", spec, "_100m_yearmod.RDS"))) # <<<<<

# # If there are INFLUENTIAL POINTS, they are usually the ones that are particularly high or low in early or late years
(yearmod_loo <- loo(yearmod_brms)) %>% plot(., label_points = TRUE)

saveRDS(yearmod_loo, here::here("brms_mods_FINAL", paste0("brms_", park, "_", spec, "_100m_yearmod_loo.RDS"))) # <<<<

dat_loo <- cbind(dat[, c("location_name", "event_date", "sum_indiv")], yearmod_loo$diagnostics$pareto_k)
View(dat_loo)


### GLMMTMB TREND MODELS ----
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
### BASE TREND MODEL & CHECK RE VARIANCES
# BASE model is trend with year trend, year RE and location intercept&slope RE 
trend_mod_base_compois <- glmmTMB(sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc | location_name), data = mod_dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8)) 
summary(trend_mod_base_compois)

trend_mod_base_poisson <- glmmTMB(sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc | location_name), data = mod_dat, family = "poisson")
summary(trend_mod_base_poisson)

trend_mod_base_pois_zi <- glmmTMB(sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc | location_name), data = mod_dat, family = "poisson", ziformula = ~1)
summary(trend_mod_base_pois_zi)


## [ONLY IF MODELS WON'T CONVERGE] TRY DROPPING SOME RANDOM EFFECTS ----# IF RE variance is practically zero, consider dropping it
trend_mod_base2 <- glmmTMB(sum_indiv ~ yr_sc + (1 + yr_sc | location_name), data = mod_dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8)) 
summary(trend_mod_base2)


# IF the RE's are highly correlated AND if location slope RE is nearly zero, drop the slope
# plot the relationship btwn intercept and slope RE

trend_mod_base3 <- glmmTMB(sum_indiv ~ yr_sc + (1|yr_c_f) + (1|location_name), data = mod_dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8)) 
summary(trend_mod_base3)

# If the RE's are highly correlated BUT LOCATION SLOPE VARIANCE IS NOT NEARLY ZERO drop the correlation
trend_mod_base4 <- glmmTMB(formula = sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc ||location_name), data = mod_dat, family = "compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))
summary(trend_mod_base4)

trend_mod_base_pois_zi_wind <- glmmTMB(sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc | location_name), data = mod_dat, family = "poisson", ziformula = ~weather_wind_comb)
summary(trend_mod_base_pois_zi_wind)

anova(trend_mod_base4, trend_mod_base3)
# Compare base models
# Exponentiate to get % of population in next year. exp(trend_coef)^(1:10) to see % over next 10 years. Need to bootstrap CI to get correct CI on trend.
trend_mod_base_list <- tibble::lst(trend_mod_base_poisson, trend_mod_base_compois)
performance::compare_performance(trend_mod_base_list, rank = TRUE)
ggstats::ggcoef_compare(tibble::lst(trend_mod0, trend_mod_base4))

(cv_trend_base_mod <- cv::cv(models(trend_mod0=trend_mod0, trend_mod0_Xhab=trend_mod0_Xhab), data = mod_dat, reps = 5))

trend_mod_base_RE <- trend_mod_base3 # <<<<<<<<<<<<< ENTER HERE, BIC-SUPPORTED MODEL

## ADD COVARIATES FROM BEST-FIT YEARMOD 
trend_mod0 <- update(trend_mod_base_RE, formula = ~. + hrs_since_rise_sc + understory_cover_sd_200_sc + odd3)
summary(trend_mod0)

trend_mod1 <- update(trend_mod_base_RE, formula = ~. +   hrs_since_rise_sc + understory_cover_sd_200_sc + prop_understory_cover_200_sc)
summary(trend_mod1)


anova(trend_mod0, trend_mod1)
trend_mod_list <- tibble::lst(trend_mod_base_compois, trend_mod_base3, trend_mod_base4, trend_mod0, trend_mod1, trend_mod2, trend_mod3, trend_mod4)
performance::compare_performance(trend_mod_list, rank = TRUE)
ggstats::ggcoef_compare(tibble::lst(trend_mod0, trend_mod1))


ggstats::ggcoef_table(trend_mod0, intercept = FALSE, table_witdhs = c(2, 1)) + labs(title = "XXX") # make sure no crazy estimates or CI's


## [AS NEEDED] HABITAT INTERACTION WITH YEAR-SLOPE
trend_mod_base_hab50 <- update(trend_mod0, formula = ~.+ yr_sc*hab_type_50_comb)
summary(trend_mod_base_hab50)

trend_mod_base_hab100 <- update(trend_mod0, formula = ~.+ yr_sc*hab_type_100_comb)
summary(trend_mod_base_hab100) # p = 0.043

trend_mod_base_hab200 <- update(trend_mod_base_RE, formula = ~.+ yr_sc*hab_type_200_comb)
summary(trend_mod_base_hab200)


### [IF NEEDED] EXPLORE RANDOM EFFECTS
# Plot correlation btwn location intercept RE & slope RE
trend_mod <- trend_mod0

re <-ranef(trend_mod)[1]$cond$location_name %>%
  janitor::clean_names()
re
ggplot(re, aes(x = intercept, y = yr_sc)) + geom_point()
ggplot(re, aes(intercept)) + geom_histogram()
ggplot(re, aes(yr_sc)) + geom_histogram()

### ASSIGN TREND_MOD 
trend_mod <- trend_mod2_prior
mod <- trend_mod

### SAVE TRENDMOD RESULTS 
saveRDS(trend_mod_list, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_100m_trendmod_modlist.RDS")))
saveRDS(trend_mod, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_100m_trendmod.RDS")))
saveRDS(mod_dat, here::here("tmb_mods_FINAL", paste0("tmb_", park, "_", spec, "_100m_trendmod_dat.RDS")))

### BRMS TREND MODEL----
# Brms can be useful for identifying highly influential points and seeing if those points are responsible for model diagnostic issues.
mod_dat <- dat
           
trendmod_brms <- brms::brm(
  formula =  sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc|location_name) + hrs_since_rise_sc + perc_forest_100_sc + odd3,
  # formula =  sum_indiv ~ yr_sc + (1 | yr_c_f) + (1 + yr_sc ||location_name) + hrs_since_rise_sc + understory_cover_sd_200_sc,
  # formula =  sum_indiv ~ yr_sc + (1 | yr_c_f) + (1|location_name) + hrs_since_rise_sc + understory_cover_sd_200_sc + weather_wind_comb,
  # family = brmsfamily("poisson"),
  family = brmsfamily("com_poisson"),
  prior = reg_priors,
  # sample_prior = "only", # <<<<<<<<<<< CHECK PRIORS
  data = mod_dat, 
  iter = 4000, warmup = 1000, chains = 4, cores = 8)
pp_check(trendmod_brms, ndraws = 500)
trendmod_brms <- add_criterion(trendmod_brms, criterion = c("loo"))

saveRDS(trendmod_brms, here::here("brms_mods_FINAL", paste0("brms_", park, "_", spec, "_100m_trendmod.RDS"))) # <<<<<

# # If there are INFLUENTIAL POINTS, they are usually the ones that are particularly high or low in early or late years
(trendmod_loo <- loo(trendmod_brms)) %>% plot(., label_points = TRUE)

saveRDS(trendmod_loo, here::here("brms_mods_FINAL", paste0("brms_", park, "_", spec, "_100m_trendmod_loo.RDS"))) # <<<<

dat_loo <- cbind(dat[, c("location_name", "event_date", "sum_indiv")], mod_loo$diagnostics$pareto_k)
View(dat_loo)




### EVALUATE BRMS models (TRENDMOD_BRMS OR MOD_BRMS) ----

mod_brms = trendmod_brms # <<<<<<<<<
summary(mod_brms) # slightly shorter summary of model results

rstan::check_hmc_diagnostics(mod_brms$fit)
bayes_R2(mod_brms)

# Coefficient plot
mcmc_plot(mod_brms) # quick coefficient plot

# Model convergence check
plot(mod_brms) # plots showing posterior distribution of parameters, and trace plots

# Posterior predictive checks
pp_check(mod_brms, type = "stat_2d", ndraws = 200)
# pp_check(mod_brms, type = "loo_pit_qq", ndraws = 200)
pp_check(mod_brms, type = "rootogram", ndraws = 200)
pp_check(mod_brms, type = "bars", ndraws = 200)


# pp_check(mod_brms, type = "scatter_avg_grouped", group = "yr_c_f", ndraws = 200) + 
  # geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
pp_check(mod_brms, type = "bars_grouped", group = "yr_c_f", ndraws = 200)
pp_check(mod_brms, type = "bars_grouped", group = "location_name", ndraws = 200)
pp_check(mod_brms, type = "bars_grouped", group = "odd3", ndraws = 200)

# If categorical predictors...
pp_check(mod_brms, type = "bars_grouped", group = "weather_wind_comb", ndraws = 200) # <<<<<<<< CHANGE, DEPENDING ON BEST MODEL

pp_check(mod_brms, type = "bars_grouped", group = "odd3", ndraws = 200) # <<<<<<<< CHANGE, DEPENDING ON BEST MODEL

## Plot conditional effects
# plot(conditional_effects(mod_brms, effects ="understory_cover_sd_200_sc:prop_understory_cover_50_sc", points=TRUE, point_args = list(width = 0.3, height = 0.1, alpha = 0.4)))

# Save conditional effects plots as list, for use with ggplot functions
(p <- plot(conditional_effects(mod_brms), points=TRUE, point_args = list(width = 0.3, height = 0.1, alpha = 0.4), plot = FALSE))

# saveRDS(p, here::here("brms_mods_FINAL", paste0("brms_", park, "_", spec, "_100m_condeffects.RDS"))) # <<<
saveRDS(p, here::here("brms_mods_FINAL", paste0("brms_", park, "_", spec, "_100m_trendmod_condeffects.RDS"))) # <<<

## DHARMA
simres <- DHARMa.helpers::dh_check_brms(mod_brms)
## NOW GO TO DHARMA CODE

### PLOT PREDICTIONS OF FINAL TRENDS MODEL ----
# Check the estimated trend
# em=emmeans::emtrends(trendmod_brms, specs = c("yr_sc", "odd3"), var = c("yr_sc"))
em=emmeans::emtrends(trendmod_brms, specs = c("yr_sc"), var = "yr_sc") # trends at mean year. Ignores a 'type' argument, it seems
summary(em, infer=c(TRUE,TRUE),null=0) %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  dplyr::mutate(signif = sign(lower_hpd) == sign(upper_hpd))


### PLOT TREND WITH YEAR ESTIMATES ----
tmb_zi_trendmod <- ifelse(trend_mod$modelInfo$allForm$ziformula=="~1", TRUE, FALSE)

# (temp_trendmod_mmeans <- ggeffects::predict_response(trend_mod, terms=c("yr_sc", "hab_type_100_comb"), margin = "marginalmeans", type = ifelse(tmb_zi_trendmod==TRUE, "zero_inflated", "fixed"))) # this is the same as emmeans. For factors, it weights them equally regardless of how many are actually in each factor level
# (temp_trendmod_mmeans <- ggeffects::predict_response(trend_mod, terms=c("yr_sc"), margin = "marginalmeans", type = ifelse(tmb_zi_trendmod==TRUE, "zero_inflated", "fixed")))

# (temp_trendmod_mmode <- ggeffects::predict_response(trend_mod, terms=c("yr_sc", "hab_type_100_comb"), margin = "mean_mode", type = ifelse(tmb_zi_trendmod==TRUE, "zero_inflated", "fixed")))    
(temp_trendmod_mmode <- ggeffects::predict_response(trendmod_brms, terms=c("yr_sc", "odd3"), margin = "mean_mode", type = ifelse(tmb_zi_trendmod==TRUE, "zero_inflated", "fixed")))  

temp_trendmod_mmode %<>% 
  as.data.frame 

FuncFormatTrendDat <- function(pred, dat) {
  
  yr_df <- dat %>% dplyr::select(yr, yr_sc) %>% distinct(.)
  
  pred_out <- pred %>% 
    as.data.frame() %>%  
    janitor::clean_names() %>%
    dplyr::left_join(yr_df, by = c("x" = "yr_sc"))
  
  if("group" %in% names(pred_out)) {
    
    pred_out %<>%
      dplyr::select(
        group,
        yr,
        trend_predicted = predicted,
        trend_conf_low = conf_low,
        trend_conf_high = conf_high)
  } else {
    pred_out %<>%
      dplyr::select(
        yr,
        trend_predicted = predicted,
        trend_conf_low = conf_low,
        trend_conf_high = conf_high)
  }
  
  return(pred_out)
}

# tmb_pred_trendmod_mmeans <- FuncFormatTrendDat(pred = temp_trendmod_mmeans, dat)
brms_pred_trendmod_mmode <- FuncFormatTrendDat(pred = temp_trendmod_mmode, dat) # <<<<<<<<<<<<< FIX--THIS ISN'T QUITE RIGHT AS A COMPARISON TO YEARMOD

(p <- ggplot() + 
    geom_errorbar(data = brms_pred_yearmod, aes(x = yr, y = yr_predicted, ymin = yr_conf_low, ymax = yr_conf_high), width = 0, color = "blue") +
    geom_point(data = brms_pred_yearmod, aes(x = yr, y = yr_predicted), size = 2.5, color = "blue") +
    # geom_ribbon( # glmmTMB TREND RESULTS AS RIBBON
    #   data = tmb_pred_trendmod_mmeans, aes(x = yr, y = trend_predicted, ymin = trend_conf_low, ymax = trend_conf_high),
    #   alpha = 0.1, fill = "orange") +
    # geom_line(data =tmb_pred_trendmod_mmeans, aes(x = yr, y= trend_predicted), color = "orange", linetype = "dashed") +
    geom_ribbon(
      data = brms_pred_trendmod_mmode, aes(x = yr, y = trend_predicted, ymin = trend_conf_low, ymax = trend_conf_high),
      alpha = 0.1, fill = "orange") +
    geom_line(data =brms_pred_trendmod_mmode, aes(x = yr, y= trend_predicted), color = "orange", linetype = "dashed") +
    scale_x_continuous(breaks = unique(sort(dat$yr))) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title = paste0("Trends in relative abundance of ", spec, " at ", park), subtitle = "(predictors held at mean values)", x = "Year", y = "Index of abundance (# of birds/point count)") +
    theme_bw() +
    facet_grid(rows = vars(group)))

saveRDS(p, here::here("brms_mods_FINAL", paste0("tmb_", park, "_", spec, "_100m_trendplot.RDS")))
