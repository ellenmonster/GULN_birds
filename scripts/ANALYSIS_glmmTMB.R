### GULN BIRDS - glmmTMB

### FINAL MODELS ----
# > Should use % forest w/in 100m as covariate, instead of the two-level habitat type. But remember that different parks may define "forest" differently
# > brms cross-validate for "model selection" (instead of AICc)
# > ORDINATION to identify points that are very "different" from expectations based on habitat? Remember there were a couple sites at one park that I was going to suggest dropping b/c numbers often inexplicably very different there
# > With year-factor models, glmmTMB RE variances seem way too small, and get those strangely huge 95%CI's for years with 0 detections (for Poisson model; then for ZIP, the estimates are NaN for 95% CI's)
# > With year-factor models, rank-deficiency when researcher changes by year (e.g., GUIS had Mark Woodrey 2013-2016, then Adam Rohnke 2017-2022)... why isn't it parsing out that researcher effect?
# > QUESTION... For zero-inflated models, can/should the zero inflation be a function of location (which is a RE)?
# > For each park, write out the pre-analysis data-cleaning rules such as dropping certain researcher, any "odd" locations that should be excluded (perhaps for a subset of sites)--get QE thoughts about when it's okay to just drop a researcher who only surveyed for a year or two of the early or late years
# > For BITH, TC29 is surveyed twice in 2018, so there is just one record that is 2018_2. In final, decide if this should be excluded from analyses?
# > Some glmmTMB RE variance estimates are near zero--these may be wrong! May need to try with a different package

### TRY THESE ----
# Just try without the other cov
# See if this is also happening on the Poisson and NB or only with zero-inflated
# With n-mixture (detection estimates) and with trend model should not have this problem, right?
### NOTES ----
# > VICK-- only use data from Daniel and only starting in 2012
# > PAAL: Dropped 3 records b/c missing weather_wind_comb. Also drop hab_type_100 b/c too many missing. Has multiple researchers, so add as predictor--BUT REALLY should consider dropping the first two years so can just have one researcher [for the low pop 100m models I did drop the first two years--so not comparable to the higher pop models]. PAAL only has non-forest habitat.
# > BITH: Has multiple researchers, so add as predictor AND REMOVE FIRST THREE YEARS B/C RESEARCHER CHANGES CORRELATED WITH YEARS. BITH only has forest habitat.
# > SAAN: Dropped 1 record b/c missing hrs_since_rise. Add subunits: SAAN1-6 ARE 'RA', AND REST ARE 'MI' SUBUNIT and yr*subunit interaction. SAAN-RA has no surveys in 2013??
# > GUIS: Add subunits and yr*subunit interaction. GUIS =-MS only has forest habitat, so rank deficient in a model with subunit and forest. Something weird with GUIS-FL b/c from 2013-2017, alternate panels and the understory sd values are very different for panel 1 (GUISFL16 & 27 lowest) vs panel 2 (GUISFL 28 & 29 very high)
# > JELA: Only has one researcher

### NOTES ----
# > For N-mixture and other models I dropped 3 VICK surveys from 2010 to keep a consistent max 2 surveys per year. I will do the same here just for comparison across approaches.
# > For species that don't have habitat as important predictor, it seems I can fit models with all parks combined b/c the habitat predictor is the only weird one that would have different meanings at different parks.
# > For purposes of trend estimate, it doesn't matter if Dharma on habitat and sd are bad--those are not going to affect trend estimate

### GENERAL STEPS ----
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
# 5. For the single-visit models, start with the final model from using all data, then use the model diagnostics to see if there are any apparent issues that warrant trying a different model


rm(list = ls())

### Load packages ----

pkgs <- c("tidyverse", "partykit", "magrittr", "DHARMa", "emmeans", "ggeffects", "tidybayes", "broom", "broom.mixed", "glmmTMB", "ggstats", "performance", "AICcmodavg", "here")
installed_pkgs <- pkgs %in% installed.packages()
if (length(pkgs[!installed_pkgs]) > 0) install.packages(pkgs[!installed_pkgs], repos = "https://cloud.r-project.org" , dep=TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

source(here::here("scripts", "ANALYSIS_STARTUP.R"))

# ### RUN SINGLE-VISIT YEARMOD TMB ----           
# park = "JELA"
# spec = "CACH"
# 
# template_tmb <- readRDS(here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_yearmod.RDS")))
# 
# # >>>> RUN DATA FORMATTING CODE TO GET DAT <<<<<<<<<<<<<<
# View(dat)
# nrow(dat) == nrow(template_tmb$frame) # check this is TRUE
# 
# dat$visit <- as.integer(sub(".*_", "", dat$yr_visit))
# # View(dat %>% dplyr::select(yr_visit, visit)) # check it
# dat %<>% dplyr::filter(visit == 1)
# View(dat)
# 
# template_tmb$call # >>> THEN PASTE BELOW
# 
# mod <- glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name) + hab_type_200_comb + understory_cover_sd_200_sc + julian_prop_c + hrs_since_rise_c + I(hrs_since_rise_c^2), data = dat, family = compois, ziformula = ~0, dispformula = ~1)
# summary(mod)
# summary(template_tmb) # These CI's should be smaller
# 
# saveRDS(mod, here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_100m_1visit_yearmod.RDS")))

### Read in formatted data ----
# Subset and scale covariates ----

park = "GUIS"
spec = "GCFL"

dat <- FuncFormatFullObs(park = park, spec = spec, limit_100m = TRUE) # <<<<<<<< 'limit_100m = TRUE' to only use data within 100m of observer
# dat %<>% dplyr::filter(subunit == "GUIS-FL"); dat$subunit <- as.character(dat$subunit)
# dat %<>% dplyr::filter(subunit == "MI"); dat$subunit <- as.character(dat$subunit)


mean(dat$sum_indiv, na.rm = TRUE)
sum(dat$sum_indiv)

table(dat$researcher, dat$yr, dat$subunit) # check how many observers

# Dataset with complete cases only
# PAAL: Dropped 3 records b/c missing weather_wind_comb
# SAAN: Dropped 1 record b/c missing hrs_since_rise
nrow(dat)
colSums(is.na(dat))

dat <- dat[complete.cases(dat),] # <<<<<<<<<< SHOULD DO THIS AFTER DETERMINING COVARIATES?

# Scale predictors
dat$prop_understory_cover_50_sc <- scale(dat$prop_understory_cover_50) %>% as.vector()
dat$prop_understory_cover_100_sc <- scale(dat$prop_understory_cover_100) %>% as.vector()
dat$prop_understory_cover_200_sc <- scale(dat$prop_understory_cover_200) %>% as.vector()
dat$understory_cover_sd_50_sc <- scale(dat$understory_cover_sd_50) %>% as.vector()
dat$understory_cover_sd_100_sc <- scale(dat$understory_cover_sd_100) %>% as.vector()
dat$understory_cover_sd_200_sc <- scale(dat$understory_cover_sd_200) %>% as.vector()


nrow(dat)
table(dat$yr)
table(dat$yr, dat$location_name)
dat %>% dplyr::group_by(subunit) %>% summarize(tot = sum(sum_indiv))

### STOP HERE FOR JUST BUILDING OFF PRIOR MODELS

### Plots for exploring covariates ----

## SAAN & GUIS--add subunit
## as needed... researcher + first_yr (BITH, PAAL, GUIS?)
names(dat)

dat_tree <- partykit::ctree(sum_indiv ~ subunit, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~ researcher, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~ hab_type_200_comb +hab_type_100_comb, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~prop_understory_cover_50 + prop_understory_cover_100 + prop_understory_cover_200 + understory_cover_sd_50 + understory_cover_sd_100 + understory_cover_sd_200, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~prop_understory_cover_50 + prop_understory_cover_100 + prop_understory_cover_200, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~understory_cover_sd_100+understory_cover_sd_200, data = dat)
plot(dat_tree)

dat_tree <- partykit::ctree(sum_indiv ~weather_wind_comb + julian_prop + hrs_since_rise, data = dat)
plot(dat_tree)

# dat %<>% dplyr::filter(subunit=="MI") # <<<<< SAAN ONLY, AS NEEDED

# Check data patterns with ggplot()
ggplot(dat, aes(x=sum_indiv)) + geom_histogram() + facet_wrap(~subunit)

ggplot(dat, aes(y = sum_indiv, x= yr_c_f)) + geom_boxplot() + facet_wrap(~subunit)
dat %>% dplyr::group_by(yr) %>% dplyr::summarise(tot = sum(sum_indiv))
table(dat$hab_type_100_comb, dat$sum_indiv, dat$subunit)

# CHANGE BTWN 100 AND 200 <<<<<<
ggplot(dat, aes(y = sum_indiv, x= hab_type_200_comb)) + geom_boxplot() + facet_wrap(~subunit)
table(dat$hab_type_100_comb, dat$sum_indiv)

# RESEARCHER & FIRST YR
ggplot(dat, aes(y = sum_indiv, x= as.factor(researcher), fill = as.factor(first_yr))) + geom_boxplot()+ facet_wrap(~subunit)

# WIND
ggplot(dat, aes(y = sum_indiv, x= weather_wind_comb)) + geom_boxplot()+ facet_wrap(~subunit)
table(dat$weather_wind_comb)

# UNDERSTORY COVER
table(dat$sum_indiv) # check which sum_indiv categories have large enough samples to be reliable
  
# understory cover
ggplot(dat, aes(y = prop_understory_cover_50, x = as.factor(sum_indiv), fill = hab_type_100_comb)) + geom_boxplot()+ facet_wrap(~subunit)
ggplot(dat, aes(y = prop_understory_cover_100, x = as.factor(sum_indiv), fill = hab_type_100_comb)) + geom_boxplot()+ facet_wrap(~subunit)
ggplot(dat, aes(y = prop_understory_cover_200, x = as.factor(sum_indiv), fill = hab_type_200_comb)) + geom_boxplot()+ facet_wrap(~subunit)

# understory cover SD
ggplot(dat, aes(y = understory_cover_sd_50, x = as.factor(sum_indiv), fill = hab_type_100_comb)) + geom_boxplot()+ facet_wrap(~subunit)
ggplot(dat, aes(y = understory_cover_sd_100, x = as.factor(sum_indiv), fill = hab_type_100_comb)) + geom_boxplot()+ facet_wrap(~subunit)
ggplot(dat, aes(y = understory_cover_sd_200, x = as.factor(sum_indiv), fill = hab_type_100_comb)) + geom_boxplot()+ facet_wrap(~subunit)

# FOR POSSIBLE INTERACTION - UNDERSTORY COVER, COLORED BY HABITAT TYPE. CHANGE BTWN 50, 100, 200
ggplot(dat, aes(y = understory_cover_sd_100, x= hab_type_100_comb)) + geom_boxplot()+ facet_wrap(~subunit)



ggplot(dat, aes(x = understory_cover_sd_100, y = sum_indiv, color = hab_type_100_comb)) + geom_smooth() + geom_jitter()+ facet_wrap(~subunit)

# JULIAN AND TIME
ggplot(dat, aes(x = julian_prop_c, y = sum_indiv)) + geom_jitter() + geom_smooth() + facet_wrap(~subunit)
ggplot(dat, aes(x = hrs_since_rise_c, y = sum_indiv)) + geom_jitter() + geom_smooth() + facet_wrap(~subunit)


### MODEL NOTES ----
## NOTES ON BITH MODEL
# > exclude habitat at 200m scale bc only one level
# > incl researcher & first_yr b/c many--researcher is confounded with year. So this one is a bit funny b/c a lot of the plots seem to suggest much higher counts in the early years, but the model is attributing those to researcher differences, not to actual high numbers... so BBS models, e.g., say NOCA is declining everywhere but BITH NOCA says nonsign. increase. Will be interesting to see NOCA trend at parks with less researcher turnover
## PAAL--exclude hab b/c only non-forest. Has two researchers.
## JELA--only one researcher

### YEAR MODELS ----
### RUN 2-VISIT, 100M LIMIT MODELS ----
# Do this after running the data formatting with 100m limit...
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

# Others?
mod0pois_zi <-glmmTMB(formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), data = dat, family = poisson, ziformula = ~1)

# AICc to compare
mod_list <- tibble::lst(mod0pois, mod0comp, mod0nb, mod0pois_zi)
performance::compare_performance(mod_list, rank = TRUE)
# aictab(mod_list, second.ord = TRUE)
ggstats::ggcoef_compare(mod_list) 

# Assign best one to mod0 (if underdispersed, it may be comp instead of pois regardless of AICc?)
mod0 <- mod0pois_zi


# Add covariates--FOR SOME REASON, SCALED COVARIATES SOMETIMES HAVE PROBLEMS BUT UNSCALED DO NOT

mod1 <- update(mod0, formula = ~. + prop_understory_cover_50_sc)
summary(mod1)

mod2 <- update(mod0, formula = ~. + hab_type_100_comb)
summary(mod2)

mod3 <- update(mod0, formula = ~. + hab_type_100_comb * prop_understory_cover_50_sc)
summary(mod3)

mod4 <- update(mod0, formula = ~. + prop_understory_cover_200_sc + I(prop_understory_cover_200_sc^2))
summary(mod4)
               
               
               
               
mod5 <- update(mod0, formula = ~. + julian_prop_c + I(julian_prop_c^2))

mod6 <- update(mod0, formula = ~.+ hab_type_100_comb)

mod7 <- update(mod0, formula = ~.+ hab_type_100_comb + understory_cover_sd_50_sc)

mod8 <- update(mod0, formula = ~.+ hab_type_100_comb + understory_cover_sd_50_sc + julian_prop_c)


mod4 <- update(mod0, formula = ~. + subunit + hab_type_100_comb)

mod5 <- update(mod0, formula = ~. + subunit + prop_understory_cover_100_sc)

mod3_nozi <- update(mod3, ziformula = ~0)




### Compare model estimates ----
# Use results to see if other combinations of covariates need to be considered, and revise model list as necessary

mod_list = tibble::lst(mod0, mod1, mod2, mod3)
performance::compare_performance(mod_list, rank = TRUE)
# aictab(mod_list, second.ord = TRUE) # use DHARMa outputs over this, to decide on final model--but confirm that best AICc-supported model and best DHARMa-supported model have nearly identical estimates, esp. for trend. second.ord means AICc instead of AIC
ggstats::ggcoef_compare(mod_list) # do this for the models within delta-2 of best-supported

# Save the modlist with ALL models
saveRDS(mod_list, here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_100m_yearmod_modlist.RDS")))

### Diagnostics on top model ----
mod <- mod2


### Check the model ----
summary(mod)
# ggstats::ggcoef_table(mod, intercept = FALSE, table_witdhs = c(2, 1)) + labs(title = "XXX")

## pp-checks
# performance::check_model(mod1)
# performance::check_predictions(mod, check_range = TRUE)
performance::check_collinearity(mod)
performance::check_autocorrelation(mod) # This may be a problem any time the model includes an interaction???
# performance::check_distribution(mod) # NOTE: This function may help to check a models' distributional family and see if the model-family probably should be reconsidered. Since it is difficult to exactly predict the correct model family, consider this function as somewhat experimental.

### Dharma checks for glmmTMB ----

simres <- DHARMa::simulateResiduals(fittedModel = mod, n = 1000, refit = FALSE, plot = TRUE) # re-simulate all levels, to test model structure as a whole. The default is resimulating on all levels (unconditional sims) instead of simulating conditional on the fitted random effects.
simres <- DHARMa::simulateResiduals(fittedModel = mod, n = 1000, refit = FALSE, re.form = NULL, plot = TRUE) #re-simulate conditional on the fitted random effects, i.e., the RE's are not re-simulated. More sensitive for detecting problems, but it doesn't test all possible options.

DHARMa::testZeroInflation(simres)
DHARMa::testDispersion(simres) # test for over/under dispersion
DHARMa::testOutliers(simres, type = "bootstrap")

# Plot residuals vs each independent factor
DHARMa::plotResiduals(simres, as.factor(dat$yr_c_f))
DHARMa::plotResiduals(simres, dat$hab_type_100_comb, rank  = FALSE)
DHARMa::plotResiduals(simres, dat$hab_type_200_comb, rank  = FALSE)
DHARMa::plotResiduals(simres, dat$prop_understory_cover_200, rank  = FALSE) # 
DHARMa::plotResiduals(simres, dat$understory_cover_sd_200, rank  = FALSE)
DHARMa::plotResiduals(simres, dat$prop_understory_cover_100, rank  = FALSE) # 
DHARMa::plotResiduals(simres, dat$understory_cover_sd_100, rank  = FALSE)
DHARMa::plotResiduals(simres, dat$prop_understory_cover_50, rank  = FALSE) # 
DHARMa::plotResiduals(simres, dat$understory_cover_sd_50, rank  = FALSE)
DHARMa::plotResiduals(simres, dat$weather_wind_comb, rank  = FALSE)
DHARMa::plotResiduals(simres, dat$julian_prop_c, rank  = FALSE)
DHARMa::plotResiduals(simres, dat$hrs_since_rise_c, rank  = FALSE)
DHARMa::plotResiduals(simres, dat$researcher, rank  = FALSE)
DHARMa::plotResiduals(simres, dat$first_yr, rank  = FALSE)


### SAVE THE FINAL 100M 2-VISIT YEAR MODEL ----
saveRDS(mod, here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_100m_yearmod.RDS")))

### CHECK THE PLOT ----

## ggemmeans to get year predictions <<<<<<<<<<<<<<<<<<<<<<< CHANGE THIS IN FINAL SUMMARIES
# ggeffects::predict_response(mod, terms = c("yr_c_f", "subunit"), margin = "mean_reference", type = "zero_inflated")) # > Note that to get emmeans, 'margin = "marginalmeans"' but then can't do zero-inflation with that
(tmb_pred_yearmod <- ggeffects::predict_response(mod, terms=c("yr_c_f"), type = "zero_inflated") %>% # type = "fixed" or "zero_inflated"
   as.data.frame() %>%
   janitor::clean_names() %>%
  dplyr::mutate(num_visit = 2))

ggplot(data = tmb_pred_yearmod) + 
  geom_errorbar(
    mapping = aes(x = x, y = predicted, ymin = conf_low, ymax = conf_high), width = 0, color = "blue") +
  geom_point(
    mapping = aes(x = x, y = predicted), size = 2.5, color = "blue") +
  facet_grid(rows = vars(group))


### NOW REDO WITH ONLY ONE VISIT PER YEAR, AND USE THE SAME MODEL (BUT CHECK THE DIAGNOSTICS) ----
dat1visit <- dat[grep("_1", dat$yr_visit), ]
sort(unique(dat1visit$yr_visit)) # check all sample years included, but only first visit each year
table(dat1visit$sum_indiv)

mean(dat1visit$sum_indiv, na.rm = TRUE)
sum(dat1visit$sum_indiv)

mod1visit <- update(mod, data = dat1visit)

# # If need to try alternative models for 1visit...
# mod1visit <-update(mod0pois, data = dat1visit)
# 
# # other models?
m2 <- update(mod1visit, formula = ~. + hrs_since_rise_c)
summary(m2)

m3 <- update(mod1visit, formula = ~. + hrs_since_rise_c)
summary(m3)

m4 <- update(mod1visit, formula = ~. +understory_cover_sd_100 + hrs_since_rise_c)
summary(m4)

mod3 <- update(mod2, ziformula = ~1)
mod4 <- update(mod1visit, ziformula = ~1)
mod5 <- update(mod1visit, family="compois", ziformula = ~0, dispformula = ~researcher, control = glmmTMBControl(parallel = 8))
mod6 <- update(mod1visit, formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), family="compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))
mod6 <- update(mod1visit, formula = sum_indiv ~ 0 + yr_c_f + (1 | location_name), family=nbinom2)
# 
# mod3 <- update(mod1visit, formula = ~. + prop_understory_cover_100_sc)
# 
# mod4 <- update(mod1visit, family="compois", ziformula = ~0, dispformula = ~1, control = glmmTMBControl(parallel = 8))
# 
# mod5 <- update(mod1visit, ziformula = ~1)
# 
# mod4 <- update(mod1visit, formula = ~. + subunit + hab_type_100_comb * prop_understory_cover_100_sc)
# 
# mod5 <- update(mod1visit, formula = ~. + subunit + hab_type_100_comb)
# 
# mod6 <- update(mod1visit, formula = ~. + subunit + prop_understory_cover_100_sc)
# 
# mod7 <- update(mod1visit, formula = ~. + subunit * hab_type_100_comb * prop_understory_cover_100_sc)
# 
# mod8 <- update(mod1visit, formula = ~. + subunit + hab_type_100_comb + prop_understory_cover_100_sc)
# 

performance::compare_performance(tibble::lst(mod1visit, m2), rank = TRUE)
# aictab(tibble::lst(mod1visit, m2, m3, m4))
# ggstats::ggcoef_compare(mod_list)
# 
# mod1visit <- m2

### Check the model ----
summary(mod1visit)
# ggstats::ggcoef_table(mod, intercept = FALSE, table_witdhs = c(2, 1)) + labs(title = "XXX")

## pp-checks
# performance::check_model(mod1visit)
# performance::check_predictions(mod1visit, check_range = TRUE)
performance::check_collinearity(mod1visit)
performance::check_autocorrelation(mod1visit)
# performance::check_distribution(mod1visit)
performance::check_overdispersion(mod1visit)

### Dharma checks for glmmTMB ----

simres1visit <- DHARMa::simulateResiduals(fittedModel = mod1visit, n = 1000, refit = FALSE, plot = TRUE) # re-simulate all levels, to test model structure as a whole. The default is resimulating on all levels (unconditional sims) instead of simulating conditional on the fitted random effects.
simres1visit <- DHARMa::simulateResiduals(fittedModel = mod1visit, n = 1000, refit = FALSE, re.form = NULL, plot = TRUE) #re-simulate conditional on the fitted random effects, i.e., the RE's are not re-simulated. More sensitive for detecting problems, but it doesn't test all possible options.

DHARMa::testZeroInflation(simres1visit)
DHARMa::testDispersion(simres1visit) # test for over/under dispersion
DHARMa::testOutliers(simres1visit, type = "bootstrap")

# Plot residuals vs each independent factor
DHARMa::plotResiduals(simres1visit, as.factor(dat1visit$yr_c_f))
DHARMa::plotResiduals(simres1visit, dat1visit$hab_type_100_comb, rank  = FALSE)
DHARMa::plotResiduals(simres1visit, dat1visit$hab_type_200_comb, rank  = FALSE)
DHARMa::plotResiduals(simres1visit, dat1visit$prop_understory_cover_200, rank  = FALSE) # 
DHARMa::plotResiduals(simres1visit, dat1visit$understory_cover_sd_200, rank  = FALSE)
DHARMa::plotResiduals(simres1visit, dat1visit$prop_understory_cover_100, rank  = FALSE) # 
DHARMa::plotResiduals(simres1visit, dat1visit$understory_cover_sd_100, rank  = FALSE)
DHARMa::plotResiduals(simres1visit, dat1visit$prop_understory_cover_50, rank  = FALSE) # 
DHARMa::plotResiduals(simres1visit, dat1visit$understory_cover_sd_50, rank  = FALSE)
DHARMa::plotResiduals(simres1visit, dat1visit$weather_wind_comb, rank  = FALSE)
DHARMa::plotResiduals(simres1visit, dat1visit$julian_prop_c, rank  = FALSE)
DHARMa::plotResiduals(simres1visit, dat1visit$hrs_since_rise_c, rank  = FALSE)
DHARMa::plotResiduals(simres1visit, dat1visit$researcher, rank  = FALSE)
DHARMa::plotResiduals(simres1visit, dat1visit$first_yr, rank  = FALSE)

### CHECK THE PLOT ----

## ggemmeans to get year predictions
(tmb_pred_yearmod1visit <- ggeffects::predict_response(mod1visit, terms=c("yr_c_f"), type = "zero_inflated") %>% # type = "fixed"
   as.data.frame() %>%
   janitor::clean_names() %>%
  dplyr::mutate(num_visit = 1))

plot_dat <- rbind(tmb_pred_yearmod, tmb_pred_yearmod1visit)

ggplot(data = plot_dat, aes(x = x, y = predicted, ymin = conf_low, ymax = conf_high, group = as.factor(num_visit), color = as.factor(num_visit))) + 
  geom_errorbar(width = 0, position = position_dodge(0.3)) +
  geom_point(size = 2.5, position = position_dodge(0.3)) +
  scale_color_manual(values = c("1" = "orange", "2" = "blue"))+
  facet_grid(rows = vars(group))

### SAVE THE FINAL 100M 1-VISIT YEAR MODEL ----
saveRDS(mod1visit, here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_100m_1visit_yearmod.RDS")))


# ### PLOT glmmTMB PREDICTIONS OF FINAL TRENDS MODEL ----
# em=emmeans::emtrends(mod, specs = c("yr_c"), var = "yr_c") # trends at mean year. Ignores a 'type' argument, it seems
# summary(em, infer=c(TRUE,TRUE),null=0) %>%
#   as.data.frame() %>%
#   janitor::clean_names() %>%
#   dplyr::mutate(signif = p_value < 0.05)
# 
# ggpredict(mod, terms=c("yr_c")) # I think ggpredict uses the most common for factor levels, while ggemmeans will do some kind of averaging. Can get separate estimates for each researcher by terms = c("yr_c", "researcher")
# x<-ggemmeans(mod, terms="yr_c")
# 
# x<- ggpredict(mod, terms = list(
#   yr_c = seq(-5, 5, by =1)
#   )
#   )
# 
# plot(x, facets = FALSE)
# 
# 
# ### TREND MODELS ----
# >>>>> CHECK IF I SHOULD HAVE ADDED DISPFORMULA =~1, OR IF THAT IS DEFAULT
# min_tmb <- glmmTMB(sum_indiv ~ yr_c + (1 | yr_c_f) + (1 | location_name) + hab_type_100_comb + julian_prop_c + weather_wind_comb, family = compois, data = dat)
# 
# min2_tmb <- glmmTMB(sum_indiv ~ yr_c + (1 | yr_c_f) + (1 | location_name) + hab_type_200_comb + understory_cover_sd_200 + julian_prop_c + hrs_since_rise_c + I(hrs_since_rise_c^2) + prop_understory_cover_200, family = compois, data = dat)
# 
# ### YEAR MODELS IF ALREADY HAVE TREND MODEL ----
# 
# park = "VICK"
# spec = "WEVI"
# mod <- readRDS(here::here("tmb_mods", paste0("tmb_", park, "_", spec, ".RDS")))
# mod$modelInfo$allForm # check formula
# 
# ### REMOVE RESEARCHER & FIRST YEAR AND YR RE AND YR FE CONTINUOUS
# yr_mod <- glmmTMB(sum_indiv ~ 0 + yr_c_f + (1 | location_name) + hab_type_200_comb + 
#                     understory_cover_sd_200_sc + julian_prop_c + hrs_since_rise_c + I(hrs_since_rise_c^2), family = compois, data = dat)
# 
# summary(yr_mod)
# 
# tmb_zi <- ifelse(yr_mod$modelInfo$allForm$ziformula=="~1", TRUE, FALSE)
# ggpredict(yr_mod, terms=c("yr_c_f"), type = ifelse(tmb_zi==TRUE, "zero_inflated", "fixed"))
# saveRDS(yr_mod, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_yearmod.RDS")))
# ##### <<<<<<<<<<<<<<<<<<<<<<<<
