### GULN BIRDS - BRMS

### PAAL --DROP THE FIRST TWO YEARS THAT ARE A DIFFERENT RESEARCHER

### CHANGE TO MAKE IN DATA FOR OTHER ANALYSES ALSO ----
# > VICK-- only use data from Daniel and only starting in 2012
# > PAAL: Dropped 3 records b/c missing weather_wind_comb. Also drop hab_type_100 b/c too many missing. Has multiple researchers, so add as predictor
# > BITH: Has multiple researchers, so add as predictor AND REMOVE FIRST THREE YEARS B/C RESEARCHER CHANGES CORRELATED WITH YEARS
# > SAAN: Dropped 1 record b/c missing hrs_since_rise. Add subunits: SAAN1-6 ARE 'RA', AND REST ARE 'MI' SUBUNIT and yr*subunit interaction
# > GUIS: Add subunits and yr*subunit interaction

### NOTES ----
# > For N-mixture and other models I dropped 3 VICK surveys from 2010 to keep a consistent max 2 surveys per year. I will do the same here just for comparison across approaches.
# > For species that don't have habitat as important predictor, it seems I can fit models with all parks combined b/c the habitat predictor is the only weird one that would have different meanings at different parks.
# > For purposes of trend estimate, it doesn't matter if Dharma on habitat and sd are bad--those are not going to affect trend estimate
rm(list = ls())

### Load packages ----

pkgs <- c("tidyverse", "brms", "partykit", "magrittr", "DHARMa", "DHARMa.helpers", "emmeans", "ggeffects", "tidybayes", "broom", "broom.mixed", "cmdstanr", "glmmTMB", "ggstats", "performance", "AICcmodavg")
installed_pkgs <- pkgs %in% installed.packages()
if (length(pkgs[!installed_pkgs]) > 0) install.packages(pkgs[!installed_pkgs], repos = "https://cloud.r-project.org" , dep=TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

options(mc.cores = parallel::detectCores(), 
        brms.backend = "cmdstanr")

### FUNCTIONS ----
FuncFormatFullObs <- function(park, spec) {
  # Format the df_full_obs for use in GLMM
  # Add a column to indicate a researcher's first year of survey
  # Combine cov. levels
  # Merge in the location data
  
  obs_dat <- df_full_obs
  
  # Add subunits for GUIS & SAAN
  obs_dat$subunit <- obs_dat$unit_code
  obs_dat$unit_code[obs_dat$subunit %in% c("GUIS-FL", "GUIS-MS")] <- "GUIS"
  obs_dat$subunit[obs_dat$unit_code == "SAAN"] <- "MI"
  obs_dat$subunit[obs_dat$location_name %in% c("SAAN01", "SAAN02", "SAAN03", "SAAN04", "SAAN05", "SAAN06")] <- "RA"
  
  # Subset the data and add covariates
  survey_dat_keep <- df_survey_conditions %>%
    dplyr::select(location_name, event_date, researcher, julian_prop, hrs_since_rise, weather_wind)
  loc_dat_keep <- df_locs %>%
    dplyr::select(location_name, hab_type_100, hab_type_200, prop_understory_cover_50, prop_understory_cover_100, prop_understory_cover_200, understory_cover_sd_50, understory_cover_sd_100, understory_cover_sd_200)
  
  # Add and combine covariates as needed
  obs_dat %<>%
    dplyr::filter(unit_code == park & species_code == spec) %>% 
    dplyr::left_join(loc_dat_keep, by = c("location_name")) %>%
    dplyr::left_join(survey_dat_keep, by = c("location_name", "event_date"))
  
  if(park == "VICK") {
    
    ## For VICK, only using Daniel;s data and only data starting 2012 so can avoid messy changes of researchers and also the very odd timing of 5 surveys (the yr_visits had long timespans that overlapped each other) in 2010
    obs_dat %<>% 
      dplyr::filter(researcher == "Twedt, Daniel" & yr >= 2012)
  }
  
  if(park == "BITH") {
    ## REMOVE first three years of data b/c each of these first three years had a different researcher and strong researcher impact
    obs_dat %<>% filter(yr >= 2017)
  }
  
  if(park == "GUIS") {
    obs_dat %<>%
      dplyr::filter(!researcher %in% c("Walker,  Jake", "Sculley,  Mike"))
  }
  
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
    
  return(mod_dat)
}

FuncSub100m <- function() {
  df_finaldat <- readRDS(here::here("Data_out", "df_finaldat.RDS"))
  df_species <- readRDS(here::here("Data_out", "df_species.RDS"))
  
  yr_visit_template <- df_finaldat %>%
    dplyr::select(location_name, event_date, yr) %>%
    dplyr::distinct() %>%
    dplyr::arrange(location_name, yr, event_date) %>%
    dplyr::group_by(location_name, yr) %>%
    dplyr::mutate(
      within_yr_survey = row_number(),
      yr_visit = paste0(yr, "_", within_yr_survey)) %>%
    dplyr::select(-within_yr_survey)
  
  # For each unit, every combination of loc-survey event and bird species (ever found in that unit)
  count_template <- merge(
    df_survey_conditions %>% dplyr::select(unit_code, location_name, event_date) %>% dplyr::distinct() %>% dplyr::inner_join(df_finaldat %>% dplyr::select(location_name, event_date) %>% dplyr::distinct(), by = c("location_name", "event_date")), # the inner join ensures that we are working with the same date range of data from the two data files b/c sometimes the sitecovs data includes more recent survey events than included in the finaldat
    df_finaldat %>% dplyr::select(unit_code, species_code) %>% dplyr::distinct() %>%
      dplyr::inner_join(df_species[c("unit_code", "species_code")] %>% distinct(), by = c("unit_code", "species_code"))) %>% # for each unit separately, all species found in each unit
    dplyr::left_join(df_species %>% dplyr::select(unit_code, species_code, scientific_name, common_name, landbird, park_breeding = present_in_np_species_classified_as_park_breeding), by = c("unit_code", "species_code")
    )
  
  df_full_obs <- df_finaldat %>%
    dplyr::filter(distance_bin_id %in% c(0:2)) %>%
    dplyr::inner_join(df_species[c("unit_code", "species_code")] %>% distinct(), by = c("unit_code", "species_code")) %>%
    dplyr::select(location_name, event_date, species_code, count) %>%
    dplyr::group_by(location_name, event_date, species_code) %>%
    dplyr::summarize(
      sum_indiv = sum(count, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% 
    dplyr::full_join(count_template) %>%
    dplyr::inner_join(yr_visit_template, by = c("location_name", "event_date")) %>%
    dplyr::arrange(location_name, event_date, species_code) %>%
    dplyr::distinct()
  df_full_obs$sum_indiv[is.na(df_full_obs$sum_indiv)] <- 0
  return(df_full_obs)
}
  
  
### Read in formatted data ----
df_survey_conditions <- readRDS(here::here("Data_out", "df_survey_conditions.RDS")) # event covariate data
df_full_obs <- readRDS(here::here("Data_out", "df_full_obs.RDS"))
df_locs <- readRDS(here::here("Data_out", "df_locs.RDS"))

# df_full_obs <- FuncSub100m() # <<<<<< FOR COMPARING TO DISTANCE SAMPLING RESULTS, USING SAME SUBSET OF DATA
# Subset and add covariates ----
park = "PAAL"
spec = "NOCA"

dat <- FuncFormatFullObs(park = park, spec = spec)

table(dat$researcher, dat$yr, dat$subunit) # check how many observers

# Dataset with complete cases only
# PAAL: Dropped 3 records b/c missing weather_wind_comb
# SAAN: Dropped 1 record b/c missing hrs_since_rise
nrow(dat)
colSums(is.na(dat))

dat <- dat[complete.cases(dat),]

# Scale predictors
dat$prop_understory_cover_50_sc <- scale(dat$prop_understory_cover_50)
dat$prop_understory_cover_100_sc <- scale(dat$prop_understory_cover_100)
dat$prop_understory_cover_200_sc <- scale(dat$prop_understory_cover_200)
dat$understory_cover_sd_50_sc <- scale(dat$understory_cover_sd_50)
dat$understory_cover_sd_100_sc <- scale(dat$understory_cover_sd_100)
dat$understory_cover_sd_200_sc <- scale(dat$understory_cover_sd_200)


nrow(dat)
table(dat$yr)
table(dat$yr, dat$location_name)
dat %>% dplyr::group_by(subunit) %>% summarize(tot = sum(sum_indiv))

### >>>> IF GUIS OR SAAN...
# dat%<>% dplyr::filter(subunit=="GUIS-MS")
# unique(dat$location_name)


### Classification trees for exploring covariates ----

## SAAN & GUIS--add subunit
## as needed... researcher + first_yr (BITH, PAAL, GUIS?)
dat_tree <- partykit::ctree(sum_indiv ~hab_type_100_comb + hab_type_200_comb +  prop_understory_cover_50 + prop_understory_cover_100 + prop_understory_cover_200 + understory_cover_sd_50 + understory_cover_sd_100 + understory_cover_sd_200 + weather_wind_comb + julian_prop + hrs_since_rise, data = dat)

plot(dat_tree)
# dat %<>% dplyr::filter(subunit=="MI") # <<<<< SAAN ONLY, AS NEEDED

### Check data patterns ----
ggplot(dat, aes(y = sum_indiv, x = yr_visit, color = location_name, group = location_name)) + geom_line() ->p
plotly::ggplotly(p)
ggplot(dat, aes(x=sum_indiv)) + geom_histogram() + facet_wrap(~subunit)
ggplot(dat, aes(x = sum_indiv)) + geom_histogram()

ggplot(dat, aes(y = sum_indiv, x= hab_type_200_comb)) + geom_boxplot() + facet_wrap(~subunit)
ggplot(dat, aes(y = sum_indiv, x= as.factor(researcher), fill = as.factor(first_yr))) + geom_boxplot()
ggplot(dat, aes(y = sum_indiv, x= weather_wind_comb)) + geom_boxplot()+ facet_wrap(~subunit)

ggplot(dat, aes(y = prop_understory_cover_50, x = as.factor(sum_indiv))) + geom_boxplot()+ facet_wrap(~subunit)
ggplot(dat, aes(y = prop_understory_cover_100, x = as.factor(sum_indiv))) + geom_boxplot()+ facet_wrap(~subunit)
ggplot(dat, aes(y = prop_understory_cover_200, x = as.factor(sum_indiv))) + geom_boxplot()+ facet_wrap(~subunit)
ggplot(dat, aes(y = understory_cover_sd_50, x = as.factor(sum_indiv))) + geom_boxplot()
ggplot(dat, aes(y = understory_cover_sd_100, x = as.factor(sum_indiv))) + geom_boxplot()
ggplot(dat, aes(y = understory_cover_sd_200, x = as.factor(sum_indiv))) + geom_boxplot()
ggplot(dat, aes(x = understory_cover_sd_100, y = understory_cover_sd_200)) + geom_point()
ggplot(dat, aes(y = understory_cover_sd_100, x = prop_understory_cover_100)) + geom_point()
ggplot(dat, aes(y = prop_understory_cover_50, x = prop_understory_cover_200)) + geom_point()
ggplot(dat, aes(y = prop_understory_cover_200, x= hab_type_200_comb)) + geom_boxplot()
ggplot(dat, aes(y = prop_understory_cover_100, x= as.factor(sum_indiv), fill = hab_type_200_comb)) + geom_boxplot() + facet_wrap(~subunit)
ggplot(dat, aes(y = prop_understory_cover_200, x= as.factor(sum_indiv), fill = hab_type_100_comb)) + geom_boxplot()
ggplot(dat, aes(x = prop_understory_cover_200)) + geom_histogram() + facet_wrap(~subunit)
ggplot(dat, aes(x = understory_cover_sd_50)) + geom_histogram()+ facet_wrap(~subunit)
ggplot(dat, aes(x = hrs_since_rise_c)) + geom_histogram()+ facet_wrap(~subunit)
ggplot(dat, aes(x = julian_prop_c)) + geom_histogram()+ facet_wrap(~subunit)

### FINDINGS NOTES ----
# > BITH has many different researchers, so that should be a predictor- for NOCA, that was the primary predictor based on ctree. Researcher is confounded with year effects.
# > BITH only has one level for hab_type_200
# > BITH CARW has one outlier count of 9


### Start with glmmTMB to save time ----

## NOTES ON BITH MODEL
# > exclude habitat bc only one level
# > incl researcher & first_yr b/c many--researcher is confounded with year. So this one is a bit funny b/c a lot of the plots seem to suggest much higher counts in the early years, but the model is attributing those to researcher differences, not to actual high numbers... so BBS models, e.g., say NOCA is declining everywhere but BITH NOCA says nonsign. increase. Will be interesting to see NOCA trend at parks with less researcher turnover
# PAAL--exclude hab b/c only non-forest. Has two researchers.

min_tmb <- glmmTMB(sum_indiv ~ yr_c + (1 | yr_c_f) + (1 | location_name) + prop_understory_cover_200 + understory_cover_sd_200, family = compois, data = dat)

min2_tmb <- glmmTMB(sum_indiv ~ yr_c + (1 | yr_c_f) + (1 | location_name) + hab_type_200_comb + understory_cover_sd_200 + julian_prop_c + hrs_since_rise_c + I(hrs_since_rise_c^2) + prop_understory_cover_200, family = compois, data = dat)

###### >>>>> YEAR MODELS... <<<<<<<<<<< ----
park = "VICK"
spec = "WEVI"
mod <- readRDS(here::here("tmb_mods", paste0("tmb_", park, "_", spec, ".RDS")))
mod$modelInfo$allForm # check formula


### REMOVE RESEARCHER & FIRST YEAR AND YR RE AND YR FE CONTINUOUS
yr_mod <- glmmTMB(sum_indiv ~ 0 + yr_c_f + (1 | location_name) + hab_type_200_comb + 
understory_cover_sd_200_sc + julian_prop_c + hrs_since_rise_c + I(hrs_since_rise_c^2), family = compois, data = dat)

summary(yr_mod)

tmb_zi <- ifelse(yr_mod$modelInfo$allForm$ziformula=="~1", TRUE, FALSE)
ggpredict(yr_mod, terms=c("yr_c_f"), type = ifelse(tmb_zi==TRUE, "zero_inflated", "fixed"))
saveRDS(yr_mod, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "_yearmod.RDS")))
##### <<<<<<<<<<<<<<<<<<<<<<<<

mod_list = tibble::lst(min_tmb, min2_tmb)
aictab(mod_list) # use DHARMa outputs over this, to decide on final model--but confirm that best AIC-supported model and best DHARMa-supported model have nearly identical estimates, esp. for trend
ggstats::ggcoef_compare(mod_list)
# , include = c(yr_c, weather_wind_comb, prop_understory_cover_200, understory_cover_sd_100, hrs_since_rise_c)

### Check the model ----
mod <- min_tmb
summary(mod)
# ggstats::ggcoef_table(mod, intercept = FALSE, table_witdhs = c(2, 1)) + labs(title = "XXX")

## pp-checks
performance::check_predictions(mod, check_range = TRUE)
performance::check_collinearity(mod)
performance::check_autocorrelation(mod)
# performance::check_distribution(mod)
### Dharma checks for glmmTMB ----

simres <- DHARMa::simulateResiduals(fittedModel = mod, n = 1000, refit = FALSE, plot = TRUE) # re-simulate all levels, to test model structure as a whole. The default is resimulating on all levels (unconditional sims) instead of simulating conditional on the fitted random effects.
simres <- DHARMa::simulateResiduals(fittedModel = mod, n = 1000, refit = FALSE, re.form = NULL, plot = TRUE) #re-simulate conditional on the fitted random effects, i.e., the RE's are not re-simulated. More sensitive for detecting problems, but it doesn't test all possible options.

DHARMa::testZeroInflation(simres)
DHARMa::testDispersion(simres) # test for over/under dispersion
DHARMa::testOutliers(simres, type = "bootstrap")

# Plot residuals vs each independent factor
DHARMa::plotResiduals(simres, as.factor(dat$yr_c))
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

### SAVE THE FINAL MODEL ----
# park = "GUIS-FL"
saveRDS(mod, here::here("tmb_mods", paste0("tmb_", park, "_", spec, ".RDS")))
# saveRDS(mod, here::here("tmb_mods", paste0("tmb_", park, "_", spec, "limit100m.RDS")))

### PLOT glmmTMB PREDICTIONS ----
em=emmeans::emtrends(mod, specs = c("yr_c"), var = "yr_c") # trends at mean year. Ignores a 'type' argument, it seems
summary(em, infer=c(TRUE,TRUE),null=0) %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  dplyr::mutate(signif = p_value < 0.05)

ggpredict(mod, terms=c("yr_c")) # I think ggpredict uses the most common for factor levels, while ggemmeans will do some kind of averaging. Can get separate estimates for each researcher by terms = c("yr_c", "researcher")
x<-ggemmeans(mod, terms="yr_c")

x<- ggpredict(mod, terms = list(
  yr_c = seq(-5, 5, by =1)
  )
  )

plot(x, facets = FALSE)

