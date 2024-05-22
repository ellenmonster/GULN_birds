library(tidyverse)
library(ggeffects)
library(ubms)
library(patchwork)

source("TEST_STARTUP.R")

### NOTE:
# Discussion on getting estimates per year... https://github.com/rvlenth/emmeans/issues/118

### ADD TITLE (PARK-BIRD), CAPTION, ACTUAL YEARS, BETTER LEGEND TITLE. ALSO MENTION THAT FROM 2019 ON, THE POINTS ARE THE SUM OF TWO SURVEYS. Then for bar graphs, show multiple species. Mention that in xx years, not all locs were surveyed in the park each year and from 2019 on, there were two surveys so that's why there are bigger circles in latter yrs but the important thing is the relative size of each circle within year

### Read in formatted data ----
df_survey_conditions <- readRDS(here::here("Data_out", "df_survey_conditions.RDS")) # event covariate data
df_full_obs <- readRDS(here::here("Data_out", "df_full_obs.RDS"))
df_locs <- readRDS(here::here("Data_out", "df_locs.RDS"))

# DONE: 
#> NOCA: 
#> WEVI: 
### User entry----
park <- "GUIS"  # "GUIS"
spec <- "NOCA" # "NOCA", "WEVI"
# spec_name <- "Northern cardinal"

mod_list <- list()
  
  ### Get the subset of raw data----
  dat <- FuncFormatFullObs(park = park, spec = spec)
  dat <- dat[complete.cases(dat),]
  table(dat$researcher, dat$yr, dat$subunit) # check how many observers
  
### Read in models----
mod_tmb_yearmod <- readRDS(here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_yearmod.RDS")))
mod_tmb <- readRDS(here::here("tmb_mods", paste0("tmb_", park, "_", spec, ".RDS")))

mod_pcount_yearmod <- readRDS(here::here("stan_pcount_mods", paste0("pcount_", park, "_", spec, "_yearmod.RDS")))
mod_pcount <- readRDS(here::here("stan_pcount_mods", paste0("pcount_", park, "_", spec, ".RDS")))

### TMB-predicted values----
tmb_zi <- ifelse(mod_tmb$modelInfo$allForm$ziformula=="~1", TRUE, FALSE) # CHECK FOR ZERO INFLATION!

## ggemmeans to get trend ribbon averaged over researchers
# >>>>>>>>>>>>> FOR ZERO-INFLATED, USE type ="zi_fixed" !!! >>>>> NO, THIS DOESN'T WORK!
(temp_tmb_pred <- ggemmeans(mod_tmb, terms=c("yr_c", "subunit"), type = "fixed"))
                            
                            # type = ifelse(tmb_zi==TRUE, "fe.zi", "fixed"))) #### ZERO-INFLATION PART DOESN' WORK!!!

temp_tmb_pred %<>% 
  as.data.frame 
temp_tmb_pred$yr <- temp_tmb_pred$x+(range(dat$yr) %>% median())
temp_tmb_pred %<>%  
  janitor::clean_names() %>%
  dplyr::mutate(yr_c = as.integer(x)) %>%
  dplyr::select(
    yr,
    yr_c,
    subunit = group,
    trend_predicted = predicted,
    trend_conf_low = conf_low,
    trend_conf_high = conf_high)

## ggemmeans to get year predictions
(temp_tmb_pred_yearmod <- ggemmeans(mod_tmb_yearmod, terms=c("yr_c_f", "subunit"), type = "fixed") %>% # # type = ifelse(tmb_zi==TRUE, "fe.zi", "fixed"))) #### ZERO-INFLATION PART DOESN' WORK!!!
  as.data.frame() %>%
    janitor::clean_names())
tmb_pred <- temp_tmb_pred_yearmod %>%
  dplyr::mutate(yr_c = as.integer(as.character(x))) %>%
  dplyr::select(
    yr_c,
    subunit = group,
    yr_predicted = predicted,
    yr_conf_low = conf_low,
    yr_conf_high = conf_high)

tmb_pred %<>% dplyr::full_join(temp_tmb_pred, by = c("yr_c", "subunit"))

tmb_pred

### Stan_pcount-predicted values----

## REMEMBER THAT YR_SC IS MEAN-CENTERED, WHILE YR_C IS MEDIAN-CENTERED!!!
mod_pcount@call$formula # <<<<< CHECK WHAT PREDICTING ON

table(mod_pcount@data@siteCovs$hab_type_200_comb)
# table(mod_pcount@data@siteCovs$researcher)
        
set_subunit = "GUIS-FL"#<<<<<<<<<<<<<<<< CHANGE SUBUNIT: GUIS-FL, GUIS-MS, MI, RA
nd = data.frame(yr_sc = sort(unique(mod_pcount@data@siteCovs$yr_sc)), 
                yr_f = sort(unique(mod_pcount@data@siteCovs$yr_f)),
                yr_c = sort(unique(mod_pcount@data@siteCovs$yr_c)),
                subunit = set_subunit, 
                first_yr = 0,
                weather_wind_comb = "calm",
                hab_type_100_comb = ifelse(park == "SAAN", "not_forest", "forest"),
                hab_type_200_comb = ifelse(park == "SAAN", "not_forest", "forest"), 
                prop_understory_cover_50_sc = 0,
                prop_understory_cover_100_sc = 0,
                prop_understory_cover_200_sc = 0,
                understory_cover_sd_50_sc = 0,
                understory_cover_sd_100_sc = 0,
                understory_cover_sd_200_sc = 0,
                julian_prop_c = 0,
                hrs_since_rise_c = 0, 
                prop_understory_cover_50 = mean(mod_pcount@data@siteCovs$prop_understory_cover_50),
                understory_cover_sd_50 = mean(mod_pcount@data@siteCovs$understory_cover_sd_50),
                understory_cover_sd_200 = mean(mod_pcount@data@siteCovs$understory_cover_sd_200),
                prop_understory_cover_200 = mean(mod_pcount@data@siteCovs$prop_understory_cover_200))
nd

# Get trend predictions from stan_pcount n-mixture
mod_pcount@call$formula

temp_pcount_pred <- ubms::predict(mod_pcount, submodel = "state", newdata = nd, re.form = NA) %>%
  janitor::clean_names()

temp_pcount_pred %<>%
  dplyr::select(
    trend_predicted = predicted,
    trend_conf_low = x2_5_percent,
    trend_conf_high = x97_5_percent)
temp_pcount_pred$yr <- sort(unique(mod_pcount@data@siteCovs$yr))
temp_pcount_pred$subunit = set_subunit
temp_pcount_pred

# Get the year estimates
mod_pcount_yearmod@call$formula
temp_pcount_pred_yearmod <- ubms::predict(mod_pcount_yearmod, submodel = "state", newdata = nd, re.form = NA) %>%
  janitor::clean_names()

temp_pcount_pred_yearmod %<>%
  dplyr::select(
    yr_predicted = predicted,
    yr_conf_low = x2_5_percent,
    yr_conf_high = x97_5_percent)

pcount_pred <- cbind(temp_pcount_pred, temp_pcount_pred_yearmod)
pcount_pred 

if(set_subunit == "RA") {
  pcount_pred$yr_predicted[pcount_pred$yr==2013] <- NA
  pcount_pred$yr_conf_low[pcount_pred$yr==2013] <- NA
  pcount_pred$yr_conf_high[pcount_pred$yr==2013] <- NA
  }

yr_vec = as.integer(pcount_pred$yr)
### FULL PLOT ----
(p <- ggplot() + 
  geom_count(
    data = subset(dat, subunit == set_subunit), 
    mapping = aes(x = yr, y = sum_indiv), alpha = 0.15) +
  scale_size_area(max_size = 8) +
   scale_size(breaks = c(5,10,15,20,25,30),labels = c(5,10,15,20,25,30), limits = c(1,30),range=c(1,8)) + # NOCA
   # scale_size(breaks = c(10,20,30,40),labels = c(10,20,30,40), limits = c(1,40),range=c(1,8)) + # WEVI
# geom_ribbon( # glmmTMB TREND RESULTS AS RIBBON
#   data = tmb_pred,
#   mapping = aes(x = yr, y = trend_predicted, ymin = trend_conf_low, ymax = trend_conf_high),
#   alpha = 0.1,
#   fill = "blue") +
# geom_line(data =tmb_pred, aes(x=yr, y= trend_predicted), color = "blue", linetype = "dotted") +
# geom_ribbon( # N-MIXTURE TREND RESULTS AS RIBBON
#   data = pcount_pred,
#   mapping = aes(x = yr, y= trend_predicted, ymin = trend_conf_low, ymax = trend_conf_high),
#   alpha = 0.1,
#   fill = "red") +
# geom_line(data = pcount_pred, aes(x=yr, y=trend_predicted), color = "red", linetype = "dotted") +
  geom_errorbar( # glmmTMB YEAR RESULTS AS POINTS W/ERROR BARS
    data = subset(tmb_pred, subunit == set_subunit),
    mapping = aes(x = yr, y = yr_predicted, ymin = yr_conf_low, ymax = yr_conf_high), width = 0, color = "blue") +
  geom_point(
    data = subset(tmb_pred, subunit == set_subunit),
    mapping = aes(x = yr, y = yr_predicted), size = 2.5, color = "blue") +
  geom_errorbar( # N-MIXTURE RESULTS AS POINTS W/ERROR BARS
    data = pcount_pred,
    mapping = aes(x = yr, y = yr_predicted, ymin = yr_conf_low, ymax = yr_conf_high), width = 0, color = "red") +
  geom_point(
    data = pcount_pred,
    mapping = aes(x = yr, y = yr_predicted), size = 2.5, color = "red") +
  scale_x_continuous(breaks = yr_vec) +
  scale_y_continuous(expand = expansion(add = c(1,1)), limits = c(0, NA)) +
  labs(title = ifelse(park == "GUIS", set_subunit, paste0(park, "-", set_subunit)), x = "Year", y = "Index of abundance (# of birds/point count)", size = "# of point counts") +
  geom_vline(xintercept = 2018.5, linewidth = 2, color = "lightgray") + theme_bw())
  # geom_text() + annotate("text", label = "Every location surveyed 2X/yr since 2019", x = 2018.6, hjust = 0, y = 15.7, size = 2.75) +
  # geom_segment(aes(x=2018.5, y=15.2, xend=2022, yend=15.2), arrow = arrow(length=unit(.2, 'cm'))) + 
  # geom_text() + annotate("text", label = "Half of sites surveyed each year before 2018", x = 2017.6, hjust = 1, y = 15.7, size = 2.75) +
  # geom_segment(aes(x=2017.5, y=15.2, xend=2014, yend=15.2), arrow = arrow(length=unit(.2, 'cm')))) 
    p
### SAVE PLOT ----
subunit <- "SAANRA"

saveRDS(p, here::here("RESULTS", "abund_plots", paste0("plot_abund_", subunit, "_", spec, ".RDS")))



# Put the plots together
# Western park units - SAAN, PAAL, BITH (NOCA=30)
# Eastern park units - VICK, JELA, GUIS (NOCA=30)
title = paste0(park, " - ", spec_name), subtitle = "gray circles = actual counts, blue = relative abundance, red = 'detection-corrected' abundance") 

