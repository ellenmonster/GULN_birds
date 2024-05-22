### MAKE SURE THE GGEMMEANS ESTIMATES AND CI'S ACCOUNT FOR ANY DISPERSION ESTIMATE AND ZERO-INFLATION, WHEN APPROPRIATE!!! <-----------

### This script generates faceted plots for comparing the annual estimates using different approaches. Each approach (e.g., glmmTMB, gDISTSAMP, etc.) is in a different facet


rm(list=ls())

library(tidyverse)
library(ggeffects)
library(ubms)
# library(patchwork)
library(here)

source(here::here("scripts", "ANALYSIS_STARTUP.R"))

### User entry----
park <- "JELA" 
subunit_name <- NULL #<<<<<<<<<<<<<<<< CHANGE SUBUNIT: NULL, GUIS-FL, GUIS-MS, MI, RA
spec <- "WEVI" 
  
### Read in models----
mod_tmb_yearmod <- readRDS(here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_yearmod.RDS")))
mod_tmb_1visit_yearmod <- readRDS(here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_1visit_yearmod.RDS")))

mod_tmb_100m_yearmod <- readRDS(here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_100m_yearmod.RDS")))
mod_tmb_100m_1visit_yearmod <- readRDS(here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_100m_1visit_yearmod.RDS")))

mod_pcount_yearmod <- readRDS(here::here("stan_pcount_mods", paste0("pcount_", park, "_", spec, "_yearmod.RDS")))
mod_pcount_100m_yearmod <- readRDS(here::here("stan_pcount_mods", paste0("pcount_", park, "_", spec, "_100m_yearmod.RDS")))

mod_dist_1visit_yearmod <- readRDS(here::here("gdistsamp_mods", paste0("gdistsamp_", park, "_", spec, "_1visit_yearmod.RDS")))
mod_dist_yearmod <- readRDS(here::here("gdistsamp_mods", paste0("gdistsamp_", park, "_", spec, "_2visit_yearmod.RDS")))

### Stan_pcount-predicted values ----
# Get the year estimates for STAN PCOUNT

FuncPredPcount <- function(mod, bin_limit) {
  if(is.null(subunit_name)) {
    pcount_2visit_pred <- ubms::plot_effects(object = mod, submodel = "state", covariate = "yr_f", draws = 2000)$data %>%
      dplyr::mutate(
        yr = as.integer(as.character(covariate)),
        lim = bin_limit,
        method = "Nmix_pcount(2)",
        visits = 2) %>%
      dplyr::select(
        lim,
        method,
        visits,
        yr,
        yr_predicted = mn,
        yr_conf_low = lower,
        yr_conf_high = upper)
  } else {
    # >>>>>> CHECK FOR RESEARCHER REFERENCE VALUE
    # Use median (not mean) and reference values below to be comparable with predictions for other methods
    nd = data.frame(yr_sc = sort(unique(mod_pcount_yearmod@data@siteCovs$yr_sc)), 
                    yr_f = sort(unique(mod_pcount_yearmod@data@siteCovs$yr_f)),
                    yr_c = sort(unique(mod_pcount_yearmod@data@siteCovs$yr_c)),
                    subunit = subunit_name,
                    first_yr = 0,
                    weather_wind_comb = "calm",
                    hab_type_100_comb = ifelse(park == "SAAN", "not_forest", "forest"),
                    hab_type_200_comb = ifelse(park == "SAAN", "not_forest", "forest"), 
                    prop_understory_cover_50_sc = median(mod_pcount_yearmod@data@siteCovs$prop_understory_cover_50_sc, na.rm = TRUE),
                    prop_understory_cover_100_sc = median(mod_pcount_yearmod@data@siteCovs$prop_understory_cover_100_sc, na.rm = TRUE),
                    prop_understory_cover_200_sc = median(mod_pcount_yearmod@data@siteCovs$prop_understory_cover_200_sc, na.rm = TRUE),
                    understory_cover_sd_50_sc = median(mod_pcount_yearmod@data@siteCovs$understory_cover_sd_50_sc, na.rm = TRUE),
                    understory_cover_sd_100_sc = median(mod_pcount_yearmod@data@siteCovs$understory_cover_sd_100_sc, na.rm = TRUE),
                    understory_cover_sd_200_sc = median(mod_pcount_yearmod@data@siteCovs$understory_cover_sd_200_sc, na.rm = TRUE),
                    julian_prop_c = median(mod_pcount_yearmod@data@obsCovs$julian_prop_c, na.rm = TRUE),
                    hrs_since_rise_c = median(mod_pcount_yearmod@data@obsCovs$hrs_since_rise_c, na.rm = TRUE))
    glimpse(nd)
    
    temp_pcount_pred_yearmod <- ubms::predict(mod_pcount_yearmod, submodel = "state", newdata = nd, re.form = NA) %>%
      janitor::clean_names()
    
    pcount_2visit_pred <- cbind(nd$yr_f, temp_pcount_pred_yearmod) %>%
      dplyr::mutate(
        yr = as.integer(as.character(`nd$yr_f`)),
        lim = bin_limit,
        method = "Nmix_pcount(2)",
        visits = 2) %>%
      dplyr::select(
        lim,
        method,
        visits,
        yr,
        yr_predicted = predicted,
        yr_conf_low = x2_5_percent,
        yr_conf_high = x97_5_percent)
    
    if(subunit_name == "RA") {
      pcount_2visit_pred$yr_predicted[pcount_2visit_pred$yr==2013] <- NA
      pcount_2visit_pred$yr_conf_low[pcount_2visit_pred$yr==2013] <- NA
      pcount_2visit_pred$yr_conf_high[pcount_2visit_pred$yr==2013] <- NA
    }
  }
  
  
  # If confidence interval is infinite, then set infinite_CI = TRUE (color point red in plot) and get rid of CI
  pcount_2visit_pred$infinite_CI = FALSE
  pcount_2visit_pred$infinite_CI[is.infinite(pcount_2visit_pred$yr_conf_high)] <- TRUE
  pcount_2visit_pred$yr_conf_high[is.infinite(pcount_2visit_pred$yr_conf_high)] <- pcount_2visit_pred$yr_conf_low[is.infinite(pcount_2visit_pred$yr_conf_high)]
  
  return(pcount_2visit_pred)
}

pcount_2visit_pred <- FuncPredPcount(mod = mod_pcount_yearmod, bin_limit = "unlimited")
pcount_100m_2visit_pred <- FuncPredPcount(mod = mod_pcount_100m_yearmod, bin_limit = "100m")
(yr_vec = as.integer(pcount_2visit_pred$yr))

### TMB-predicted values----
# ggemmeans to get year predictions
FuncTMBYrPredict <- function(mod, meth, num_visit, bin_limit) {
  # setting 'typical = "median"' to be comparable to predictions from the other methods
  if(is.null(subunit_name)) {
    temp_tmb_yearmod_pred <- ggemmeans(mod, terms=c("yr_c_f"), typical = "median", type = "fixed") 
  } else {
    temp_tmb_yearmod_pred <- ggemmeans(mod, terms=c("yr_c_f", "subunit"), typical = "median", type = "fixed") %>%
      dplyr::filter(group == subunit_name)
  }# type = ifelse(tmb_zi==TRUE, "fe.zi", "fixed"))) #### ZERO-INFLATION PART DOESN'T WORK!!!
  
  tmb_pred <- temp_tmb_yearmod_pred %>%
     as.data.frame() %>%
     janitor::clean_names() %>%
    dplyr::mutate(
      yr_c = as.numeric(as.character(x)),
      lim = bin_limit,
      method = meth,
      visits = num_visit) %>%
      dplyr::select(
        lim,
        method,
        visits,
        yr_c,
        yr_predicted = predicted,
        yr_conf_low = conf_low,
        yr_conf_high = conf_high) %>%
    dplyr::arrange(yr_c)
    
  tmb_pred$yr <- yr_vec
  tmb_pred$yr_c <- NULL
  tmb_pred$infinite_CI = FALSE
  # tmb_pred %<>% dplyr::full_join(temp_tmb_pred, by = c("yr_c"))
  
  # If confidence interval is infinite, then set infinite_CI = TRUE (color point red in plot) and get rid of CI
  tmb_pred$infinite_CI[is.infinite(tmb_pred$yr_conf_high)] <- TRUE
  tmb_pred$yr_conf_high[is.infinite(tmb_pred$yr_conf_high)] <- tmb_pred$yr_conf_low[is.infinite(tmb_pred$yr_conf_high)]
  
  tmb_pred
}

(tmb_1visit_pred <- FuncTMBYrPredict(mod_tmb_1visit_yearmod, meth = "GLMM(1)", num_visit = 1, bin_limit = "unlimited"))
(tmb_2visit_pred <- FuncTMBYrPredict(mod_tmb_yearmod, meth = "GLMM(2)", num_visit = 2, bin_limit = "unlimited"))

(tmb_100m_1visit_pred <- FuncTMBYrPredict(mod_tmb_100m_1visit_yearmod, meth = "GLMM(1)", num_visit = 1, bin_limit = "100m"))
(tmb_100m_2visit_pred <- FuncTMBYrPredict(mod_tmb_100m_yearmod, meth = "GLMM(2)", num_visit = 2, bin_limit = "100m"))

### Distance-predicted values----
FuncDistYrPred <- function(mod, meth, num_visit, yrs = yr_vec) {
  if(is.null(subunit_name)) {
    temp_pred <- mod$gdist_best_mod_pred
  } else {
    temp_pred <- mod$gdist_best_mod_pred %>%
      dplyr::filter(subunit == subunit_name)
  }
  
  names(temp_pred) <- tolower(names(temp_pred))
  temp_pred$yr <- yrs
  pred <- temp_pred %>%
    dplyr::mutate(
      lim = "100m",
      method = meth,
      visits = num_visit) %>%
    dplyr::select(
      lim,
      method,
      visits,
      yr,
      yr_predicted = predicted,
      yr_conf_low = lower,
      yr_conf_high = upper)
  
  # If confidence interval is infinite, then set infinite_CI = TRUE (color point red in plot) and get rid of CI
  pred$infinite_CI = FALSE
  pred$infinite_CI[is.infinite(pred$yr_conf_high)|pred$yr_conf_high>10000] <- TRUE
  pred$yr_conf_high[is.infinite(pred$yr_conf_high)|pred$yr_conf_high>10000] <- pred$yr_conf_low[is.infinite(pred$yr_conf_high)|pred$yr_conf_high>10000]
  pred
}

(dist_1visit_pred <- FuncDistYrPred(mod_dist_1visit_yearmod, meth = "Distance(1)", num_visit = 1))
(dist_2visit_pred <- FuncDistYrPred(mod_dist_yearmod, meth = "Distance(2)", num_visit = 2))

### Combine data----
plot_dat <- rbind.data.frame(
  tmb_1visit_pred,
  tmb_2visit_pred,
  tmb_100m_1visit_pred,
  tmb_100m_2visit_pred,
  # dist_1visit_pred,
  # dist_2visit_pred,
  pcount_2visit_pred,
  pcount_100m_2visit_pred
  )
plot_dat$visits <- as.factor(plot_dat$visits)
plot_dat
### FULL PLOT ----

if(is.null(subunit_name)) {subunit_name <- park}
(p <- ggplot(plot_dat, aes(x = yr, y = yr_predicted, color = lim)) +
   geom_point(aes(shape = infinite_CI), size =1.5, position = position_dodge(0.75)) + # >>>> DODGE WHEN INCLUDING 100M
   geom_errorbar(aes(ymin = yr_conf_low, ymax = yr_conf_high, linetype = visits), width = 0, position = position_dodge(0.75)) + # >>>> DODGE WHEN INCLUDING 100M
   scale_linetype_manual(values = c("2" = "solid","1" = "twodash"), guide = "none") +
   scale_shape_manual(values = c(19, 21), guide = "none") +
  scale_x_continuous(breaks = yr_vec) +
  scale_y_continuous(expand = expansion(add = c(1,1)), limits = c(0, NA)) +
    scale_color_manual(values = c("black", "#009E73")) +
  labs(title = paste(subunit_name, spec), subtitle = "Error bars are 95%CI; dashed means single visit. Open circles have infinite CI's. Black = limited to 100m radius.", x = "Year", y = "Index of abundance (# of birds/point count)") +
  theme_bw() +
  facet_wrap(vars(method)))
    
### SAVE PLOT ----

saveRDS(plot_dat, here::here("model_compare_plots", paste0("plotdat_abund_compare_", subunit_name, "_", spec, ".RDS")))
