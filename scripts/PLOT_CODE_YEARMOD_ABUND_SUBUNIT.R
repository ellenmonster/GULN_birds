### Modified from PLOT_CODE_ABUND on 5/15/24
# Uses predict_response so can incorporate ZI
# Only compares year estimates for glmmTMB 1 vs 2 visit and Stan pcount (2 visit)
# All models are limited to 100m
# No trend models included
# NOTE that the stan_pcount estimates hold other covariates to their MEDIANS or factor reference levels; glmmTMB "mean_reference" sets non-focal predictors to the MEANS and reference levels (or "most common" value for character vectors). In FINAL, need to make it consistent.

### Load packages ----

pkgs <- c("tidyverse", "ggeffects", "ubms", "patchwork", "here")
installed_pkgs <- pkgs %in% installed.packages()
if (length(pkgs[!installed_pkgs]) > 0) install.packages(pkgs[!installed_pkgs], repos = "https://cloud.r-project.org" , dep=TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

source(here::here("scripts", "ANALYSIS_STARTUP.R"))

### NOTE:
# Discussion on getting estimates per year... https://github.com/rvlenth/emmeans/issues/118

### ADD TITLE (PARK-BIRD), CAPTION, ACTUAL YEARS, BETTER LEGEND TITLE. ALSO MENTION THAT FROM 2019 ON, THE POINTS ARE THE SUM OF TWO SURVEYS. Then for bar graphs, show multiple species. Mention that in xx years, not all locs were surveyed in the park each year and from 2019 on, there were two surveys so that's why there are bigger circles in latter yrs but the important thing is the relative size of each circle within year

### FUNCTIONS
FuncTMB <- function(mod, mod_type, num_visits) {
  # ggeffects::predict_response(mod, terms = c("yr_c_f", "subunit"), margin = "mean_reference", type = "zero_inflated")) # > Note that to get emmeans, 'margin = "marginalmeans"' but then can't do zero-inflation with that
  pred <- ggeffects::predict_response(model = mod, margin = "mean_reference", terms=c("yr_c_f"), type = ifelse(tmb_zi, "zero_inflated", "fixed")) %>% 
    as.data.frame() %>%
    janitor::clean_names() %>%
    dplyr::mutate(
      yr_c = as.numeric(as.character(x)),
      inf_ci = is.infinite(conf_high),
      n_visits = num_visits,
      mod = mod_type) %>%
    dplyr::select(
      mod,
      n_visits,
      yr_c,
      yr_predicted = predicted,
      yr_conf_low = conf_low,
      yr_conf_high = conf_high,
      inf_ci)
  
  pred$yr <- pred$yr_c+(range(dat$yr) %>% median())
  pred$yr_conf_high[pred$inf_ci == TRUE] <- pred$yr_conf_low[pred$inf_ci == TRUE]
  
  return(pred)
}

### User entry----
park <- "BITH"  # "BITH", "JELA", "PAAL", "VICK", 
spec <- "MODO" 
# spec_name <- "Northern cardinal"
  
### Get the subset of raw data----
dat <- FuncFormatFullObs(park = park, spec = spec)
  
### Read in 100m models----
mod_tmb_yearmod <- readRDS(here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_100m_yearmod.RDS")))
mod_tmb_1visit_yearmod <- readRDS(here::here("tmb_mods_yearmod", paste0("tmb_", park, "_", spec, "_100m_1visit_yearmod.RDS")))

mod_stanpcount_yearmod <- readRDS(here::here("stan_pcount_mods", paste0("pcount_", park, "_", spec, "_100m_yearmod.RDS")))

### TMB-predicted values----
tmb_zi <- ifelse(mod_tmb_yearmod$modelInfo$allForm$ziformula=="~1", TRUE, FALSE) # CHECK FOR ZERO INFLATION!

## glmmTMB-predicted values----
pred_tmb1 <- FuncTMB(mod = mod_tmb_1visit_yearmod, mod_type = "GLMM (index of abundance)", num_visits = 1)
pred_tmb2 <- FuncTMB(mod = mod_tmb_yearmod, mod_type = "GLMM (index of abundance)", num_visits = 2)

### Stan_pcount-predicted values----

pred_stan <- plot_marginal(mod_stanpcount_yearmod, "state")$data %>%
   dplyr::mutate(
     yr_c = as.numeric(as.character(covariate)),
     inf_ci = FALSE,
     mod = "N-MIXTURE (superpopulation abundance est'd from repeat visits)",
     n_visits = 2) %>%
   dplyr::select(
     mod,
     n_visits,
     yr_c,
     yr_predicted = mn,
     yr_conf_low = lower,
     yr_conf_high = upper,
     inf_ci)
pred_stan$yr <- pred_stan$yr_c+(range(dat$yr) %>% median())

### Plot----

# dataframe for plot
pdat <- rbind.data.frame(pred_tmb1, pred_tmb2, pred_stan)

### FULL PLOT ----
(p <- ggplot(data = pdat, aes(group = n_visits)) + 
  geom_errorbar(
    mapping = aes(x = yr, y = yr_predicted, ymin = yr_conf_low, ymax = yr_conf_high, linetype = as.character(n_visits), color = inf_ci), width = 0, position = position_dodge(0.3)) +
  geom_point(
    mapping = aes(x = yr, y = yr_predicted, color = inf_ci), size = 2.5, position = position_dodge(0.3)) +
   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), name = "Infinite 95%CI?") +
   scale_linetype_manual(values = c("2" = "solid", "1" = "dotted"), name = "Max # of visits/yr") +
  labs(title = paste(toupper(park), toupper(spec), "within 100m radius"), x = "Year", y = "Birds per point count   (NOTE: Y-axis scales differ)") +
 theme_bw() +
  facet_wrap(~mod, ncol = 1, scales = "free_y"))