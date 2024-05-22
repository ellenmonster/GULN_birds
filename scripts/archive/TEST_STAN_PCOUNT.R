### GULN BIRDS - NMIXTURE IN UBMS
### CARW, NOCA, 
### NOTES:
# > Distinguishes 0-counts from non-visits. Non-visits are listed as NA, just to have the same number of columns per site (one col per visit)

### CHECK ON THIS:
# > Is it really treating missing data (NA) correctly and different from 0??

### TO DO:
# > https://rdrr.io/cran/unmarked/man/detFuns.html   plotting detection functions in 'unmarked'

### USEFUL:
# https://cran.r-project.org/web/packages/ubms/vignettes/ubms.html#understanding-the-ubms-output-summary

rm(list = ls())

### Load packages ----
library(tidyverse)
library(magrittr)
library(ubms)

source(here::here("scripts", "TEST_STARTUP.R"))


### Read in formatted data ----
df_full_obs <- readRDS(here::here("Data_out", "df_full_obs.RDS"))

# Subset and add covariates ----
park = "VICK"
spec = "NOCA"

dat <- FuncFormatFullObs(park = park, spec = spec)
table(dat$researcher, dat$yr, dat$subunit) # check how many observers

# Dataset with complete cases only
# PAAL: Dropped 3 records b/c missing weather_wind_comb
# SAAN: Dropped 1 record b/c missing hrs_since_rise
nrow(dat)
colSums(is.na(dat))

dat <- dat[complete.cases(dat),]
nrow(dat)

### N-MIXTURE ANALYSIS ----
table(dat$yr_visit) # make sure nothing funny

# NOTE some assumptions of N-mixture models. One is that we are not overcounting--we are only undercounting. That is, in each survey, we detect only a subset of the birds that are present. If some birds are being double-counted (seems likely), this could violate assumptions. Also assuming that probability of detection is homogeneous. N-mixture can substantially overestimate abundance, in the face of unmodeled heterogeneity in the data.

dat_Nmix <- dat

### Format data for 'ubms' ----
  
(primary_surv <- sort(unique(dat_Nmix$yr))) # these are the years of survey

(secondary_surv <- 1:max(as.numeric(sub(".*_", "", dat_Nmix$yr_visit)), na.rm = TRUE)) # sequence from one to maximum number of repeat survey visits (over entire dataset)

template <- expand_grid(location_name = sort(unique(dat_Nmix$location_name)), yr = primary_surv, visit = secondary_surv) %>%
  dplyr::mutate(yr_visit = paste(yr, visit, sep = "_")) %>%
  dplyr::select(-yr, -visit)

expand_dat <- template %>%
  left_join(dat_Nmix, by = c("location_name", "yr_visit")) %>%
  dplyr::select(-yr)

table(expand_dat$yr_visit) # check that every year has the same number of repeat visits and the same number of locs. The additional records are NA, indicating missing data.

### Response variable
y_st <- FuncMatrix_Nmix(dat = expand_dat, num_surv = max(secondary_surv), cov_vec = c("sum_indiv"), response_mat = TRUE, stack = TRUE) # convert it to stacked format, specifying that each year has 2 surveys

dim(y_st)
head(y_st)

### FOR GUIS AND SAAN INCL SUBUNIT!!!
site_covs <- dat_Nmix %>%
  dplyr::select(location_name, hab_type_100_comb, hab_type_200_comb, prop_understory_cover_50, prop_understory_cover_100, prop_understory_cover_200, understory_cover_sd_50, understory_cover_sd_100, understory_cover_sd_200) %>%
  dplyr::distinct() %>%
  dplyr::mutate_if(is.character, as.factor)

# table(site_covs$hab_type) # are there enough sites in each habitat type to use this; what about in habitat type comb?
# table(site_covs$hab_type_comb)

site_covs_st <- do.call("rbind", replicate(
  length(unique(na.omit(dat_Nmix$yr))), site_covs, simplify = FALSE)) # stack the data, one rep per year
site_covs_st$yr <- rep(sort(unique(dat_Nmix$yr)), each = length(unique(na.omit(dat_Nmix$location_name)))) # for year trend
site_covs_st$yr_c <- scale(site_covs_st$yr, center = TRUE, scale = FALSE) 
site_covs_st$yr_sc <- scale(site_covs_st$yr) %>% as.numeric()
site_covs_st$yr_f <- as.factor(site_covs_st$yr) # for year random effect
dim(site_covs_st)
head(site_covs_st)

### >>>>>Observation covariates -- INCLUDE RESEARCHER AND FIRST YEAR WHERE NEEDED!!!
obs_covs_st <- FuncMatrix_Nmix(dat = expand_dat, num_surv = max(secondary_surv), cov_vec = c("researcher", "first_yr", "julian_prop_c", "hrs_since_rise_c", "weather_wind_comb"), response_mat = FALSE, stack = TRUE) # don't use weather_noise_comb b/c a lot of missing data
names(obs_covs_st)

umf <- unmarkedFramePCount(
  y=y_st, 
  siteCovs=site_covs_st, # this is where site and year are included
  obsCovs=obs_covs_st)
summary(umf)

### START THE MODELS ----
# Only allows Poisson distribution
# Site and year random effects, plus year fixed effect
# Each model takes ~30+ minutes!

### STAN_PCOUNT

# This one uses habitat variables in both the detection and abundance equations
m_tmb <- stan_pcount( # detection covariates first
  ~researcher + first_yr + weather_wind_comb + prop_understory_cover_200 +  understory_cover_sd_200 + julian_prop_c + hrs_since_rise_c + I(hrs_since_rise_c^2)
  ~yr_sc + (1|location_name) + (1|yr_f) + prop_understory_cover_200 +  understory_cover_sd_200,
  data = umf, 
  chains=4, 
  iter=4000, 
  thin = 2,
  cores = 4) # cov of detection, followed by abundance


# FOR SUBUNITS...
m_tmb <- stan_pcount( # detection covariates first
  ~julian_prop_c + prop_understory_cover_200:hab_type_200_comb:subunit + prop_understory_cover_200 + hab_type_200_comb
  ~yr_c*subunit  + (1|location_name) + (1|yr_f) +prop_understory_cover_200:hab_type_200_comb:subunit + prop_understory_cover_200 + hab_type_200_comb,
  data = umf, 
  chains=4, 
  iter=4000, 
  thin = 2,
  cores = 4) # cov of detection, followed by abundance

### CHECK MOD----
mod <- m_tmb
# gof(mod) not reliable
mod

mod@call$formula
plot_marginal(mod, "state")
plot_marginal(mod, "det")


saveRDS(mod, paste0("pcount_", park, "_", spec, ".RDS"))



### TREND ESTIMATE ---- 
# NOTE: IF THE ESTIMATE FOR YR_EFFECT IS SCALED, IT GIVES THE CHANGE PER 1 SD CHANGE IN YEAR. TO CONVERT TO CHANGE IN YEAR, DIVIDE BY THE SD OF THE YEAR COL. YOU CAN FIND THAT BY: sd(umf@siteCovs$yr). IN GENERAL, TO UNSCALE A SCALED VECTOR, MULTIPLY BY SD AND THEN ADD THE MEAN.
names(mod)
yr_effect <- ((rstan::extract(mod, "beta_state[yr_sc]")[[1]])/sd(umf@siteCovs$yr)) %>% exp(.) # the mean of this vector is the model estimate, then divide by year sd to unscale it, then exponentiate it
(trend_estim <- quantile(yr_effect, probs = c(0.025, 0.5, 0.975))) # this gives the same 95% CI that the model output shows
## I SHOULD BE ABLE TO CALCULATE THIS STRAIGHT FROM THE MODEL OUTPUT ALSO B/C 95%CI PROVIDED-- yes, confirmed that is true

# FOR SUBUNITS
names(mod)
yr_effect <- (rstan::extract(mod, "beta_state[yr_c]")[[1]]) %>% exp(.)

yr_effect <- (rstan::extract(mod, "beta_state[yr_c]")[[1]] + rstan::extract(mod, "beta_state[yr_c:subunitRA]")[[1]]) %>% exp(.)
(trend_estim <- quantile(yr_effect, probs = c(0.025, 0.5, 0.975))) 
 

### CHECK THE MODEL ----
# looking at the summary of best_mod, we conclude MCMC chains have converged if all R^<1.05.The rule of thumb is to have n_eff > 100 * number of chains (400). Otherwise, an easy solution would be to re-run this model with more iterations.
summary(mod, "state")
summary(mod, "det")

# Look at traceplots to check for convergence
ubms::traceplot(mod, pars=c("beta_state", "beta_det"))

### Model fit
# Plot residuals against fitted values. When the submodel has a binomial response (e.g., detection models), regular residual plots are not typically informative. Instead, the residuals and fitted values are divided into bins based on fitted value and the averages are plotted. For a count response (e.g., Poisson), Pearson residuals are calculated.
#  For binned residual plots, the shaded area represents plus/minus two standard deviations around the mean residual. If the model fits the data well, you would expect 95% of the binned residual points to fall within the shaded area.
ubms::plot_residuals(mod, submodel = "state")
ubms::plot_residuals(mod, submodel = "det")

# # Plot residuals against a covariate
plot_residuals(mod, submodel = "state", covariate = "yr_sc")
plot_residuals(mod, submodel = "state", covariate = "understory_cover_sd_50")

plot_residuals(mod, submodel = "det", covariate = "understory_cover_sd_200")
plot_residuals(mod, submodel = "det", covariate = "julian_prop_c")
plot_residuals(mod, submodel = "det", covariate = "hrs_since_rise_c")

# Instead of AIC, models are compared using leave-one-out cross-validation (LOO) (Vehtari, Gelman, and Gabry 2017) via the loo package. Based on this cross-validation, the expected predictive accuracy (elpd) for each model is calculated. The model with the largest elpd (fit_global) performed best. The elpd_diff column shows the difference in elpd between a model and the top model; if this difference is several times larger than the standard error of the difference (se_diff), we are confident the model with the larger elpd performed better. LOO model weights, analogous to AIC weights, are also calculated.

(mod_loo <- ubms::loo(mod, cores = 5))
plot_loo <- plot(mod_loo, label_points = TRUE)

## Posterior distributions
(plot_posterior <- ubms::plot_posteriors(mod)) # check to make sure approx. normally distributed

# Extracting individual parameters. Extract the posterior for the effect of a covariate
# names(best_mod)

### Check estimates
# In order to get the actual random intercept values, we use the ranef function.
# generate summary statistics for each random intercept. I think this gives the estimate for each level. That is, it adds the random effect to the intercept.
# ubms::ranef(mod, submodel="state", summary=TRUE)

# You can use the posterior_predict function to simulate new datasets, which you can use to calculate your own fit statistic. The following command generates 100 simulated datasets.
sim_y <- posterior_predict(mod, "y", draws=2000)
dim(sim_y)

# The output is a matrix with dimensions draws x observations (in site-major order). As an example, we can calculate the proportion of zeros in each simulated dataset

prop0 <- apply(sim_y, 1, function(x) mean(x==0, na.rm=TRUE))
# and compare that to the proportion of zeros in the actual dataset.
max_val <- apply(sim_y, 1, function(x) max(x, na.rm = TRUE))

actual_prop0 <- sum(getY(mod)==0, na.rm = T)/sum(!is.na(getY(mod)))
actual_max <- max(getY(mod), na.rm = TRUE)

hist0_df <- data.frame(prop0 = prop0)
(hist_prop0 <- ggplot(hist0_df, aes(x = prop0)) + 
  geom_histogram(bins = 12, fill = "white", color = "black") +
  labs(title = paste0(spec, "@", park, ": % of survey zero-counts, expected vs. actual (red line)"), subtitle = deparse(mod@call$formula, width.cutoff = 500L), x = "Proportion of zero counts in surveys", y = "# of model-generated samples") + 
  geom_vline(xintercept = actual_prop0, color = "red") + 
  theme_bw())

(max_plot <- ggplot(data.frame(max_val = max_val), aes(x = max_val)) + 
    geom_histogram(bins = 12, fill = "white", color = "black") +
    labs(title = paste0(spec, "@", park, ": Maximum count, expected vs. actual (red line)"), subtitle = deparse(mod@call$formula, width.cutoff = 500L), x = "Maximum count in surveys", y = "# of model-generated samples") + 
    geom_vline(xintercept = actual_max, color = "red") + 
    theme_bw())


### Model Inference
# Marginal covariate effects

## Plot marginal effects and 95% credible intervals of covariates
# THIS ONE IS USEFUL
# When creating marginal_effects for a particular predictor (or interaction of two predictors), one has to choose the values of all other predictors to condition on. By default, the median is used for continuous variables and the reference category is used for factors. Random effects are ignored
# This seems to pretty much give the same plots as plot_effects()
mod@call$formula
plot_marginal(mod, "state")
plot_marginal(mod, "det")

# (plot_marginal_abund <- plot_marginal(mod, "state") + labs(title = paste0(spec, "@", park, ": Predicted trend in abundance"), subtitle = deparse(mod@call$formula, width.cutoff = 500L), x = "Year (centered)", y = "Average abundance per plot"))
# 
# if("grob" %in% class(plot_marginal(mod, "det"))) {
#   tg <- ggpubr::text_grob(paste0(spec, "@", park, ": Predicted covariate effects on detection"), size =13)
#   sg <- ggpubr::text_grob(deparse(mod@call$formula, width.cutoff = 500L), size =10)
#   g = c(list(tg),list(sg),plot_marginal_cov)
#   (plot_marginal_cov <- gridExtra::grid.arrange(grobs=list(tg, sg, plot_marginal_cov),heights=c(0.5,0.5,10)))
# } else {
#   (plot_marginal_cov <- plot_marginal(mod, "det")) +
#     labs(title = paste0(spec, "@", park, ": Predicted covariate effects on detection"), subtitle = deparse(mod@call$formula, width.cutoff = 500L), x = "Year (centered)", y = "Average abundance per plot")
# }


## # Predict parameter values
# THIS IS NICE
# Get predicted mean abundance (lambda) estimates for each site-year (has same number of rows as original data frame)
 #>>> DOES ESTIMATED DETECTION PROBABILITY CHANGE YEAR-TO-YEAR? 
save_prefix <- paste(park, spec, sep = "_")
pred <- predict(mod, submodel = "state")
names(pred) <- c("plot_abund", "plot_abund_SD", "plot_abund_low", "plot_abund_high")
predict_dat <- cbind(
  mod@data@siteCovs,
  data.frame(obs_counts = rowSums(mod@data@y, na.rm = TRUE)),
  pred)

actual1 <- cbind(mod@data@siteCovs, obs_count = mod@data@y[,1], visit = 1)
actual2 <- cbind(mod@data@siteCovs, obs_count = mod@data@y[,2], visit = 2)
# actual_dat <- rbind(actual1, actual2) %>%
  dplyr::filter(!is.na(obs_count))
actual_dat <- actual1%>%
  dplyr::filter(!is.na(obs_count))

(plot_abund <- ggplot(predict_dat, aes(x = yr, y = plot_abund)) +
  labs(title = paste0("Est'd mean abundance/plot: ", save_prefix), subtitle = toString(mod@call$formula), x = "Year", y = "Est'd mean # of birds/plot") +
  geom_line() + 
  geom_point(size = 0.8) + 
  geom_ribbon(aes(x = yr, ymin = plot_abund_low, ymax = plot_abund_high, fill = hab_type_200_comb), alpha = 0.3) +  # <<<<<<<
  geom_jitter(data = actual_dat, aes(x = yr, y = obs_count), size = 0.8, color = "red") +
  theme_bw() + 
  facet_wrap("location_name"))
  plotly::ggplotly(plot_abund)
  











### COMPILE ALL RESULTS INTO DATAFRAME AND LIST----
# filenames <- list.files("~/NPS PROJECTS_CURRENT/GULN_birds/R/Model_fits/Stan_Pcount", pattern="*.RData", full.names=FALSE)
# 
# PCOUNT_elpd_list <- list()
# PCOUNT_models_list <- list()
# for (i in filenames) {

  i=filenames[13] # <<<<<<<<<<<<
  cat(i)
  load(paste0("~/NPS PROJECTS_CURRENT/GULN_birds/R/Model_fits/Stan_Pcount/", i))
  
  PCOUNT_models_list[[i]] <- list()
  
  txt <- str_split_1(gsub(".RData", "", i), "_")
  
  save_prefix <- paste(park, spec, sep = "_")
  
  ### Compare models ----
  # First we combine the models into a fitList. Then we generate a model selection table:

  mod_list <- listN(m_full, m_null, m1, m3, m4, m5, m6, m7, m8)
  
  
  
  (compare_mod <- round(modSel(fitList(mod_list)), 3))
  
  # Model formulas
  mod_form_df <- data.frame(mod_formula = rep(NA, nrow(compare_mod)))
  for(j in 1:nrow(compare_mod)) {
    mod_name <- rownames(compare_mod)[j]
    mod_form_df$mod_formula[j] <- get(mod_name)@call$formula %>% deparse(., width.cutoff = 500L)
  }
  
  # ELPD table
  elpd_df <- compare_mod %>%
    dplyr::select(-elpd, -nparam, -se_diff) %>%
    dplyr::rename(delta_aic_elpd = elpd_diff, mod_weight = weight) %>%
    dplyr::mutate(
      rdata_file = gsub(".RData", "", i),
      method = "STAN_PCOUNT",
      species_code = txt[2],
      unit_code = txt[3],
      dist_bound = "unlimited",
      mod_mixture = "P"
    )
  elpd_df$model_rank <- 1:nrow(elpd_df)

  temp_loo_df <- data.frame(mod_gof = as.character(), chi_pvalue = as.numeric(), prop0_boot_p = as.numeric(), plot_abund = as.numeric(), plot_abund_lcl = as.numeric(), plot_abund.ucl = as.numeric(), trend_signif = as.character())
  # Now build the rest of the data frame model by model
  for(j in 1:nrow(compare_mod)) {
    cat(j, "// ")
    mod <- get(rownames(compare_mod)[j])
    
    PCOUNT_models_list[[i]][[rownames(compare_mod)[j]]] <- mod
    
     
  # Build the rest of the summary data frame, model by model
    (mod_gof <- gof(mod, draws=1000, quiet=TRUE))
    (mod_loo <- ubms::loo(mod, cores = 5))
    bad_k <- sum(mod_loo$pointwise[,5]>0.7) # this is the number of bad pareto k
    bad_k_perc <- round(bad_k/dim(mod_loo$pointwise)[1]*100,1)
    
    sim_y <- posterior_predict(mod, "y", draws=2000)
    prop0 <- apply(sim_y, 1, function(x) mean(x==0, na.rm=TRUE)) # and compare that to the proportion of zeros in the actual dataset.
    actual_prop0 <- sum(getY(mod)==0, na.rm = T)/sum(!is.na(getY(mod)))
    hist0_df <- data.frame(prop0 = prop0)
    
      pred <- predict(mod, submodel = "state")
  names(pred) <- c("plot_abund", "plot_abund_SD", "plot_abund_low", "plot_abund_high")
  predict_dat <- cbind(
    mod@data@siteCovs,
    data.frame(obs_counts = rowSums(mod@data@y, na.rm = TRUE)),
    pred)
  
  yr_effect <- ((rstan::extract(mod, "beta_state[yr_sc]")[[1]])/sd(umf@siteCovs$yr)) %>% exp(.) # the mean of this vector is the model estimate, then divide by year sd to unscale it, then exponentiate it
  (trend_estim <- quantile(yr_effect, probs = c(0.025, 0.5, 0.975)))
    
  temp_loo_df[j, "mod_gof"] <- paste0("Number of data fits with Pareto k> 0.7 (bad) = ", bad_k, " (", bad_k_perc, "%)")
  temp_loo_df[j, "chi_pvalue"] <- mod_gof@post_pred_p
  temp_loo_df[j, 3:6] <- c(
    round(sum(actual_prop0 > hist0_df$prop0)/length( hist0_df$prop0), 2),
    mean(predict_dat$plot_abund, na.rm = TRUE) %>% round(., 2),
    mean(predict_dat$plot_abund_low, na.rm = TRUE) %>% round(., 2),
    mean(predict_dat$plot_abund_high, na.rm = TRUE) %>% round(., 2))
  temp_loo_df[j, "trend_signif"] <- 
    case_when(
        sum(trend_estim > 1)==3 ~ "UP", 
        sum(trend_estim > 1)==0 ~ "DOWN",
        .default = "not signif.")
  }# run through all models
      
    temp_loo <- cbind(temp_loo_df, mod_form_df, elpd_df) %>%
      dplyr::select(rdata_file, method, species_code, unit_code, dist_bound, model_rank, delta_aic_elpd, mod_weight, mod_mixture, mod_formula, everything()) %>%
      dplyr::mutate(across(where(is.double), \(x) round(x, 3)))
    
    PCOUNT_elpd_list[[i]] <- temp_loo
    rm(list = names(mod_list))
    rm(list = c("compare_mod", "elpd_df", "hist0_df", "mod", "mod_form_df", "mod_list", "pred", "predict_dat", "temp_loo_df", "umf", "mod_gof", "mod_loo", "temp_loo", "sim_y"))
}
PCOUNT_elpd_df <- do.call("rbind.data.frame", PCOUNT_elpd_list)

saveRDS(PCOUNT_elpd_df,  file = "~/NPS PROJECTS_CURRENT/GULN_birds/R/Model_fits/Stan_Pcount/PCOUNT_elpd_df.RDS")
write_csv(PCOUNT_elpd_df,  file = "~/NPS PROJECTS_CURRENT/GULN_birds/R/Model_fits/Stan_Pcount/PCOUNT_elpd_df.csv")
saveRDS(PCOUNT_models_list,  file = "~/NPS PROJECTS_CURRENT/GULN_birds/R/Model_fits/Stan_Pcount/PCOUNT_models_list.RDS")


# >>>>>>>>>>> END HERE 


# # # You can also supply newdata as a data.frame.
# nd <- data.frame(location_name = "TC01", hab_type_comb = "forest", yr_sc = 0, yr_f = "2017", julian_prop = 0.389, start_frac_time = 0.328, weather_wind_comb = "windy")
# predict(best_mod, submodel="state", newdata=nd)

# ## One of the advantages of BUGS/JAGS is that you can directly model latent parameters, such as the true unobserved occupancy state of a site z. Using the posterior_predict function in ubms, you can generate an equivalent posterior distribution of latent abundance z (niter x nsite-year)
# # NOTE: THESE ARE NOT AFFECTED BY SCALE
# 
# zpost <- ubms::posterior_predict(best_mod, "z", draws=100) # "y" means the observed outcome and "z" means the latent unobserved state
# hist(rowMeans(zpost, na.rm = TRUE))
# mean(rowMeans(zpost, na.rm = TRUE))
# 
# ypost <- ubms::posterior_predict(best_mod, "y", draws=100)
# y_means <- colMeans(ypost, na.rm = T)
# # The output has one row per posterior draw, and one column per site. The posterior of z can be useful for post-hoc analyses.
# 
# # For example, suppose you wanted to test for a difference in mean occupancy probability between sites 1-50 and sites 51-100:
# 
# group1 <- rowMeans(zpost[,1:50], na.rm=TRUE)
# group2 <- rowMeans(zpost[,51:100], na.rm=TRUE)
# 
# plot_dat <- rbind(
#   data.frame(
#     group="group1",
#     occ=mean(group1),
#     lower=quantile(group1, 0.025),
#     upper=quantile(group1, 0.975)),
#   data.frame(group="group2",
#              occ=mean(group2),
#              lower=quantile(group2, 0.025),
#              upper=quantile(group2, 0.975)))
# 
# # Now plot the posterior distributions of the two means:
# ggplot(plot_dat, aes(x=group, y=occ)) +
#   geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
#   geom_point(size=3) +
#   ylim(0,1) +
#   labs(x="Group", y="Occupancy + 95% UI") +
#   theme_bw() +
#   theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
#         axis.text=element_text(size=12), axis.title=element_text(size=14))
# 
# ## Look at the Stan code used to fit the model
# # Note that because a common Stan file is used to fit several types of models (to improve compilation times), there are many parts of the output that will not be used in a given model. Thus this code is much more complicated than what you would need to  fit your own version of an N-mixture model in Stan
# cat(get_stancode(best_mod))

  