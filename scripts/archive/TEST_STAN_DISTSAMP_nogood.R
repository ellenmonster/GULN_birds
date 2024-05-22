### N-MIXTURE PLUS DISTANCE SAMPLING
# https://groups.google.com/g/unmarked/c/OesoTLrawgE

### NOTE!!! Only using a subset of data that are used in the other analyses! Only using the data collected within the first 100m!
### >>>>>>>>>>>> JELA NOCA & WEVI DATA--WITH EACH YEAR THE % OF DATA IN THE OUTERMOST BIN INCREASES!! <<<<<<<<<<<<< possibly paying attention further out with time, or distance estimation is increasingly off, or they are becoming more liberal with attributing it to the last distance class. PAAL NOCA SEEMS TO HAVE SIMILAR ISSUE, NOT AS BAD. VICK is good. SAAN has a high% in last bin for first year for both NOCA & WEVI. BITH, GUIS, no obvious issue

### WHAT I AM DOING HERE--ALWAYS CUT OFF THE DATA AT 100M. FOR SPECIES THAT REQUIRE COMBINING BINS OR THAT OTHERWISE CAN'T BE USED WITH STANDARD DISTANCE ANALYSIS, SKIP THEM.

rm(list = ls())

### Load packages ----
library(tidyverse)
library(magrittr)
library(ubms)

source(here::here("scripts", "TEST_STARTUP.R"))
  
### Load files----
df_finaldat <- readRDS(here::here("Data_out", "df_finaldat.RDS"))

### USER ENTRY ----
park <- "VICK"
spec <- "NOCA"

max_bin_bound <- 100 # <<<<<< MUST enter some value--cannot use Inf
combine_bins = FALSE # TRUE to combine first and second bin intervals

### Format the data ----
dat <- FuncFormatFinalDat(park = park, spec = spec)

# # CHECK THE ISSUE WITH INCREASING % OF COUNTS IN LAST BIN!!!
# x <- dat %>%
#   dplyr::mutate(
#     distbin_0to2 = distbin_0 + distbin_1 + distbin_2
#     ) %>%
#   dplyr::select(yr, distbin_0to2, distbin_3) %>%
#   dplyr::group_by(yr) %>%
#   dplyr::summarize(
#     sum_0to2bins = sum(distbin_0to2),
#     sum_3bin = sum(distbin_3)
#   ) %>%
#   dplyr::mutate(prop_last_bin = sum_3bin/(sum_3bin+sum_0to2bins))
# x
# plot(x$prop_last_bin~x$yr)
table(dat$researcher, dat$yr, dat$subunit) # check how many observers

# Dataset with complete cases only
# PAAL: Dropped 3 records b/c missing weather_wind_comb
# SAAN: Dropped 1 record b/c missing hrs_since_rise
nrow(dat)
colSums(is.na(dat))

dat <- dat[complete.cases(dat),]
nrow(dat)

dat_Nmix <- dat

### N-MIXTURE ANALYSIS ----
### Format data for 'unmarked' and 'ubms' ----

### Response variable
y_st <- dplyr::select(dat_Nmix, starts_with("distbin")) %>% as.matrix()

dim(y_st)
head(y_st)

### Abundance and detection covariates <<<<<<<<< MAY NEED TO ADD RESEARCHER & SUBSET, DEPENDING ON PARK
covs_st <- dat_Nmix %>%
  dplyr::select(yr, yr_sc, yr_f, location_name, researcher, first_yr, julian_prop_c, hrs_since_rise_c, weather_wind_comb,  hab_type_100_comb, hab_type_200_comb, prop_understory_cover_50, prop_understory_cover_100, prop_understory_cover_200, understory_cover_sd_50, understory_cover_sd_100, understory_cover_sd_200)

covs_st$hrs_since_rise_sc <- scale(covs_st$hrs_since_rise_c, center = FALSE, scale = TRUE)
summary(covs_st)

### Create frame for stan_distsamp----
umf <- unmarkedFrameDS(
  y=y_st,
  siteCovs=covs_st,
  dist.breaks=c(0,.025,.050,.1),
  unitsIn="km", 
  survey="point")
summary(umf)


### TREND MODEL ----
m_tmb <- stan_distsamp( # detection covariates first
  ~hrs_since_rise_c + julian_prop_c + weather_wind_comb + hab_type_200_comb + prop_understory_cover_200 +  I(prop_understory_cover_200^2) + understory_cover_sd_50 +  I(understory_cover_sd_50^2)
  ~yr_sc  + (1|location_name) + (1|yr_f) + hab_type_200_comb + prop_understory_cover_200 +  I(prop_understory_cover_200^2),
  data = umf, 
  keyfun = "halfnorm",
  output = "abund",
  unitsOut = "ha",
  chains=4, 
  iter=4000, 
  thin = 2,
  cores = 4)

mod <- m_tmb
saveRDS(mod, here::here("stan_distsamp_mods", paste0("distsamp_", park, "_", spec, ".RDS")))

### YEARMOD ----
# Since the detection covariates are only for determining shape of distance-detection curve, and not for estimating availability (prob. of calling), julian_prop and hrs_since_rise should be covariates for ABUNDANCE, not for the DETECTION CURVE--unless we think that hrs_since_rise affects visual detections (which should be true, but only for hrs while still dark and supposedly most detections are acoustic)
vick_noca <- stan_distsamp( # detection covariates first, then abundance covariates
  formula = ~hrs_since_rise_sc + hab_type_200_comb + prop_understory_cover_200 + 
    understory_cover_sd_50 + I(understory_cover_sd_50^2) ~ yr_f + (1 | location_name) + hab_type_200_comb + prop_understory_cover_200 + I(prop_understory_cover_200^2),
  data = umf, 
  keyfun = "halfnorm",
  output = "abund",
  chains=4, 
  iter=6000, 
  thin = 3,
  cores = 20)

mod <- vick_noca
saveRDS(mod, here::here("stan_distsamp_mods", paste0("distsamp_", park, "_", spec, "yearmod.RDS")))


### CHECK MOD----
# looking at the summary of best_mod, we conclude MCMC chains have converged if all R^<1.05.The rule of thumb is to have n_eff > 100 * number of chains (400). Otherwise, an easy solution would be to re-run this model with more iterations.

# gof(mod) # This is not reliable
summary(mod, "state")

mod@call$formula
plot_marginal(mod, "det")
plot_marginal(mod, "state")

### TREND ESTIMATE ---- 
# NOTE: IF THE ESTIMATE FOR YR_EFFECT IS SCALED, IT GIVES THE CHANGE PER 1 SD CHANGE IN YEAR. TO CONVERT TO CHANGE IN YEAR, DIVIDE BY THE SD OF THE YEAR COL. YOU CAN FIND THAT BY: sd(umf@siteCovs$yr). IN GENERAL, TO UNSCALE A SCALED VECTOR, MULTIPLY BY SD AND THEN ADD THE MEAN.
names(mod)
yr_effect <- ((rstan::extract(mod, "beta_state[yr_sc]")[[1]])/sd(umf@siteCovs$yr)) %>% exp(.) # the mean of this vector is the model estimate, then divide by year sd to unscale it, then exponentiate it
(trend_estim <- quantile(yr_effect, probs = c(0.025, 0.5, 0.975))) # this gives the same 95% CI that the model output shows

# # FOR SUBUNITS
# names(mod)
# yr_effect <- (rstan::extract(mod, "beta_state[yr_c]")[[1]]) %>% exp(.)
# 
# yr_effect <- (rstan::extract(mod, "beta_state[yr_c]")[[1]] + rstan::extract(mod, "beta_state[yr_c:subunitGUIS-MS]")[[1]]) %>% exp(.)
# (trend_estim <- quantile(yr_effect, probs = c(0.025, 0.5, 0.975))) 

### CHECK THE MODEL ----

# Look at traceplots to check for convergence
ubms::traceplot(mod, pars=c("beta_state", "beta_det"))

### Model fit
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

(mod_loo <- ubms::loo(mod, cores = 20))
plot_loo <- plot(mod_loo, label_points = TRUE)

## Posterior distributions
(plot_posterior <- ubms::plot_posteriors(mod)) # check to make sure approx. normally distributed

# Extracting individual parameters. Extract the posterior for the effect of a covariate
# names(best_mod)

### Check estimates
# generate summary statistics for each random intercept. I think this gives the estimate for each level. That is, it adds the random effect to the intercept.
# ubms::ranef(mod, submodel="state", summary=TRUE)

# You can use the posterior_predict function to simulate new datasets, which you can use to calculate your own fit statistic. The following command generates 100 simulated datasets.
sim_y <- posterior_predict(mod, draws=2000)
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
