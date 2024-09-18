### GULN BIRDS - BRMS

# rm(list = ls())

### Load packages ----
library(tidyverse)
library(brms)
library(magrittr)

### Functions ----
FuncLOOTable_BRM <- function(mod_list, AICwts_vec, stackwts_vec) {
# Function to generate a LOO summary table
#
# Args:
#   mod_list:  Named list of fitted regression models (with 'add_ic' already applied. Add prior to function so can check and reloo as necessary
#   AICwts_vec:  Pass in a named vector of the weights, as estimated by 'model_weights' with weights = "loo" (AIC weights computed based on LOO-IC)
#   stackwts_vec:  Pass in a named vector of the weights, as estimated by 'model_weights' with weights = "loo2" (stacked weights computed)
#
# Returns:
#   IC summary table with weights
#
temp_list <- lapply(mod_list, function(x) as.numeric(c(x$criteria$loo$estimates[c("looic"),], x$criteria$loo$estimates[c("p_loo"),])))
temp_df <- data.frame(do.call("rbind", temp_list))
temp_summary <- cbind(names(mod_list), temp_df) %>%
mutate_if(is.factor, as.character)
rownames(temp_summary) <- NULL
colnames(temp_summary) = c("model", "looic", "se(looic)", "p-loo", "se(p-loo)")
AICwts_df <- enframe(AICwts_vec)
stackwts_df <- enframe(stackwts_vec)
loo_summary <- temp_summary %>%
full_join(AICwts_df, by = c("model" = "name")) %>%
dplyr::rename(AICweight = value) %>%
full_join(stackwts_df, by = c("model" = "name")) %>%
dplyr::rename(stackingweight = value) %>%
arrange(looic) %>%
mutate_if(is.numeric, round, 2)
return(loo_summary)
}

### BRMS ANALYSIS ----
# spec <- "TUTI"



## BRMS models ----
# Simplest model to try
m0_brm <- brm(sum_indiv ~ yr_base*unit_code + (1|location_name),
              data = mod_dat_complete,
              family = poisson,
              iter = 4000, warmup = 1000, chains = 4, cores = 4,
              file = here::here("Model_fits", "TUTI_m0")
              )
m0_brm <- add_criterion(m0_brm, criterion = c("loo", "waic"), force_save = TRUE)

# ... habitat type
m1_brm <- update(m0_brm, formula. = ~ . + hab_type_200_comb, newdata = mod_dat_complete, file = here::here("Model_fits", "TUTI_m1"))
m1_brm <- add_criterion(m1_brm, criterion = c("loo", "waic"), force_save = TRUE)

# ... habitat type with PRIORS
m1_priors_brm <- brm(sum_indiv ~ yr_base*unit_code + hab_type_200_comb + (1|location_name),
                     data = mod_dat_complete,
                     family = poisson,
                     prior = c(prior(normal(0, 100), class = Intercept), 
                               prior(normal(0, 1), class = b)), 
                     iter = 4000, warmup = 1000, chains = 4, cores = 4,
                     file = here::here("Model_fits", "TUTI_m1_priors_brm")
)
m1_priors_brm <- add_criterion(m1_priors_brm, criterion = c("loo", "waic"), force_save = TRUE)



## Evaluate BRMS models ----
#### USER INPUT ON MODEL CHECKS BELOW... <<<<<<<<<<<<<<<<<<<<
mod_brms = trendmod_brms # <<<<<<<<<<<<< SELECT BEST MODEL
summary(mod_brms) # slightly shorter summary of model results

# Coefficient plot
mcmc_plot(mod_brms) # quick coefficient plot


# Model convergence check
plot(mod_brms) # plots showing posterior distribution of parameters, and trace plots


# Posterior predictive checks
# pp_check(mod_brms, type = "dens_overlay", ndraws = 200)
pp_check(mod_brms, type = "stat_2d", ndraws = 200)
# pp_check(mod_brms, type = "loo_pit_qq", ndraws = 200)
pp_check(mod_brms, type = "rootogram", ndraws = 200)
pp_check(mod_brms, type = "bars_grouped", group = "hab_type_50_comb") # <<<<<<<< CHANGE, DEPENDING ON BEST MODEL
pp_check(mod_brms, type = "violin_grouped", group = "hab_type_100_comb") # <<<<<<<< CHANGE, DEPENDING ON BEST MODEL

## Plot conditional effects
plot(conditional_effects(mod_brms, effects = "hrs_since_rise_c"), points=TRUE, point_args = list(width = 0.3, height = 0.1, alpha = 0.4))

# Save conditional effects plots as list, for use with ggplot functions
p <- plot(conditional_effects(mod_brms), points=TRUE, point_args = list(width = 0.3, height = 0.1, alpha = 0.4), plot = FALSE)
p$unit_code + theme_bw(base_size = 14)
p$`yr_base:unit_code` + theme_bw(base_size = 14) + facet_wrap(vars(unit_code))

## Look for high leverage points
plot(loo(mod_brms), label_points = TRUE)

## Select among BRMS models ----
# To get a table with model weights (assumes you've added the IC to the model, as shown above)...
mod_list <- list(m0_brm = m0_brm, m1_priors_brm = m1_priors_brm, m3_brm = m3_brm, m6_brm = m6_brm)
AICwts_vec = model_weights(m0_brm, m1_priors_brm, m3_brm, m6_brm,  weights = "loo") # assumes you've used 'add_ic' to append LOO-IC to each model
stackwts_vec = model_weights(m0_brm, m1_priors_brm, m3_brm, m6_brm, weights = "stacking")

# mod_list <- list(m0_brm = m0_brm, m1_brm = m1_brm, m1_priors_brm = m1_priors_brm, m2_brm = m2_brm, m3_brm = m3_brm, m4_brm = m4_brm, m5_brm = m5_brm, m6_brm = m6_brm, m7_brm = m7_brm, m8_brm = m8_brm)
# AICwts_vec = model_weights(m0_brm, m1_brm, m1_priors_brm, m2_brm, m3_brm, m4_brm, m5_brm, m6_brm, m7_brm, m8_brm,  weights = "loo") # assumes you've used 'add_ic' to append LOO-IC to each model
# stackwts_vec = model_weights(m0_brm, m1_brm, m1_priors_brm, m2_brm, m3_brm, m4_brm, m5_brm, m6_brm, m7_brm, m8_brm, weights = "stacking")

FuncLOOTable_BRM(mod_list = mod_list, AICwts_vec = AICwts_vec, stackwts_vec = stackwts_vec)

## BRMS UPSHOT ----
# Best-supported model included hab_type_200_comb (forest vs. non-forest), next supported also added noise but it's far less supported. This is consistent with classification tree results
# <<<<<<<<<<<<<< PICK UP FROM HERE--BRMS PERFORMS BETTER WITH THE PRIORS, BUT NEED TO CENTER/SCALE EVERYTHING AND SELECT WEAKLY INFORMATIVE PRIORS. INCLUDE TEMPERATURE AS A POLYNOMIAL.

