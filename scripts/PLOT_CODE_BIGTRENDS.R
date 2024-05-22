library(tidyverse)
# library(ggeffects)
# library(ubms)
library(patchwork)
library(emmeans)
library(bbsBayes2)

### GULN PARK TRENDS ----
### Read in files 
tmb_files <- list.files(path = here::here("tmb_mods"), full.names = FALSE)

trends_df <- data.frame(
  source = as.character(), 
  park = as.character(), 
  subunit = as.character(),
  spec = as.character(), 
  annual_perc_change = as.numeric(), 
  conf_low = as.numeric(), 
  conf_high = as.numeric(), 
  signif = as.numeric())

for(i in 1:length(tmb_files)) {
  mod <- readRDS(here::here("tmb_mods", tmb_files[i]))
  
  cat(tmb_files[i])
  
  spl <- stringr::str_split(gsub(".RDS", "", tmb_files[i]), pattern = "_")
  park <- spl[[1]][2]
  spec <- spl[[1]][3]
  
  ## Get the data
  dat <- FuncFormatFullObs(park = park, spec = spec)
  dat <- dat[complete.cases(dat),]
  # Scale predictors
  dat$prop_understory_cover_50_sc <- scale(dat$prop_understory_cover_50)
  dat$prop_understory_cover_100_sc <- scale(dat$prop_understory_cover_100)
  dat$prop_understory_cover_200_sc <- scale(dat$prop_understory_cover_200)
  dat$understory_cover_sd_50_sc <- scale(dat$understory_cover_sd_50)
  dat$understory_cover_sd_100_sc <- scale(dat$understory_cover_sd_100)
  dat$understory_cover_sd_200_sc <- scale(dat$understory_cover_sd_200)
  
  ## Get the trend
  em=emmeans::emtrends(mod, specs = ifelse(park %in% c("GUIS", "SAAN") & "subunit" %in% names(mod$frame), c("subunit"), c("yr_c")), var = c("yr_c"), data = dat)
  trend <- summary(em, infer=c(TRUE,TRUE),null=0) %>%
    as.data.frame() %>%
    janitor::clean_names() %>%
    dplyr::mutate(signif = p_value < 0.05)
  
  if(!park %in% c("GUIS", "SAAN")) {trend$subunit = park}
  if(park %in% c("GUIS", "SAAN") & !"subunit" %in% names(mod$frame)) {
    temp <- dplyr::left_join(mod$frame, dat[, c("subunit", "location_name")] %>% distinct(), by= "location_name")
    trend$subunit <- unique(temp$subunit)
  }
  trend %<>%
    dplyr::select(subunit, yr_c_trend, asymp_lcl, asymp_ucl, signif) %>%
    dplyr::mutate(
      source = "GULN park unit", 
      park = park, 
      spec = spec,
      annual_perc_change = 100*yr_c_trend,
      conf_low = 100*asymp_lcl,
      conf_high = 100*asymp_ucl) %>%
    dplyr::select(source, park, subunit, spec, annual_perc_change, conf_low, conf_high, signif)
  
  # Get mean fitted for each subunit-spec
  fit <- mean(fitted(mod), na.rm = TRUE)
  fit_df <- cbind(dat[, c("subunit")], fit)
  subunit_mean <- fit_df %>%
    dplyr::group_by(subunit) %>%
    dplyr::summarize(estim_mean_per_surv = mean(fit))
  
  # Add the mean to trend
  trend %<>% dplyr::left_join(subunit_mean, by = "subunit")
  
  trends_df <- rbind(trends_df, trend)
}
trends_df$subunit[trends_df$subunit=="MI"] <- "SAAN-MI"
trends_df$subunit[trends_df$subunit=="RA"] <- "SAAN-RA"

View(trends_df)
saveRDS(trends_df, here::here("RESULTS", "GULN_TREND_RESULTS.RDS"))

### GET BBS RESULTS----
bbs_files <- list.files(path = here::here("bbs_mods"), pattern = "_spatial_fattail", full.names = FALSE)

bbs_trends_df <- data.frame(
  source = as.character(), 
  park = as.character(), 
  subunit = as.character(),
  spec = as.character(), 
  annual_perc_change = as.numeric(), 
  conf_low = as.numeric(), 
  conf_high = as.numeric(), 
  signif = as.numeric())

for(i in 1:length(bbs_files)) {
  mod <- readRDS(here::here("bbs_mods", bbs_files[i]))
  
  spl <- stringr::str_split(bbs_files[i], pattern = "_")
  spec <- spl[[1]][3]
  
  # Get trends
  ind <- generate_indices(model_output = mod, hpdi = TRUE) # Calculate annual indices of relative abundance by year for different regions.
  
  (plot_ind <- plot_indices(indices = ind, 
                            add_observed_means = TRUE, # optional argument to show raw observed mean
                            add_number_routes = TRUE))
  
  temp_trend <- generate_trends(ind, quantiles = c(0.025, 0.05, 0.95, 0.975), slope = TRUE, hpdi = TRUE, prob_decrease = c(0,25))
  trend_est <- temp_trend$trends %>%
    dplyr::filter(region_type == "stratum") %>%
    dplyr::mutate(
      source = "BBS BCR",
      spec = spec,
      signif = sign(slope_trend_q_0.025) == sign(slope_trend_q_0.975)) %>%
    dplyr::select(source, park = strata_included, subunit = strata_included, spec, annual_perc_change = slope_trend, conf_low = slope_trend_q_0.025, conf_high = slope_trend_q_0.975, signif)
  
  bbs_trends_df <- rbind(bbs_trends_df, trend_est)
  }
View(bbs_trends_df)
bbs_trends_df$estim_mean_per_surv <- NA
saveRDS(bbs_trends_df, here::here("RESULTS", "BBS_TREND_RESULTS.RDS"))


all_trends <- rbind(trends_df, bbs_trends_df)

name_template <- data.frame(rbind(
  c("AMCR", "American crow"),
  c("BLJA", "Blue jay"),
  c("CARW", "Carolina wren"),
  c("EATO", "Eastern towhee"),
  c("GTGR", "Great-tailed grackle"),
  c("MODO", "Mourning dove"),
  c("NOCA", "Northern cardinal"),
  c("NOMO", "Northern mockingbird"),
  c("RBWO", "Red-bellied woodpecker"),
  c("TUTI", "Tufted titmouse"),
  c("WEVI", "White-eyed vireo"),
  c("YBCU", "Yellow-billed cuckoo")))
names(name_template) <- c("spec", "species_fullname")
all_trends %<>% dplyr::left_join(name_template, by = "spec")

all_trends$Trend <- "Not Signif."
all_trends$Trend[all_trends$signif == TRUE & all_trends$annual_perc_change> 0] <- "Signif. Increase"
all_trends$Trend[all_trends$signif == TRUE & all_trends$annual_perc_change< 0] <- "Signif. Decline"
all_trends$Trend <- factor(all_trends$Trend, levels = c("Signif. Increase", "Not Signif.", "Signif. Decline"))


all_trends %<>%
  dplyr::group_by(spec) %>% 
  mutate(rank = rank(-annual_perc_change, ties.method = "first")) %>%
  dplyr::arrange(spec, rank)


saveRDS(all_trends, here::here("RESULTS", "COMBINED_TREND_RESULTS.RDS"))

### CREATE THE BIG SUMMARY PLOT ----
all_trends <- readRDS(here::here("RESULTS", "COMBINED_TREND_RESULTS.RDS"))
all_trends$estim_mean_per_surv[all_trends$source == "BBS BCR"] <- 1.2
all_trends$subunit <- factor(all_trends$subunit, levels = c("BCR20", "BCR36", "SAAN-MI", "SAAN-RA", "BCR21", "PAAL", "BCR37", "BITH", "BCR25", "VICK", "JELA", "BCR26", "GUIS-MS", "GUIS-FL", "BCR27"))
all_trends %<>% dplyr::filter(estim_mean_per_surv > 0.2)


plot_dat <- all_trends %>% dplyr::filter(spec %in% group2) # #<<<<<<<<<<<
p<-ggplot(plot_dat) +
  geom_point(aes(x = subunit, y = annual_perc_change, size = estim_mean_per_surv, fill = Trend, shape = source, color = Trend)) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high, x = subunit, color = Trend), width = 0, size = 0.65) + 
  scale_shape_manual(values = c("GULN park unit" = 16, "BBS BCR" = 2)) +
  scale_fill_manual(values = c("Signif. Increase" = "#2166ac", "Not Signif." = "darkgray", "Signif. Decline" = "#b2182b")) +
  scale_color_manual(values = c("Signif. Increase" = "#2166ac", "Not Signif." = "darkgray", "Signif. Decline" = "#b2182b")) +
  scale_size_area(name = "# of birds/point count", max_size =5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_text(data = subset(plot_dat, subunit %in% c("SAAN-MI", "BITH", "JELA", "GUIS-FL")), aes(label = subunit, x = subunit, y = conf_high), fontface = "bold", vjust = -1.2, hjust = .5, size = 3) + # label above
  geom_text(data = subset(plot_dat, subunit %in% c("SAAN-RA", "PAAL", "VICK", "GUIS-MS")), aes(label = subunit, x = subunit, y = conf_low), fontface = "bold", vjust = 1.4, hjust = .5, size = 3) + # label below
   # labs(title = "White-eyed vireo trends in GULN park units and surrounding areas", subtitle = "Parks and BBS BCR's are ordered west (left) to east. For GULN, area of circle is scaled to mean # of birds/point count\nShowing park units with >0.2 birds/point count, averaged over all years", x = "Estimated annual % change (bars are 95%CI)", y = "Estimated % change per year") +
  labs(title = "Species with MOSTLY non-significant trends in GULN park units", subtitle = "Parks and BCR's are ordered west (left) to east. For GULN, area of circle is scaled to mean # of birds/point count\nShowing park units with >0.2 birds/point count, averaged over all years", x = "Estimated annual % change (bars are 95%CI)", y = "Estimated % change per year") +
  scale_y_continuous(expand = expansion(add = c(5, 7))) +
  scale_x_discrete(expand = expansion(add = c(1, 1))) +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 13),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
  facet_wrap(~species_fullname, nrow = 2)

p
levels = c("BCR20", "BCR36", "PAAL", "SAAN", "BCR21", "BCR37", "BITH", "BCR25", "BCR26", "VICK", "JELA", "BCR27", "GUIS")

# Need to set the relative abundacnce scale!
group1 <- c("AMCR", "RBWO", "TUTI")
group2 <- c("MODO", "NOMO", "BLJA", "CARW")


### >>>>>>> PLOT TRENDS FOR A SPECIFIC BCR <<<<<
### FULL PLOT ----

### First, read in the model as 'mod' ###
spec_name = "American crow"
ind <- generate_indices(model_output = mod, hpdi = TRUE) # Calculate annual indices of relative abundance by year for different regions.

trend_results <- ind$indices %>%
  dplyr::filter(region_type == "stratum") %>%
  dplyr::mutate(stratum_label = paste0(strata_included, " (", n_routes_total, " routes)"))

## Put all plots on a page
(p_all <- ggplot(trend_results) + 
 geom_errorbar( # 
   data = trend_results,
   mapping = aes(x = year, y = index, ymin = index_q_0.025, ymax = index_q_0.975), width = 0, color = "blue") +
   geom_point(
     data = trend_results,
     mapping = aes(x = year, y = index), size = 2.5, color = "blue") +
   scale_x_continuous(breaks = scales::breaks_pretty()) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
   labs(title = paste0("Breeding bird survey results for ", spec_name), x = "Year", y = "Index of abundance (mean/route)", subtitle = "Indices are estimated using overdispersed Poisson mixed model\nEach plot is a Bird Conservation Region. Y-axes vary across plots.") + 
    theme_bw()) + 
  facet_wrap(~stratum_label, scales = "free_y")

## Just generate a plot for one stratum
stratum <- "BCR20" # <<<<<<<<<<<<<<<

trend_sub <- trend_results %>% dplyr::filter(strata_included == stratum)

(p_one <- ggplot(trend_sub) + 
    geom_errorbar( # 
      data = trend_sub,
      mapping = aes(x = year, y = index, ymin = index_q_0.025, ymax = index_q_0.975), width = 0, color = "blue") +
    geom_point(
      data = trend_sub,
      mapping = aes(x = year, y = index), size = 2.5, color = "blue") +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    labs(title = paste0("Breeding bird survey results for ", spec_name, " (", stratum, ")"), x = "Year", y = "Index of abundance (mean/route)", subtitle = "Indices are estimated using overdispersed Poisson mixed model") + 
    theme_bw())
