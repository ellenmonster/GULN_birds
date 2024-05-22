### SAVE THIS TO WORKING COMPUTER

### GDISTSAMP DETECTION FUNCTION PLOTS
bin_areas_m2 <- c(
  bin0 = pi*25^2,
  bin1 = (pi*50^2) - pi*(25^2),
  bin2 = (pi*100^2) - pi*(50^2)
)
bin_areas_ha = bin_areas_m2/10000

### EXAMINE DETECTION FUNCTION PLOTS FOR SINGLE VISIT, NO COVS ----
bin_sums <- colSums(umf_gds@y, na.rm = TRUE)

adj_factor <- .75 # <<<<<<<<<<< EYEBALL THIS
count_per_ha = bin_sums/bin_areas_ha
count_dens_df <- data.frame(
  y = count_per_ha,
  y_scaled = (count_per_ha/count_per_ha[1])*adj_factor,
  x = c(median(c(0, 0.025)), median(c(0.025, 0.05)), median(c(0.05, 0.1)))
)

# Generate model-estimated detection function
# Get model estimates
det_coef <- mod@estimates[2]@estimates

# If there are covariates on detection, modify as necessary
(det_shape <- det_coef[1] %>% exp())

# HALF-NORM
plot(function(x) unmarked::gxhn(x, det_shape), 0, 0.1, main = paste0(park, " ", spec, " HALF-NORMAL DISTANCE FUNCTION"), xlab="Distance (km)", ylab="Detection probability", lwd = 1.5, cex.lab=1, cex.axis=1)
points(count_dens_df$y_scaled ~count_dens_df$x, pch = 16, col = "blue")

# EXP
plot(function(x) unmarked::gxexp(x, rate = det_shape), 0, 0.1, main = paste0(park, " ", spec, " EXPONENTIAL DISTANCE FUNCTION"), xlab="Distance (km)", ylab="Detection probability", lwd = 1.5, cex.lab=1, cex.axis=1)
points(count_dens_df$y_scaled ~count_dens_df$x, pch = 16, col = "blue")

# HAZARD
(det_scale <- mod@estimates[3]@estimates %>% exp())
plot(function(x) unmarked::gxhaz(x, shape = det_shape, scale = det_scale), 0, 0.1, main = paste0(park, " ", spec, " HAZARD DISTANCE FUNCTION"), xlab="Distance (km)", ylab="Detection probability", lwd = 1.5, cex.lab=1, cex.axis=1)
points(count_dens_df$y_scaled ~count_dens_df$x, pch = 16, col = "blue")

### EXAMINE DETECTION FUNCTION PLOTS FOR SINGLE VISIT, BY A COV ----
plot_dat <- cbind(umf_gds@y, umf_gds@yearlySiteCovs["researcher"])
# plot_dat <- cbind(umf_gds@y, umf_gds@siteCovs[c("hab_type_200_comb", "prop_understory_cover_200")])
glimpse(plot_dat)
names(plot_dat) <- c("bin0", "bin1", "bin2", "researcher")
(bin_sums <- plot_dat %>%
  dplyr::filter(!is.na(researcher)) %>%
  dplyr::group_by(researcher) %>%
  dplyr::summarize(
    bin0_sum = sum(bin0, na.rm = TRUE),
    bin1_sum = sum(bin1, na.rm = TRUE),
    bin2_sum = sum(bin2, na.rm = TRUE)
  ))

count_per_ha <- bin_sums %>%
  dplyr::mutate(
    bin0_dens = bin0_sum/bin_areas_ha[1],
    bin1_dens = bin1_sum/bin_areas_ha[2],
    bin2_dens = bin2_sum/bin_areas_ha[3]
  )

# Legend text
legend_vals <- count_per_ha %>%
  dplyr::select(researcher, bin0_sum, bin1_sum, bin2_sum) %>%
  dplyr::mutate(
    N = bin0_sum + bin1_sum + bin2_sum,
    plot_legend = paste0(researcher, " (N=", N, ")")
  ) %>%
  pull(plot_legend)

# Plot data
adj_factor <- 0.95 # <<<<<<<<<<< EYEBALL THIS

count_dens_df <- list()
for(i in 1:length(count_per_ha$researcher)) {
  count_dens_df[[count_per_ha$researcher[i]]] <- 
    data.frame(
      y = as.numeric(count_per_ha[i, c("bin0_dens", "bin1_dens", "bin2_dens")]),
      y_scaled = (as.numeric(count_per_ha[i, c("bin0_dens", "bin1_dens", "bin2_dens")])/max(as.numeric(count_per_ha[i, c("bin0_dens", "bin1_dens", "bin2_dens")])))*adj_factor, 
      x = c(median(c(0, 0.025)), median(c(0.025, 0.05)), median(c(0.05, 0.1)))
    )
}

# Get model estimates
(det_coef <- mod@estimates[2]@estimates)

# If there are covariates on detection, modify as necessary
# # For JELA WEVI only...
# det_shape <- exp(-4.0190744 + 0.354*2.0670017)
# det_shape2 <- exp(-4.0190744 + 0.5239949 + 0.238*2.0670017)

det_shape <- det_coef[1] %>% exp()
det_shape2 <- (det_coef[1] + det_coef[2]) %>% exp()
det_shape3 <- (det_coef[1] + det_coef[3]) %>% exp()

# HALF-NORMAL
plot(function(x) unmarked::gxhn(x, sigma = det_shape), 0, 0.1, main = paste0(park, " ", spec, " HALF-NORMAL DISTANCE FUNCTION"), xlab="Distance (km)", ylab="Detection probability", lwd = 1.5, cex.lab=1, cex.axis=1)
plot(function(x) unmarked::gxhn(x, sigma = det_shape2), 0, 0.1, col = "orange", lwd = 1.5, add = TRUE)
plot(function(x) unmarked::gxhn(x, sigma = det_shape3), 0, 0.1, col = "green", lwd = 1.5, add = TRUE)

# EXP
plot(function(x) unmarked::gxexp(x, rate = det_shape), 0, 0.1, main = paste0(park, " ", spec, " EXPONENTIAL DISTANCE FUNCTION"), xlab="Distance (km)", ylab="Detection probability", lwd = 1.5, cex.lab=1, cex.axis=1)
plot(function(x) unmarked::gxexp(x, rate = det_shape2), 0, 0.1, col = "orange", lwd = 1.5, add = TRUE)
plot(function(x) unmarked::gxexp(x, rate = det_shape3), 0, 0.1, col = "green", lwd = 1.5, add = TRUE)

# ...THEN ADD THE OBSERVED DATA
points(count_dens_df[[1]]$y_scaled ~count_dens_df[[1]]$x, pch = 16, col = "black")
points(count_dens_df[[2]]$y_scaled ~count_dens_df[[2]]$x, pch = 16, col = "orange")
points(count_dens_df[[3]]$y_scaled ~count_dens_df[[3]]$x, pch = 16, col = "green")
legend('topright', legend_vals, col=c("black", "orange", "green"), bty = "n", lty=1, lwd = 1.5, cex=1)

# SAVE AS JPEG WIDTH 656, HEIGHT 488 (e.g., gdistsamp_BITH_NOCA_1visit_yearmod_detectfunc_adj95) ----






# # Hazard----
# if(gdist_best_mod@keyfun == "hazard") {
#   det_scale <- gdist_best_mod@estimates[4]@estimates %>% exp()
#   
#   (gdist_detect_plot <- plot(function(x) unmarked::gxhaz(x, shape = det_shape, scale = det_scale), 0, max_bin_bound/1000, main = paste0("HAZARD Function: ", save_prefix), xlab="Distance (km)", ylab="Detection probability", cex.lab=0.7, cex.axis=0.7, las=1))
#   
#   # Hazard with COVARIATE (SAVES AS FUNCTION)
#   gdist_detect_plot <- function() {
#     det_shape <- 0.0007122781
#     det_shape2 <- 0.02675571 
#     det_shape3 <- 0.001653192 
#     det_shape4 <- 0.004203732 
#     det_shape5 <- 0.005407489 
#     det_scale <- 0.9352508  
#     plot(function(x) unmarked::gxhaz(x, shape = det_shape, scale = det_scale), 0, max_bin_bound/1000, main = paste0("HAZARD Function: ", save_prefix), xlab="Distance (km)", ylab="Detection probability", cex.lab=0.7, cex.axis=0.7, las=1)
#     plot(function(x) unmarked::gxhaz(x, shape = det_shape2, scale = det_scale), 0, max_bin_bound/1000, col = "green", add = TRUE)
#     plot(function(x) unmarked::gxhaz(x, shape = det_shape3, scale = det_scale), 0, max_bin_bound/1000, col = "orange", add = TRUE)
#     plot(function(x) unmarked::gxhaz(x, shape = det_shape4, scale = det_scale), 0, max_bin_bound/1000, col = "brown", add = TRUE)
#     plot(function(x) unmarked::gxhaz(x, shape = det_shape5, scale = det_scale), 0, max_bin_bound/1000, col = "red", add = TRUE)
#     
#     legend('topright', c("Developed", "Forest", "Mixed", "Open", "Shrub"), col=c("black", "green", "orange", "brown", "red"), lty=1, cex=1)
#   }
#   gdist_detect_plot()
# }
# 
# # Exponential ----
# if(gdist_best_mod@keyfun == "exp") {
# 
# (gdist_detect_plot <- plot(function(x) unmarked::gxexp(x, rate = det_shape), 0, max_bin_bound/1000, main = paste0("NEG. EXP. Function: ", save_prefix), xlab="Distance (km)", ylab="Detection probability", cex.lab=0.7, cex.axis=0.7, las=1))
# 
# # Exp with COVARIATE (SAVES AS FUNCTION)
# gdist_detect_plot <- function() {
#   det_shape <-0.03170303
#   det_shape2 <-0.01237553 
#   plot(function(x) unmarked::gxexp(x, rate = det_shape), 0, max_bin_bound/1000, main = paste0("NEG. EXP. Function: ", save_prefix), xlab="Distance (km)", ylab="Detection probability", cex.lab=0.7, cex.axis=0.7, las=1)
#   plot(function(x) unmarked::gxexp(x, rate = det_shape2), 0, max_bin_bound/1000, col = "green", add = TRUE)
#   legend('topright', c("Calm", "Windy"), col=c("black", "green"), lty=1, cex=1)
# }
# }