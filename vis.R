library(coda)
library(nimble)
library(tidyverse)
library(tidyterra)
library(terra)
library(gridExtra)

source("isdm_fn.R")

generate_spatpred <- function(result, species, suffix, inat_betas, thin) {
  source("prep_dat_forISDM.R")
  
  scaling_factors <- result$scaling_factors
  #get_scaling_factors(species, seed = 0159876, subset = 1)
  
  covar_values_scaled <- as.data.frame(grid)
  for (i in 1:nrow(scaling_factors)) {
    covar_values_scaled[, paste0(scaling_factors$covar[i], "_scaled")] <- 
      (covar_values_scaled[, scaling_factors$covar[i]] - scaling_factors$mean[i]) /
      scaling_factors$sd[i]
  }
  covar_values_scaled$Intercept <- 1
  
  if (inat_betas) {
    target_str <- "inat_beta["
    target_str_sp <- "inat_spat_ranef["
  } else {
    target_str <- "lambda_beta["
    target_str_sp <- "spat_ranef["
  }
  
  # Figure out if we're dealing with an explicitly spatial model
  spatial_model_fit <- any(grepl("spat_ranef", result$summary_list_spatial[[1]]$param))
  
  if (spatial_model_fit) {
    # 
    # thin_ratio <- result$thin2 / result$thin
    # 
    # # Predict occupancy
    # samples <- result$samples_list %>% 
    #   lapply(function(x) {
    #     as.data.frame(as.matrix(x[[1]]))
    #   }) %>% 
    #   bind_rows() %>% 
    #   select(lambda_intercept, 
    #          all_of(paste0(target_str, 1:length(result$occ_covars), "]"))) %>% 
    #   as.matrix()
    # 
    # samples_spatial <- result$samples_list %>% 
    #   lapply(function(x) {
    #     as.data.frame(as.matrix(x[[2]]))
    #   }) %>% 
    #   bind_rows() %>% 
    #   select(all_of(paste0(target_str_sp, 1:max(spec_range_cells[, gridcol]), "]"))) %>% 
    #   as.matrix()
    # 
    # 
    # ### Start the loop
    # occ_covars <- c("Intercept", result$occ_covars)
    # logit_occ_sum <- logit_occ_sumofsquares <- numeric(nrow(covar_values_scaled))
    # covar_mtx <- as.matrix(covar_values_scaled[, occ_covars])
    # spatcells <- unname(unlist(covar_values_scaled[, gridcol]))
    # 
    # pb <- progress::progress_bar$new(total = nrow(samples)/thin_ratio)
    # for (i in 1:(nrow(samples)/thin_ratio)) {
    #   
    #   pb$tick()
    #   logit_occ_sum <- logit_occ_sum + 
    #     covar_mtx %*% samples[i*thin_ratio, ] +
    #     samples_spatial[i, spatcells]
    #   
    #   logit_occ_sumofsquares <- logit_occ_sumofsquares + 
    #     (covar_mtx %*% samples[i*thin_ratio, ] +
    #        samples_spatial[i, spatcells])^2
    # }
    # covar_values_scaled$logit_occ_mean <- logit_occ_sum / (nrow(samples)/thin_ratio)
    # covar_values_scaled$logit_occ_sd <- sqrt(
    #   logit_occ_sumofsquares / (nrow(samples)/thin_ratio) - (covar_values_scaled$logit_occ_mean)^2
    # )
    
  } else {
    
    # Predict occupancy
    samples <- result$samples_list %>% 
      lapply(function(x) {
        as.data.frame(as.matrix(x[[1]]))
      }) %>% 
      bind_rows() %>% 
      select(lambda_intercept, 
             all_of(paste0(target_str, 1:length(result$occ_covars), "]"))) %>% 
      as.matrix()
    
    ### Start the loop
    occ_covars <- c("Intercept", result$occ_covars)
    logit_occ_samples <- matrix(nrow = nrow(samples) / thin,
                                             ncol = nrow(covar_values_scaled))
    covar_mtx <- as.matrix(covar_values_scaled[, occ_covars])
    
    pb <- progress::progress_bar$new(total = nrow(samples)/thin)
    for (i in 1:(nrow(samples)/thin)) {
      
      pb$tick()
      logit_occ_samples[i, ] <- covar_mtx %*% samples[i*thin, ]
    }
    
    covar_values_scaled$logit_occ_mean <- colMeans(logit_occ_samples)
    covar_values_scaled$logit_occ_sd <- apply(logit_occ_samples, 2, sd)
    
  }
  return(covar_values_scaled)
}

plot_spatpred <- function(result, species, suffix, write.out = TRUE, inat_betas = FALSE,
                          thin = 50) {
  
  covar_values_scaled <- generate_spatpred(result, species, suffix, inat_betas, thin)
  logit_occ_mean <- covar_values_scaled$logit_occ_mean
  logit_occ_sd <- covar_values_scaled$logit_occ_sd
  #### Plot of occupancy
  
  occ_mean <- expit(logit_occ_mean)
  
  occ_mean_grid <- grid
  occ_mean_grid$occ_mean <- occ_mean
  occ_mean_grid$logit_occ_sd <- logit_occ_sd
  
  p1 <- ggplot() +
    geom_spatvector(aes(fill = occ_mean), data = occ_mean_grid) +
    scale_fill_viridis_c("Occupancy", na.value = NA)
  
  ### Plot of uncertainty

  p2 <- ggplot() +
    geom_spatvector(aes(fill = logit_occ_sd), data = occ_mean_grid) +
    scale_fill_viridis_c("Occupancy", na.value = NA)
  
  
  bigplot <- gridExtra::arrangeGrob(p1, p2, nrow = 1, 
                                    top = paste0(species, suffix))
  if (write.out) {
    ggsave(plot = bigplot, paste0("plots/mod_compare_", species, suffix, ".jpg"),
           width = 12, height = 5)
  }
  return(list(
    p1 = p1,
    p2 = p2,
    arranged = bigplot
  ))
}

plot_spatpred_difference <- function(result3, result4, species, suffix, 
                                     write.out = TRUE, thin = 50) {
  
  covar_values_scaled_Joint <- generate_spatpred(result3, species, suffix, inat_betas = F, thin = thin)
  covar_values_scaled_CT    <- generate_spatpred(result4, species, suffix, inat_betas = F, thin = thin)
  covar_values_scaled_iNat  <- generate_spatpred(result4, species, suffix, inat_betas = T, thin = thin)
  
  #### Diff. in occu estimate
  
  occ_grid <- grid
  occ_grid$occ_mean_CTmJoint <- expit(covar_values_scaled_CT$logit_occ_mean) - 
    expit(covar_values_scaled_Joint$logit_occ_mean)
  occ_grid$occ_mean_iNatmJoint <- expit(covar_values_scaled_iNat$logit_occ_mean) - 
    expit(covar_values_scaled_Joint$logit_occ_mean)
  occ_grid$sd_CTmJoint <- covar_values_scaled_CT$logit_occ_sd - 
    covar_values_scaled_Joint$logit_occ_sd
  occ_grid$sd_iNatmJoint <- covar_values_scaled_iNat$logit_occ_sd - 
    covar_values_scaled_Joint$logit_occ_sd
  
  p1 <- ggplot() +
    geom_spatvector(aes(fill = occ_mean_iNatmJoint), data = occ_grid) +
    ggtitle("Diff. (iNat minus joint)") +
    scale_fill_gradient2("dOccu", mid = "lightgray", na.value = NA) +
    theme_minimal()
  p2 <- ggplot() +
    geom_spatvector(aes(fill = occ_mean_CTmJoint), data = occ_grid) +
    ggtitle("Diff. (CT minus joint)") +
    scale_fill_gradient2("dOccu", mid = "lightgray", na.value = NA) +
    theme_minimal()
  
  p3 <- ggplot() +
    geom_spatvector(aes(fill = sd_iNatmJoint), data = occ_grid) +
    ggtitle("Diff. in cloglog-scale SD (iNat minus joint)") +
    scale_fill_gradient2("dSD", mid = "lightgray", na.value = NA) +
    theme_minimal()
  p4 <- ggplot() +
    geom_spatvector(aes(fill = sd_CTmJoint), data = occ_grid) +
    ggtitle("Diff. in cloglog-scale SD (CT minus joint)") +
    scale_fill_gradient2("dSD", mid = "lightgray", na.value = NA) +
    theme_minimal()
  
  
  
  
  bigplot <- arrangeGrob(p1, p2, p3, p4, nrow = 2, top = species)
  if (write.out) {
    ggsave(plot = bigplot, paste0("plots/diffplot_", species, suffix, ".jpg"),
           width = 12, height = 10)
  }
  return(list(
    p1 = p1,
    p2 = p2,
    p3 = p3,
    p4 = p4,
    arranged = bigplot
  ))
}


#### Loop over specs, execute ####
source("prep_dat_forISDM.R")
covar_ests_list <- list()

for (i in 1:length(target_species_clean)) {
  res3 <- readRDS(paste0("intermediate/integration_results/test_kfold_", 
                         target_species_clean[i], "_T3_weakdet_DCcomp_samples.RDS"))
  res4 <- readRDS(paste0("intermediate/integration_results/test_kfold_", 
                         target_species_clean[i], "_T4_separate_DCcomp_samples.RDS"))
  
  t3   <- plot_spatpred(res3, target_species_clean[i], write.out = T, suffix = "_T3")
  t4   <- plot_spatpred(res4, target_species_clean[i], write.out = T, suffix = "_T4_CT_only")
  t4.2 <- plot_spatpred(res4, target_species_clean[i], write.out = T, suffix = "_T4_iNat_only", 
                        inat_betas = TRUE)
  
  comp <- gridExtra::arrangeGrob(t3$p1 + ggtitle("Shared"),
                                 t4$p1 + ggtitle("CT only"),
                                 t4.2$p1 + ggtitle("iNat only"), nrow = 1,
                                 top = target_species_clean[i])
  ggsave(plot = comp, 
         paste0("plots/compare_joint_separate_",  target_species_clean[i], ".jpg"), 
         width = 12, height = 4)
  
  diffplot <- plot_spatpred_difference(res3, res4, target_species_clean[i], "_diff", thin = 50)
}

#### Loop over specs. again to make covar plots ####
covar_df_list <- list()
theta_df_list <- list()
brier_df_list <- list()
for (i in 1:length(target_species_clean)) {
  res3 <- readRDS(paste0("intermediate/integration_results/test_kfold_", 
                         target_species_clean[i], "_T3_weakdet_DCcomp_samples.RDS"))
  res4 <- readRDS(paste0("intermediate/integration_results/test_kfold_", 
                         target_species_clean[i], "_T4_separate_DCcomp_samples.RDS"))
  
  covars_joint <- res3$summary_list[[1]] %>% 
    filter(grepl("lambda_beta", param)) %>% 
    mutate(type = "Joint")
  covars_CT <- res4$summary_list[[1]] %>% 
    filter(grepl("lambda_beta", param)) %>% 
    mutate(type = "CT only")
  covars_iNat <- res4$summary_list[[1]] %>% 
    filter(grepl("inat_beta", param)) %>% 
    mutate(type = "iNat only")
  
  covar_df_list[[i]] <- bind_rows(covars_joint, covars_CT, covars_iNat) %>% 
    mutate(species = target_species_clean[i])
  
  theta_df_list[[i]] <- res3$summary_list[[1]] %>% 
    filter(grepl("theta", param)) %>% 
    mutate(type = "Joint") %>% 
    mutate(species = target_species_clean[i])
  
  brier_df_list[[i]] <- bind_rows(
      res3$summary_list[[1]] %>% filter(grepl("brier", param)) %>% mutate(type = "Joint"),
      res4$summary_list[[1]] %>% filter(grepl("brier", param)) %>% mutate(type = "CT only")
    ) %>% 
    mutate(species = target_species_clean[i])
}

modtype_colors <- c(
  "iNat only" = "#e41a1c",
  "CT only" = "#377eb8",
  "Joint" = "#984ea3"
)

covar_plot <- covar_df_list %>% 
  bind_rows() %>% 
  ggplot() + 
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(parname, mean, ymin = `2.5%`, ymax = `97.5%`, col = type),
                  position = position_jitter(width = 0.2)) +
  facet_wrap(~species, nrow = 5) +
  coord_flip() + 
  theme_minimal() + 
  scale_color_manual("Model type", values = modtype_colors)

ggsave("plots/covar_summary.jpg", plot = covar_plot, width = 6.5, height = 7)

theta_plot <- theta_df_list %>% 
  bind_rows() %>% 
  ggplot() + 
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(param, mean, ymin = `2.5%`, ymax = `97.5%`, col = type),
                  position = position_jitter(width = 0.2)) +
  facet_wrap(~species, nrow = 5) +
  coord_flip() + 
  theme_minimal() + 
  scale_color_manual("Model type", values = modtype_colors)

ggsave("plots/theta_summary.jpg", plot = theta_plot, width = 6.5, height = 4)

brier_plot <- brier_df_list %>% 
  bind_rows() %>% 
  select(mean, species, type) %>% 
  pivot_wider(names_from = type, values_from = mean) %>%
  ggplot() + 
  geom_hline(yintercept = 0) +
  geom_point(aes(species, `CT only` - Joint),
                  position = position_jitter(width = 0.2)) +
  coord_flip() + 
  theme_minimal() +
  xlab("Improvement in Brier score (OOS predictive performance) due to joint model")

