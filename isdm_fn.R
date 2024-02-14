


get_monitors <- function(modtype, spatial_model) {
  
  monitors <- c(
    "lambda_intercept",
    "p_intercept",
    "lambda_beta",
    "p_gamma",
    "brier_score"
  )
  
  if (spatial_model) {
    monitors <- c(monitors, "spat_magnitude")
  }
  
  if (modtype == "T2_correlated") {
    monitors <- c(monitors, "theta0", "theta1", "theta2")
  } else if (modtype == "T3_weakdet") {
    monitors <- c(monitors, "theta0", "theta1", "p_alpha")
  } else if (modtype == "T4_separate") {
    monitors <- c(monitors, "p_alpha", "inat_beta", "inat_intercept")
    if (spatial_model) {
      monitors <- c(monitors, "inat_spat_magnitude")
    }
  }
  
  return(monitors)
}


dc_inits <- function(modtype, occ_covars, det_covars, spatial_model,
                      this_adj_info = NULL) {
  inits_list <- list(
    # Always
    lambda_intercept = 0.2,
    p_intercept = 0.5,
    lambda_beta = rnorm(length(occ_covars), 0, sd = 0.001),
    p_gamma = rnorm(length(det_covars), 0, sd = 0.001),
    dummy = 1
  )
  
  if (spatial_model) {
    inits_list[["spat_magnitude"]] <- 1
    inits_list[["spat_ranef"]] <- rnorm(length(this_adj_info$num), 0, 0.1)
  }
  
  if (modtype == "T2_correlated") {
    inits_list[["theta0"]] <- 0
    inits_list[["theta1"]] <- 0.1
    inits_list[["theta2"]] <- 1
  } else if (modtype == "T3_weakdet") {
    inits_list[["theta0"]] <- 0
    inits_list[["theta1"]] <- 0.1
    inits_list[["p_alpha"]] <- rnorm(length(occ_covars), 0, sd = 0.001)
  } else if (modtype == "T4_separate") {
    inits_list[["p_alpha"]] <- rnorm(length(occ_covars), 0, sd = 0.001)
    inits_list[["inat_beta"]] <- rnorm(length(occ_covars), 0, sd = 0.001)
    inits_list[["inat_intercept"]] <- 1
    if (spatial_model) {
      inits_list[["inat_spat_magnitude"]] <- 0.2
      inits_list[["inat_spat_ranef"]] <- rnorm(length(this_adj_info$num), 0, 0.001)
    }
  }
  
  return(inits_list)
}



fit_integrated_model_with_CV <- function(nfolds, 
                                         kfold_mode = "kfold",
                                         folds_by = "camera", 
                                         modtype,
                                         modstrength,
                                         species,
                                         ni, 
                                         nb, 
                                         nc, 
                                         nt, 
                                         nt2,
                                         seed,
                                         subset = 1,
                                         subset_inat = 1,
                                         spatial_model = TRUE,
                                         inat_only_strong = FALSE,
                                         suffix = "") {
  set.seed(seed)
  start_time <- Sys.time()
  
  # Check that args are valid
  stopifnot(folds_by %in% c("camera", "subproject"))
  stopifnot(kfold_mode %in% c("kfold", "single"))
  stopifnot(modtype %in% c("T1_strong", "T2_correlated", "T3_weakdet", "T4_separate"))
  stopifnot(species %in% target_species_clean)
  spec_index <- which(target_species_clean == species)


  if (spatial_model) {
    stop("Spatial model is not implemented")
    # model_code <- feb_proposal_code_wSpat
  } else {
    if (inat_only_strong) {
      model_code <- inat_only_code
    } else {
      model_code <- feb_proposal_code
    }
  }
  
  ### Load in and format the data
  source("prep_dat_forISDM.R")
  
  # Scale occupancy covariates
  covar_df <- as.data.frame(grid)
  covars_toscale <- c("m_water", "m_traffic", "m_canopy", "m_popdens", "m_income")
  scaling_factors <- data.frame(
    covar = covars_toscale, mean = NA, sd = NA
  )
  for (i in 1:length(covars_toscale)) {
    scaling_factors$mean[i] <- mean(covar_df[, covars_toscale[i]])
    scaling_factors$sd[i]   <- sd(covar_df[, covars_toscale[i]])
    covar_df[[paste0(covars_toscale[i], "_scaled")]] <- as.numeric(
      (covar_df[, covars_toscale[i]] - scaling_factors$mean[i]) / 
        scaling_factors$sd[i]
    )
  }
  
  dethist <- dethist_list[[spec_index]]
  
  dethist$yday_scaled <- as.numeric(scale(dethist$yday_start))
  dethist$yday_scaled_sq <- dethist$yday_scaled^2
  

  occ_covars <- c("m_water_scaled", 
                  "m_traffic_scaled", 
                  "m_canopy_scaled", 
                  "m_popdens_scaled", 
                  "m_income_scaled")
  det_covars <- c("yday_scaled", "yday_scaled_sq")
  
  siteID_transl <- data.frame(
    deployment_id = unique(dethist$deployment_id)
  ) %>% 
    mutate(site_ID = row_number())
  
  
  y.all        <- dethist %>% 
    left_join(siteID_transl, by = "deployment_id") %>% 
    arrange(site_ID) %>% 
    rename(value = observed)
  
  inds_to_keep <- which(!is.na(y.all$value))
  y.all <- y.all[inds_to_keep,]
  
  occ.covs.all <- left_join(
      distinct(dethist, deployment_id, block_ID),
      covar_df, by = "block_ID"
    ) %>% 
    left_join(siteID_transl, by = "deployment_id")
  
  occ.covs.all$Intercept <- 1
  
  det.covs.all <- y.all[, c("site_ID", det_covars)]
  det.covs.all$Intercept <- 1
  
  stopifnot(nrow(det.covs.all) == nrow(y.all))
  stopifnot(nrow(occ.covs.all) == max(y.all$site_ID))
  
  ### Drop sites with only one obs. (for NIMBLE, for now)
  sites_to_drop <- count(y.all, site_ID) %>% filter(n == 1) 
  
  if (any(rowSums(is.na(occ.covs.all[, gsub("_scaled", "", occ_covars)])) > 0)) {
    sites_to_drop <- sites_to_drop %>% 
      bind_rows(data.frame(
        site_ID = occ.covs.all$site_ID[rowSums(is.na(occ.covs.all[, gsub("_scaled", "", occ_covars)])) > 0],
        n = 1
      ))
  }
  
  warning(paste0("Dropping ", sum(rowSums(is.na(occ.covs.all[, gsub("_scaled", "", occ_covars)])) > 0), 
                 " cameras due to missing covars."))
  
  if (subset < 1) {
    testSubset <- sample_frac(occ.covs.all, size = 1-subset)
    
    sites_to_drop <- 
      bind_rows(sites_to_drop, data.frame(site_ID = testSubset$site_ID, n = 1)) %>% 
      distinct()
  }
  
  if (nrow(sites_to_drop) > 0) {
    
    new_site_ID_df <- data.frame(
      old_site_ID = occ.covs.all$site_ID,
      included = !(occ.covs.all$site_ID %in% sites_to_drop$site_ID),
      new_site_ID = NA
    )
    new_site_ID_df$new_site_ID[new_site_ID_df$included] <-
      1:sum(new_site_ID_df$included)
    
    occ.covs.all <- occ.covs.all %>% 
      left_join(new_site_ID_df, by = c("site_ID" = "old_site_ID")) %>% 
      filter(!is.na(new_site_ID)) %>% 
      mutate(site_ID = new_site_ID)
    
    y.all <- y.all %>% 
      left_join(new_site_ID_df, by = c("site_ID" = "old_site_ID")) %>% 
      filter(!is.na(new_site_ID)) %>% 
      mutate(site_ID = new_site_ID)
    
    det.covs.all <- det.covs.all %>% 
      left_join(new_site_ID_df, by = c("site_ID" = "old_site_ID")) %>% 
      filter(!is.na(new_site_ID)) %>% 
      mutate(site_ID = new_site_ID)
    
  }
  stopifnot(nrow(det.covs.all) == nrow(y.all))
  stopifnot(nrow(occ.covs.all) == max(y.all$site_ID))
  
  ### Set up the holdout data
  if (folds_by == "camera") {
    folds_df <- data.frame(
      camera_ID = sort(unique(y.all$site_ID)),
      fold = sample(rep(1:nfolds, length(unique(y.all$site_ID)))[
        1:length(unique(y.all$site_ID))
      ])
    )
  } else if (folds_by == "subproject") {
    stop("Need to re-run integration_data_prep to produce subproject_name")
    folds_df <- data.frame(
      camera_ID = 1:nrow(y.all),
      subproject_ID = occ.covs.all$subproject_name,
      fold = NA
    )
    
    fold_ID_vec <- sample(rep(1:nfolds, length(unique))[1:nrow(y)])
    
    
  } else {
    stop('Invalid "folds_by" value')
  }
  
  
  inat_grid_cts <- left_join(inat_grid_cts, covar_df) %>% 
    filter(n > 0)
  
  if (subset_inat < 1) {
    stop("subset_inat < 1 not implemented")
  }

  #### Loop over folds.
  # We need to rebuild the model each time because the likelihood structure changes
  # (unavoidable due to diff. number of obs. per camera deployment)
  samples_list <- list()
  summary_list <- list()
  summary_list_spatial <- list()
  
  nfolds_for_loop <- ifelse(kfold_mode == "single", 1, nfolds)
  ### Start loop over folds
  for (fold_i in 1:nfolds_for_loop) {
    holdout_IDs <- folds_df$camera_ID[folds_df$fold == fold_i]
    
    inmod_inds_y   <- which(!(y.all$site_ID %in% holdout_IDs))
    holdout_inds_y <- which(y.all$site_ID %in% holdout_IDs)
    
    y_inmod <- y.all[inmod_inds_y, ]
    y_holdout <- y.all[holdout_inds_y, ]
    
    det_x_inmod <-   det.covs.all[inmod_inds_y, ]
    det_x_holdout <- det.covs.all[holdout_inds_y, ]
    
    # occ_x_inmod <- occ.covs.all[!(occ.covs.all$site_ID %in% holdout_IDs), ]
    # occ_x_holdout <- occ.covs.all[occ.covs.all$site_ID %in% holdout_IDs, ]
    
    # Get the start/end vectors
    siteID_vec <- unique(y_inmod$site_ID)
    start_vec <- end_vec <- numeric(length(siteID_vec))
    
    for (i in 1:length(siteID_vec)) {
      this_id <- siteID_vec[i]
      start_vec[i] <- min(which(y_inmod$site_ID == this_id))
      end_vec[i] <-   max(which(y_inmod$site_ID == this_id))
    }
    
    siteID_vec_holdout <- unique(y_holdout$site_ID)
    start_vec_holdout <- end_vec_holdout <- numeric(length(siteID_vec_holdout))
    
    for (i in 1:length(siteID_vec_holdout)) {
      this_id <- siteID_vec_holdout[i]
      start_vec_holdout[i] <- min(which(y_holdout$site_ID == this_id))
      end_vec_holdout[i] <-   max(which(y_holdout$site_ID == this_id))
    }
    
    # We need:
    # adj/num, the adjacency info
    # spatcell, a vector mapping each row of xdat to elements of num
    # spatcell_inat, a vector mapping each element of inat_ct to the correct scale4 cell
    
    # Build the model
    mod <- nimbleModel(code = model_code,
                       constants = list(
                         nobs = nrow(y_inmod),
                         start = start_vec,
                         end = end_vec,
                         nobs_holdout = nrow(y_holdout),
                         start_holdout = start_vec_holdout,
                         end_holdout = end_vec_holdout,
                         modtype = modtype,
                         
                         nLamBeta = length(occ_covars),
                         nGamma = length(det_covars),
                         
                         ncameras_indat = length(unique(y_inmod$site_ID)),
                         ncameras_holdout = length(unique(y_holdout$site_ID)),
                         deployment_ID = y_inmod$site_ID,
                         deployment_ID_holdout = y_holdout$site_ID, # Should max at ncameras_all
                         nINatCell = nrow(inat_grid_cts)
                         
                       ),
                       data = list(
                         inat_ct     = as.numeric(unlist(inat_grid_cts[, species])),
                         inat_effort = as.numeric(unlist(inat_grid_cts[, "n"])),
                         xdat_g  =     as.matrix(inat_grid_cts[, occ_covars]),
                         
                         xdat =         as.matrix(occ.covs.all[, occ_covars]),
                         y =            as.numeric(y_inmod$value),
                         y_holdout =    as.numeric(y_holdout$value),
                         wdat =         as.matrix(det_x_inmod[, det_covars]),
                         wdat_holdout = as.matrix(det_x_holdout[, det_covars])
                       ),
                       inits = dc_inits(modtype, occ_covars, det_covars,
                                         spatial_model, this_adj_info), 
                       calculate = F)
    cmod <- compileNimble(mod)
    # browser()
    
    mcmcConf <- configureMCMC(mod)

    if (inat_only_strong) {
      mcmcConf$setMonitors(c("inat_beta", "inat_intercept"))
      mcmcConf$setMonitors2("dummy")
      
    } else {
      mcmcConf$setMonitors(get_monitors(modtype, spatial_model))
      monitors2 <- c("dummy")
      if (spatial_model) {
        monitors2 <- c(monitors2, "spat_ranef")
        if (modtype == "T4_separate") monitors2 <- c(monitors2, "inat_spat_ranef")
      }
      
      if (modtype == "T3_weakdet") {
        mcmcConf$removeSampler(c("lambda_intercept", "theta0", "theta1"))
        mcmcConf$addSampler(target = c("lambda_intercept", "theta0", "theta1"),
                            type = "AF_slice")
      }
      
      mcmcConf$setMonitors2(monitors2) # add back in spatial effect when ready
    }
    
    mcmc <- buildMCMC(mcmcConf)
    cmcmc <- compileNimble(mcmc)
    
    mcmc_start_time <- Sys.time()
    # browser()
    samples_list[[fold_i]] <- 
      runMCMC(cmcmc, niter = ni, nburnin = nb, thin = nt, 
              nchains = nc, thin2 = nt2,
              samplesAsCodaMCMC = TRUE, 
              inits = dc_inits(modtype, occ_covars, det_covars,
                                spatial_model, this_adj_info))
    
    summary_list[[fold_i]] <- MCMCvis::MCMCsummary(samples_list[[fold_i]]$samples)
    summary_list[[fold_i]]$param <- rownames(summary_list[[fold_i]])
    summary_list[[fold_i]]$fold <- fold_i
    
    summary_list_spatial[[fold_i]] <- MCMCvis::MCMCsummary(samples_list[[fold_i]]$samples2)
    summary_list_spatial[[fold_i]]$param <- rownames(summary_list_spatial[[fold_i]])
    summary_list_spatial[[fold_i]]$fold <- fold_i
  }
  
  param_inds <- colnames(samples_list[[1]]$samples[[1]])
  brier_index <- which(param_inds == "brier_score")
  logsc_index <- which(param_inds == "logarithmic_score")
  
  brier_score <- as.numeric(unlist(
    lapply(samples_list, function(x) lapply(x$samples, function(y) as.numeric(y[, brier_index])))
  )) %>% mean()
  # logarithmic_score <- as.numeric(unlist(
  #   lapply(samples_list, function(x) lapply(x$samples, function(y) as.numeric(y[, logsc_index])))
  # )) %>% mean()
  logarithmic_score <- NA
  end_time <- Sys.time()
  
  lambda_covar_indices <- c(grep("lambda_beta", param_inds), 
                            grep("p_alpha", param_inds),
                            grep("inat_beta", param_inds))
  det_covar_indices <- grep("p_gamma", param_inds)
  
  for (i in 1:length(summary_list)) {
    summary_list[[i]]$parname <- NA
    summary_list[[i]]$parname[lambda_covar_indices] <- gsub("_scaled", "", occ_covars)
    summary_list[[i]]$parname[det_covar_indices] <- gsub("_scaled", "", det_covars)
  }
  
  saveRDS(list(
    time_taken = end_time - start_time,
    mcmc_time_taken = end_time - mcmc_start_time,
    thin = nt, thin2 = nt2,
    samples_list = samples_list,
    summary_list = summary_list,
    summary_list_spatial = summary_list_spatial,
    nfolds = nfolds, folds_by = folds_by,
    scaling_factors = scaling_factors,
    modtype = modtype, 
    species = species,
    subset = subset, subset_inat = subset_inat,
    brier_score = brier_score, 
    call = match.call(),
    occ_covars = occ_covars, det_covars = det_covars
  ), paste0("intermediate/integration_results/test_kfold_", 
            species, "_", modtype, suffix, "_samples.RDS"))
  
}







