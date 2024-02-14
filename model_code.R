library(nimbleEcology)


#### calcBrierScore_oneSite ####
calcBrierScore_oneSite <- nimbleFunction(
  run = function(y = double(1),
                 psi = double(0),
                 p = double(1),
                 len = double(0)){
    
    x <- y[1:len] * (1 - psi * p[1:len])^2 +
      (1 - y[1:len]) * (psi * p[1:len])^2
    
    return(sum(x))
    returnType(double(0))
  })

#### calcIntensity ####
calcIntensity_strong <- nimbleFunction(run = function(
    lambda_intercept = double(0),
    lambda_beta = double(1),
    xdat = double(2),
    spatial_ranef = double(0)) {
  
  log_lambda <- log(lambda_intercept) + xdat %*% lambda_beta + spatial_ranef
  
  mu <- sum(exp(log_lambda))
  
  return(mu)
  returnType(double(0))
})


calcIntensity_weak <- nimbleFunction(run = function(
    lambda_intercept = double(0),
    theta0 = double(0),
    theta1 = double(0),
    lambda_beta = double(1),
    xdat = double(2),
    spatial_ranef = double(0)) {
  
  log_lambda <- log(lambda_intercept) + xdat %*% lambda_beta
  
  mu <- exp(theta0 + theta1 * sum(exp(log_lambda)))
  
  return(mu)
  returnType(double(0))
})



#### February proposal model ####
# 3 possible types: T1_strong, T2_correlated, T3_weakdet, T4_separate

feb_proposal_code <- nimbleCode({
  
  ### iNaturalist data distribution
  for (g in 1:nINatCell) {
    inat_ct[g] ~ dpois(inat_effort[g] * mu[g])
    
    ### iNaturalist grid cell intensities
    if (modtype == "T1_strong") {
      log(mu[g]) <- lambda_intercept + 
        inprod(lambda_beta[1:nLamBeta], xdat_g[g, 1:nLamBeta])
    } else if (modtype %in% c("T2_correlated", "T3_weakdet")) {
      log(mu[g]) <- theta1 * (lambda_intercept + 
        inprod(lambda_beta[1:nLamBeta], xdat_g[g, 1:nLamBeta])) + theta0
    } else {
      log(mu[g]) <- inat_intercept + inprod(inat_beta[1:nLamBeta], xdat_g[g, 1:nLamBeta])
    }
  }
  
  ### Detection probabilities
  if (modtype %in% c("T1_strong", "T2_correlated")) {
    # Shared detection beta models
    for (o in 1:nobs) {
      cloglog(p[o]) <- 
        cloglog(p_intercept) + 
        theta2 * inprod(lambda_beta[1:nLamBeta], xdat[deployment_ID[o], 1:nLamBeta]) +
        inprod(p_gamma[1:nGamma], wdat[o, 1:nGamma])
    }
    for (o in 1:nobs_holdout) {
      cloglog(p_holdout[o]) <- 
        cloglog(p_intercept) + 
        theta2 * inprod(lambda_beta[1:nLamBeta], xdat[deployment_ID_holdout[o], 1:nLamBeta]) +
        inprod(p_gamma[1:nGamma], wdat_holdout[o, 1:nGamma])
    }
  } else {
    # Separate detection beta models
    for (o in 1:nobs) {
      cloglog(p[o]) <- 
        cloglog(p_intercept) + 
        inprod(p_alpha[1:nLamBeta], xdat[deployment_ID[o], 1:nLamBeta]) +
        inprod(p_gamma[1:nGamma], wdat[o, 1:nGamma])
    }
    for (o in 1:nobs_holdout) {
      cloglog(p_holdout[o]) <- 
        cloglog(p_intercept) + 
        inprod(p_alpha[1:nLamBeta], xdat[deployment_ID_holdout[o], 1:nLamBeta]) +
        inprod(p_gamma[1:nGamma], wdat_holdout[o, 1:nGamma])
    }
  }
  
  
  # Occupancy probabilities
  for (i in 1:ncameras_indat) {
    cloglog(psi[i]) <- log(lambda_intercept) +
      inprod(lambda_beta[1:nLamBeta], 
             xdat[deployment_ID[start[i]], 1:nLamBeta])
  }
  for (i in 1:ncameras_holdout) {
    cloglog(psi_holdout[i]) <- log(lambda_intercept) +
      inprod(lambda_beta[1:nLamBeta], 
             xdat[deployment_ID_holdout[start_holdout[i]], 1:nLamBeta])
  }
  
  
  # Occupancy model distributions
  for (i in 1:ncameras_indat) {
    y[start[i]:end[i]] ~ dOcc_v(
      probOcc = psi[i], 
      probDetect = p[start[i]:end[i]],
      len = end[i] - start[i] + 1
    )
  }
  for (i in 1:ncameras_holdout) {
    y_holdout[start_holdout[i]:end_holdout[i]] ~ dOcc_v(
      probOcc = psi_holdout[i], 
      probDetect = p_holdout[start_holdout[i]:end_holdout[i]],
      len = end_holdout[i] - start_holdout[i] + 1
    )
    
    pp_brier_bysite[i] <- calcBrierScore_oneSite(
      y = y_holdout[start_holdout[i]:end_holdout[i]],
      psi = psi_holdout[i],
      p = p_holdout[start_holdout[i]:end_holdout[i]],
      len = end_holdout[i] - start_holdout[i] + 1
    )
  }
  
  ### Brier score  
  brier_score <- sum(pp_brier_bysite[1:ncameras_holdout])
  
  ### Priors we always need
  p_intercept      ~ dunif(0, 1)
  lambda_intercept ~ dgamma(shape = 1, scale = 50)
  for (i in 1:nLamBeta) {
    lambda_beta[i] ~ dnorm(0, sd = 2.75)
  }
  for (i in 1:nGamma) {
    p_gamma[i] ~ dnorm(0, sd = 2.75)
  }
  
  ### Modtype-specific priors
  
  # iNat correlation params
  if (modtype %in% c("T2_correlated", "T3_weakdet")) {
    theta0 ~ dnorm(0, sd = 2.75)
    theta1 ~ dnorm(0, sd = 2.75)
  }
  
  # Detection/intensity correlation param
  if (modtype %in% c("T1_strong")) {
    theta2 <- 1
  } else if (modtype == "T2_correlated") {
    theta2 ~ dnorm(1, sd = 2.75)
  }
  
  # Separate params for detection
  if (modtype %in% c("T3_weakdet", "T4_separate")) {
    for (i in 1:nLamBeta) {
      p_alpha[i] ~ dnorm(0, sd = 2.75)
    }
  }
  
  # Separate iNat params for everything
  if (modtype == "T4_separate") {
    inat_intercept ~ dgamma(shape = 1, scale = 25)
    
    for (i in 1:nLamBeta) {
      inat_beta[i] ~ dnorm(0, sd = 2.75)
    }
  }
  
  # Dummy variable so we always have something in thin2
  dummy ~ dunif(0, 1)
})


#### iNat-only model ####
inat_only_code <- nimbleCode({
  
  ### iNaturalist data distribution
  for (g in 1:nINatCell) {
    inat_ct[g] ~ dpois(inat_effort[g] * mu[g])
    
    log(mu[g]) <- inat_intercept + inprod(inat_beta[1:nLamBeta], xdat_g[g, 1:nLamBeta])
  }
  
  # Separate iNat params for everything
  inat_intercept ~ dgamma(shape = 1, scale = 25)
  
  for (i in 1:nLamBeta) {
    inat_beta[i] ~ dnorm(0, sd = 2.75)
  }
  
  dummy ~ dunif(0,1)
})


#### Feb. model with spatial effects ####

# feb_proposal_code_wSpat <- nimbleCode({
#   
#   ### iNaturalist data distribution
#   for (g in 1:nINatCell) {
#     inat_ct[g] ~ dpois(inat_effort[g] * mu[g])
#     
#     
#     ### iNaturalist grid cell intensities
#     if (modtype == "T1_strong") {
#       mu[g] <- calcIntensity_strong(
#         lambda_intercept = lambda_intercept,
#         lambda_beta = lambda_beta[1:nLamBeta],
#         xdat = xdat_allS2[S2toS3_start[g]:S2toS3_end[g], 1:nLamBeta],
#         spatial_ranef = spat_ranef[spatcell_inat[g]] * spat_magnitude
#       )
#     } else if (modtype %in% c("T2_correlated", "T3_weakdet")) {
#       mu[g] <- calcIntensity_weak(
#         lambda_intercept = lambda_intercept,
#         theta0 = theta0, theta1 = theta1,
#         lambda_beta = lambda_beta[1:nLamBeta],
#         xdat = xdat_allS2[S2toS3_start[g]:S2toS3_end[g], 1:nLamBeta],
#         spatial_ranef = spat_ranef[spatcell_inat[g]] * spat_magnitude
#       )
#     } else {
#       mu[g] <- calcIntensity_strong(
#         lambda_intercept = inat_intercept,
#         lambda_beta = inat_beta[1:nLamBeta],
#         xdat = xdat_allS2[S2toS3_start[g]:S2toS3_end[g], 1:nLamBeta],
#         spatial_ranef = inat_spat_ranef[spatcell_inat[g]] * inat_spat_magnitude
#       )
#     }
#   }
#   
#   ### Detection probabilities
#   if (modtype %in% c("T1_strong", "T2_correlated")) {
#     # Shared detection beta models
#     for (o in 1:nobs) {
#       cloglog(p[o]) <- 
#         cloglog(p_intercept) + 
#         theta2 * inprod(lambda_beta[1:nLamBeta], xdat[deployment_ID[o], 1:nLamBeta]) +
#         inprod(p_gamma[1:nGamma], wdat[o, 1:nGamma]) +
#         spat_ranef[spatcell[deployment_ID[o]]] * spat_magnitude
#     }
#     for (o in 1:nobs_holdout) {
#       cloglog(p_holdout[o]) <- 
#         cloglog(p_intercept) + 
#         theta2 * inprod(lambda_beta[1:nLamBeta], xdat[deployment_ID_holdout[o], 1:nLamBeta]) +
#         inprod(p_gamma[1:nGamma], wdat_holdout[o, 1:nGamma]) +
#         spat_ranef[spatcell[deployment_ID_holdout[o]]] * spat_magnitude
#     }
#   } else {
#     # Separate detection beta models
#     for (o in 1:nobs) {
#       cloglog(p[o]) <- 
#         cloglog(p_intercept) + 
#         inprod(p_alpha[1:nLamBeta], xdat[deployment_ID[o], 1:nLamBeta]) +
#         inprod(p_gamma[1:nGamma], wdat[o, 1:nGamma])
#     }
#     for (o in 1:nobs_holdout) {
#       cloglog(p_holdout[o]) <- 
#         cloglog(p_intercept) + 
#         inprod(p_alpha[1:nLamBeta], xdat[deployment_ID_holdout[o], 1:nLamBeta]) +
#         inprod(p_gamma[1:nGamma], wdat_holdout[o, 1:nGamma])
#     }
#   }
#   
#   
#   # Occupancy probabilities
#   for (i in 1:ncameras_indat) {
#     cloglog(psi[i]) <- log(lambda_intercept) +
#       inprod(lambda_beta[1:nLamBeta], 
#              xdat[deployment_ID[start[i]], 1:nLamBeta]) +
#       spat_ranef[spatcell[deployment_ID[start[i]]]] * spat_magnitude
#   }
#   for (i in 1:ncameras_holdout) {
#     cloglog(psi_holdout[i]) <- log(lambda_intercept) +
#       inprod(lambda_beta[1:nLamBeta], 
#              xdat[deployment_ID_holdout[start_holdout[i]], 1:nLamBeta]) +
#       spat_ranef[spatcell[deployment_ID_holdout[start[i]]]] * spat_magnitude
#   }
#   
#   
#   # Occupancy model distributions
#   for (i in 1:ncameras_indat) {
#     y[start[i]:end[i]] ~ dOcc_v(
#       probOcc = psi[i], 
#       probDetect = p[start[i]:end[i]],
#       len = end[i] - start[i] + 1
#     )
#   }
#   for (i in 1:ncameras_holdout) {
#     y_holdout[start_holdout[i]:end_holdout[i]] ~ dOcc_v(
#       probOcc = psi_holdout[i], 
#       probDetect = p_holdout[start_holdout[i]:end_holdout[i]],
#       len = end_holdout[i] - start_holdout[i] + 1
#     )
#     
#     pp_brier_bysite[i] <- calcBrierScore_oneSite(
#       y = y_holdout[start_holdout[i]:end_holdout[i]],
#       psi = psi_holdout[i],
#       p = p_holdout[start_holdout[i]:end_holdout[i]],
#       len = end_holdout[i] - start_holdout[i] + 1
#     )
#   }
#   
#   ### Brier score  
#   brier_score <- sum(pp_brier_bysite[1:ncameras_holdout])
#   
#   ### Priors we always need
#   p_intercept      ~ dunif(0, 1)
#   lambda_intercept ~ dgamma(shape = 1, scale = 50)
#   
#   for (i in 1:nLamBeta) {
#     lambda_beta[i] ~ dnorm(0, sd = 2.75)
#   }
#   for (i in 1:nGamma) {
#     p_gamma[i] ~ dnorm(0, sd = 2.75)
#   }
#   
#   # Spatial priors
#   spat_magnitude   ~ T(dnorm(0, sd = 2.5), 1e-6, 20)
#   spat_ranef[1:nSpatCell] ~ dcar_normal(adj = adj[1:adj_len],
#                                         num = num[1:nSpatCell],
#                                         tau = 1, zero_mean = 1)
#   
#   
#   ### Modtype-specific priors
#   
#   # iNat correlation params
#   if (modtype %in% c("T2_correlated", "T3_weakdet")) {
#     theta0 ~ dnorm(0, sd = 2.75)
#     theta1 ~ dnorm(0, sd = 2.75)
#   }
#   
#   # Detection/intensity correlation param
#   if (modtype %in% c("T1_strong")) {
#     theta2 <- 1
#   } else if (modtype == "T2_correlated") {
#     theta2 ~ dnorm(1, sd = 2.75)
#   }
#   
#   # Separate params for detection
#   if (modtype %in% c("T3_weakdet", "T4_separate")) {
#     for (i in 1:nLamBeta) {
#       p_alpha[i] ~ dnorm(0, sd = 2.75)
#     }
#   }
#   
#   # Separate iNat params for everything
#   if (modtype == "T4_separate") {
#     inat_intercept ~ dgamma(shape = 1, scale = 25)
#     
#     for (i in 1:nLamBeta) {
#       inat_beta[i] ~ dnorm(0, sd = 2.75)
#     }
#     
#     inat_spat_magnitude ~ T(dnorm(0, sd = 2.5), 1e-6, 20)
#     inat_spat_ranef[1:nSpatCell] ~ dcar_normal(adj = adj[1:adj_len],
#                                                num = num[1:nSpatCell],
#                                                tau = 1, zero_mean = 1)
#   }
#   # Dummy variable so we always have something in thin2
#   dummy ~ dunif(0, 1)
# })

