library(tidyverse)
library(parallel)

source("isdm_fn.R")
source("model_code.R")
source("prep_dat_forISDM.R")


preamble <- '
library(nimbleEcology)
library(tidyverse)

source("isdm_fn.R")
source("model_code.R")
source("prep_dat_forISDM.R")


i <- '

outfiles <- c()
ct <- 0
# Run spatial models
for (i in 1:length(target_species)) {
  ct <- ct + 1
  this_outfile <- paste0("temp/source/script", ct, ".R")
  outfiles <- c(outfiles, this_outfile)
  
  cat(paste0(preamble, i, '
  
  fit_integrated_model_with_CV(
    nfolds = 5               ,
    kfold_mode = "single"    ,
    folds_by = "camera"      ,
    modtype = "T4_separate",
    species = target_species_clean[i],
    spatial_model = F,
    ni = 5000, nb = 1000,
    nc = 3, nt = 2, nt2 = 20,
    seed = 36405, subset = 1,
    inat_only_strong = TRUE,
    subset_inat = 1, suffix = "_iNatOnly"
  )
  '), file = this_outfile)
}


lapply(outfiles, source)

param_df_list <- list()
for (i in 1:length(target_species_clean)) {
  res_inatOnly <- readRDS(paste0("intermediate/integration_results/test_kfold_",
                          target_species_clean[i], 
                          "_T4_separate_iNatOnly_samples.RDS"))
  res4 <-  readRDS(paste0("intermediate/integration_results/test_kfold_",
                   target_species_clean[i], 
                   "_T4_separate_DCcomp_samples.RDS"))
  
  param_df_list[[i]] <- bind_rows(
    res4$summary_list[[1]] %>% 
      filter(grepl("inat_beta", param)) %>% 
      mutate(species = target_species_clean[i], type = "T4"),
    res4$summary_list[[1]] %>% 
      filter(grepl("inat_beta", param)) %>% 
      mutate(species = target_species_clean[i], type = "Strong")
  )
}
bind_rows(param_df_list) %>% 
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(param, mean, ymin = `2.5%`, ymax = `97.5%`, col = type),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~species, nrow = 5) +
  coord_flip() + 
  theme_minimal()
