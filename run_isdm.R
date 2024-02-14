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
    modtype = "T3_weakdet",
    species = target_species_clean[i],
    spatial_model = F,
    ni = 5000, nb = 1000,
    nc = 3, nt = 2, nt2 = 20,
    seed = 36405, subset = 1,
    subset_inat = 1, suffix = "_DCcomp"
  )
  '), file = this_outfile)
  
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
    ni = 15000, nb = 5000,
    nc = 3, nt = 2, nt2 = 20,
    seed = 36405, subset = 1, 
    subset_inat = 1, suffix = "_DCcomp"
  )
  '), file = this_outfile)
}


cl <- makeCluster(4)

parLapply(cl, outfiles, source)

