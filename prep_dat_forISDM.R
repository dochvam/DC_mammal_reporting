library(tidyverse)
library(terra)
library(tidyterra)
library(lubridate)

redo_inat <- F
redo_dethist <- F

# Prepare data for iSDM analysis. We need to end up with:
#   Covariates for each grid cell (coming in as polygons)
#   detection histories for each camera
#   iNat counts of each species, and totals, in each grid cell

#### Prepare assets ####

specs_to_drop <- c(
  "Camera Misfire", "Homo sapiens", "Other Bird species",
  "Vehicle", "Camera Trapper", "Unknown Ground-Dove", "Turdus migratorius",
  "Unknown Animal", "No Animal", "Bicycle", "Unknown Bird", "Cyanocitta cristata",
  "Meleagris gallopavo", "Corvus brachyrhynchos", "Colaptes auratus",
  "Coragyps atratus", "Cathartes aura", "Reptile species", "Owl species",
  "Animal Not on List", "Raptor Species", "Melanerpes erythrocephalus",
  "Canis familiaris", "Equus caballus", "Unknown Small Rodent"
)

ct_dat <- readxl::read_xlsx("data/final_DCCC_data.xlsx", sheet = 2) %>% 
  filter(Actual.Long > -77.3) %>% 
  filter(!Species.Name %in% specs_to_drop)

ct_speccts <- ct_dat %>% 
  count(Species.Name, Common.Name) %>% 
  arrange(-n) 

target_species <- ct_speccts$Species.Name[1:10]
target_species_clean <- gsub(" ", "_", target_species)


# Read in the grid
grid <- vect("data/l.grid.shp") %>% 
  mutate(block_ID = as.numeric(as.factor(BlockName)))

#### Get counts of iNaturalist data in each cell ####
if (!file.exists("intermediate/inat_grid_cts.csv") | redo_inat) {
  inat_dat <- readxl::read_xlsx("data/DC_mammal_iNat_data_1262023.xlsx")
  inat_dat <- inat_dat %>% 
    filter(!coordinates_obscured | !is.na(private_latitude),
           is.na(positional_accuracy) | positional_accuracy <= 1000) %>% 
    mutate(lat = ifelse(coordinates_obscured, private_latitude, latitude),
           lon = ifelse(coordinates_obscured, private_longitude, longitude)) %>% 
    select(-longitude, -latitude, -private_longitude, -private_latitude) %>% 
    filter(place_county_name == "District of Columbia")
  inat_points <- inat_dat %>% 
    vect(geom = c("lon", "lat"), crs = "+proj=longlat",
         keepgeom = TRUE) %>% 
    project(crs(grid))
  
  inat_dat$block_ID <- extract(grid, inat_points)$block_ID
  
  inat_dat <- filter(inat_dat, !is.na(block_ID))
  
  effort_ct <- inat_dat %>% count(block_ID)
  
  inat_grid_cts <- data.frame(
    block_ID = unique(grid$block_ID)
  ) %>% left_join(effort_ct)
  
  for (i in 1:length(target_species)) {
    thisspec_obs <- inat_dat %>% 
      filter(scientific_name == target_species[i]) %>% 
      count(block_ID)
    colnames(thisspec_obs)[2] <- target_species_clean[i]
    if (nrow(thisspec_obs) == 0) stop("Spec. missing from iNat")
    
    inat_grid_cts <- left_join(inat_grid_cts, thisspec_obs)
  }
  
  inat_grid_cts[is.na(inat_grid_cts)] <- 0

  write_csv(inat_grid_cts, "intermediate/inat_grid_cts.csv")
} else {
  inat_grid_cts <- read_csv("intermediate/inat_grid_cts.csv")
}


#### Make camera metadata and detection histories ####
source("dat_helper.R")
if (!file.exists("intermediate/CT") | redo_dethist) {
  
  deployments <- ct_dat %>% 
    group_by(Deployment.ID, 
             longitude = Actual.Long, 
             latitude = Actual.Lat,
             year = year(Date),
             .groups = "drop") %>% 
    summarize(start_date = as.Date(min(Date)), end_date = as.Date(max(Date))) %>% 
    mutate(depl_length = end_date - start_date,
           deployment_id = Deployment.ID) %>% 
    filter(depl_length > 1)
  
  ct_pts <- vect(deployments, geom = c("longitude", "latitude"),
                 crs = "+proj=longlat")
  deployments$block_ID <- extract(grid, ct_pts)$block_ID
  deployments <- deployments %>% 
    filter(!is.na(block_ID))
  
  sequences <- ct_dat %>% 
    select(deployment_id = Deployment.ID, obs_time = Begin.Time, Species.Name)
  
  dethist_list <- list()
  for (i in 1:length(target_species)) {
    dethist_list[[i]] <- make_detection_hist(
      species = target_species[i], deployments, sequences, 
      detection_window_days = 1, detection_gap_days = 0)
  }
  
  saveRDS(dethist_list, "intermediate/dethist_list.RDS")
} else {
  dethist_list <- readRDS("intermediate/dethist_list.RDS")
}


