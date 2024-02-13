library(tidyverse)
library(terra)
library(mgcv)

# Choose a projection (AEA)
main_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# Read in iNat data
inat_dat <- readxl::read_xlsx("data/DC_mammal_iNat_data_1262023.xlsx")
inat_dat <- inat_dat %>% 
  filter(!coordinates_obscured | !is.na(private_latitude),
         is.na(positional_accuracy) | positional_accuracy <= 1000) %>% 
  mutate(lat = ifelse(coordinates_obscured, private_latitude, latitude),
         lon = ifelse(coordinates_obscured, private_longitude, longitude)) %>% 
  select(-longitude, -latitude, -private_longitude, -private_latitude) %>% 
  filter(place_county_name == "District of Columbia")

inat_speccts <- inat_dat %>% 
  count(scientific_name) %>% 
  arrange(-n)


# Read in camera trap data
specs_to_drop <- c(
  "Camera Misfire", "Homo sapiens", "Other Bird species",
  "Vehicle", "Camera Trapper", "Unknown Ground-Dove", "Turdus migratorius",
  "Unknown Animal", "No Animal", "Bicycle", "Unknown Bird", "Cyanocitta cristata",
  "Meleagris gallopavo", "Corvus brachyrhynchos", "Colaptes auratus",
  "Coragyps atratus", "Cathartes aura", "Reptile species", "Owl species",
  "Animal Not on List", "Raptor Species", "Melanerpes erythrocephalus"
)

ct_dat <- readxl::read_xlsx("data/final_DCCC_data.xlsx", sheet = 2) %>% 
  filter(Actual.Long > -77.3) %>% 
  filter(!Species.Name %in% specs_to_drop)
ct_speccts <- ct_dat %>% 
  count(Species.Name, Common.Name) %>% 
  arrange(-n) 

# Encounter rates at each coordinate
ct_st_counts <- ct_dat %>% 
  mutate(year = year(Date), month = month(Date)) %>% 
  count(Actual.Lat, Actual.Long, #year, month, 
        Species.Name) %>% 
  group_by(Actual.Lat, Actual.Long, #year, month
           ) %>% 
  mutate(n_encounters = sum(n), this_enc_rate = n / n_encounters)

USMap <- geodata::gadm(country = "USA", path = "../example_data_combination/data/geodata")
dcMap <- USMap[USMap$NAME_1 == "District of Columbia",] %>% 
  project(main_crs)



#### Align datasets on a grid #####

# Make a spatial dataset for ct data
grid <- rast(ext(dcMap), resolution = 2500)
values(grid) <- sample(1:ncell(grid))
crs(grid) <- main_crs

ct_points <- ct_dat %>% 
  vect(geom = c("Actual.Long", "Actual.Lat"), crs = "+proj=longlat",
       keepgeom = TRUE) %>% 
  project(main_crs)

ct_points$cell <- extract(grid, ct_points)$lyr.1

ct_convex_hull <- terra::convHull(ct_points)
grid_pts <- as.points(grid)
grid_pts$in_convhull <- !is.nan(terra::extract(ct_convex_hull, grid_pts)[,2])
grid_pts <- grid_pts[grid_pts$in_convhull, ]

ct_cellcounts <- ct_points %>% 
  as.data.frame() %>% 
  count(cell, Species.Name)
ct_effort <- ct_points %>% 
  as.data.frame() %>% 
  count(cell) %>% 
  rename(effort = n)


inat_points <- inat_dat %>% 
  vect(geom = c("lon", "lat"), crs = "+proj=longlat",
       keepgeom = TRUE) %>% 
  project(main_crs)

inat_points$cell <- extract(grid, inat_points)$lyr.1

inat_cellcounts <- inat_points %>% 
  as.data.frame() %>% 
  count(cell, scientific_name) %>% 
  filter(!is.na(cell))
inat_effort <- inat_points %>% 
  as.data.frame() %>% 
  count(cell) %>% 
  filter(!is.na(cell)) %>% 
  rename(effort = n)


cell_coords <- as.data.frame(grid, xy = TRUE)
colnames(cell_coords)[3] <- "cell"

# Filter to cells with any effort
cell_coords <- cell_coords %>% 
  filter(cell %in% c(inat_effort$cell, ct_effort$cell))

# Filter to cells with sampling in both
cells_wboth <- cell_coords %>% 
  filter(cell %in% inat_effort$cell, cell %in% ct_effort$cell)


#### Fit a GAM for one species ####

fit_gam <- function(species) {
  
  inat_cts_thisspec <- left_join(
    inat_effort, by = "cell",
    inat_cellcounts %>%  filter(scientific_name == species)
  ) %>% 
    mutate(nsuccess = ifelse(is.na(n), 0, n),
           nfail = effort - nsuccess) %>% 
    left_join(cell_coords, by = "cell")
  ct_cts_thisspec <- left_join(
    ct_effort, by = "cell",
    ct_cellcounts %>%  filter(Species.Name == species)
  ) %>% 
    mutate(nsuccess = ifelse(is.na(n), 0, n),
           nfail = effort - nsuccess) %>% 
    left_join(cell_coords, by = "cell")
  
  
  inat_gam <- gam(cbind(nsuccess, nfail) ~ te(x, y, k = 5),
                  data = inat_cts_thisspec, 
                  family = quasibinomial(link = "logit"),
                  optimizer = c("outer", "bfgs"),
                  control = list(maxit = 100000))
  ct_gam   <- gam(cbind(nsuccess, nfail) ~ te(x, y, k = 5),
                  data = ct_cts_thisspec, 
                  family = quasibinomial(link = "logit"),
                  optimizer = c("outer", "bfgs"),
                  control = list(maxit = 100000))
  
  check_ct   <- k.check(ct_gam)
  check_inat <- k.check(inat_gam)
  
  # inat_pred <- predict(inat_gam, newdata = cell_coords, se.fit = TRUE)
  # ct_pred <-   predict(ct_gam, newdata = cell_coords, se.fit = TRUE)
  
  rmvn <- function(n,mu,sig) { ## MVN random deviates
    L <- mroot(sig);m <- ncol(L);
    t(mu + L%*%matrix(rnorm(m*n),m,n)) 
  }
  
  Xp_ct    <- predict(ct_gam, newdata = cells_wboth, type="lpmatrix") 
  Xp_inat  <- predict(inat_gam, newdata = cells_wboth, type="lpmatrix") 
  # Xp_ct    <- predict(ct_gam, newdata = cell_coords, type="lpmatrix") 
  # Xp_inat  <- predict(inat_gam, newdata = cell_coords, type="lpmatrix") 
  
  br_ct   <- rmvn(10000, coef(ct_gam),ct_gam$Vp) ## 1000 replicate param. vectors
  br_inat <- rmvn(10000, coef(inat_gam),inat_gam$Vp) ## 1000 replicate param. vectors
  
  res1 <- rep(0,10000)
  res2 <- rep(0,10000)
  res3 <- rep(0,10000)
  for (i in 1:10000) { 
    pr_ct   <- Xp_ct   %*% br_ct[i,] ## replicate predictions
    pr_inat <- Xp_inat %*% br_inat[i,] ## replicate predictions
    res1[i] <- median(pr_ct) ## median eButterfly prediction
    res2[i] <- median(pr_inat) ## median iNat prediction
    res3[i] <- median(pr_inat - pr_ct) ## median difference
  }
  
  
  target_df <- data.frame(
    species = species, 
    inat_ct = sum(inat_cts_thisspec$nsuccess),
    ct_ct = sum(ct_cts_thisspec$nsuccess)
  )
  
  target_df$difference   <- median(res3)
  target_df$diff_SE <- sd(res3)
  
  target_df$inat_medPred <- median(res2)
  target_df$inat_medResp <- nimble::expit(median(res2))
  target_df$inat_SE <- sd(res2)
  target_df$ct_medPred <- median(res1)
  target_df$ct_medResp <- nimble::expit(median(res1))
  target_df$ct_SE <- sd(res1)
  
  target_df$ct_GAM_pval <- min(check_ct[, 4])
  target_df$inat_GAM_pval <- min(check_inat[, 4])
  
  return(target_df)
}


target_species <- unique(ct_speccts$Species.Name[ct_speccts$n > 100])
target_species <- target_species[!target_species %in% c(
  "Canis familiaris", "Equus caballus", "Unknown Small Rodent"
)]

pb <- progress::progress_bar$new(total = length(target_species))
result_list <- list()

for (i in 1:length(target_species)) {
  result_list[[i]] <- fit_gam(species = target_species[i])
  pb$tick()
}

result_df <- bind_rows(result_list)



#### Visualize results ####

result_df$LB <- result_df$difference - 1.96 * result_df$diff_SE
result_df$UB <- result_df$difference + 1.96 * result_df$diff_SE
result_df <- result_df %>% mutate(
  sig_type = ifelse(
    sign(LB) != sign(UB), "Non-significant",
    ifelse(LB > 0, "Overreported", "Underreported")
  )
)


output_table <- result_df %>% 
  mutate(CI_text = paste0(round(difference, digits = 2), "(", 
                          round(LB, digits = 2), ", ",
                          round(UB, digits = 2), ")")) %>% 
  select(species, CT_count = ct_ct, iNat_count = inat_ct, Finding = sig_type, 
         OR_index = CI_text)



output_table


#### Vis (exploratory) ####
ct_st_counts %>% 
  filter(Species.Name %in% ct_speccts$Species.Name[ct_speccts$n > 100]) %>%
  ggplot() +
  tidyterra::geom_spatvector(data = dcMap) +
  geom_point(aes(Actual.Long, Actual.Lat, color = this_enc_rate)) +
  theme_void() +
  facet_wrap(~Species.Name) +
  scale_color_viridis_c()



ct_points %>% 
  ggplot() +
  tidyterra::geom_spatvector(data = dcMap) +
  tidyterra::geom_spatvector(data = as.polygons(grid), fill = NA) +
  tidyterra::geom_spatvector(aes(color = cell)) 

ggplot() +
  tidyterra::geom_spatvector(data = dcMap) +
  tidyterra::geom_spatvector(data = inat_points) 
  # theme_void()






# Draft CT gam
mod_data <- ct_st_counts %>% filter(Species.Name == "Sciurus carolinensis")
fit <- mgcv::gam(cbind(n, n_encounters - n) ~ te(Actual.Long, Actual.Lat),
          data = ct_st_counts, family = binomial)





cell_coords$inat_pred <- inat_pred
cell_coords$ct_pred <- ct_pred
cell_coords %>% 
  ggplot() +
  geom_tile(aes(x, y, fill = ct_pred)) +
  tidyterra::geom_spatvector(data = dcMap, fill = NA) +
  scale_fill_viridis_c()
