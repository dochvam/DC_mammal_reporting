fill_detection_history <- function(empty_df, this_sequences) {
  filled_df <- empty_df
  filled_df$observed <- FALSE
  filled_df$yday_start <- yday(filled_df$obs_start_date)
  filled_df$yday_end <- yday(filled_df$obs_end_date)
  
  this_sequences <- this_sequences %>% 
    filter(deployment_id %in% filled_df$deployment_id)
  
  indat <- logical(nrow(this_sequences))
  for (i in 1:nrow(this_sequences)) {
    filled_df$observed[
      filled_df$deployment_id == this_sequences$deployment_id[i] & 
        filled_df$yday_start <= this_sequences$obs_yday[i] &
        filled_df$yday_end > this_sequences$obs_yday[i]
    ] <- TRUE
    indat[i] <- any(
      filled_df$deployment_id == this_sequences$deployment_id[i] & 
        filled_df$yday_start <= this_sequences$obs_yday[i] &
        filled_df$yday_end > this_sequences$obs_yday[i]
    )
  }
  
  return(filled_df)
}


make_empty_deployment_df <- function(deployments, detection_window_days,
                                     detection_gap_days) {
  
  df_list <- list()
  
  for (i in 1:nrow(deployments)) {
    
    this_date <- deployments$start_date[i]
    date_seq_start <- Date()
    date_seq_end <- Date()
    date_seq_len <- numeric()
    while(this_date < deployments$end_date[i]) {
      date_seq_start <- c(date_seq_start, this_date)
      date_seq_end <- c(date_seq_end, min(this_date + detection_window_days,
                                          deployments$end_date[i]))
      
      this_date <- this_date + detection_window_days + detection_gap_days
    }
    date_seq_len <- date_seq_end - date_seq_start
    
    df_list[[i]] <- data.frame(
      deployment_id = deployments$deployment_id[i],
      block_ID = deployments$block_ID[i],
      longitude = deployments$longitude[i],
      latitude = deployments$latitude[i],
      year = deployments$year[i],
      obs_start_date = date_seq_start,
      obs_end_date = date_seq_end,
      obs_len = date_seq_len
    )
  }
  
  return(bind_rows(df_list))
}


make_detection_hist <- function(species, 
                                deployments,
                                sequences,
                                detection_window_days, 
                                detection_gap_days) {
  
  this_sequences <- sequences %>% 
    filter(Species.Name == species) %>% 
    mutate(obs_hour = hour(obs_time),
           obs_year = year(obs_time),
           obs_yday = ifelse(
             obs_hour < 12, yday(obs_time), yday(obs_time) + 1
           ))
  
  
  this_detection_history <- 
    make_empty_deployment_df(deployments, detection_window_days, 
                             detection_gap_days) %>% 
    mutate(obs_yday = yday(obs_start_date)) %>% 
    fill_detection_history(this_sequences)
  
  return(this_detection_history)
}