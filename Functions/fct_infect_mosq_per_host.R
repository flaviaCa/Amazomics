# ---------------------------------------  function_infect_mosq_per_host.R  ---------------------------------------------------------------------------------------------------------------------

# Assign infected mosquitos to available host to infect.
# If more than one mosquito is infecting the same host, it defines it as a super-infection
# input 1 new_host: hosts to be infected 
# input 2 infection_df: infected mosquitos ready to infect human hosts in this iteration
# output is a data frame with rows: host IDs to infect, columns: infection characteristics (super-infection or single infection, polyclonal or monoclonal of which gentoype) 

fct_infect_mosq_per_host <- function(new_host,infection_df){
  
  # create a vector of all mosquito available for the selected hosts
  all_mosq <- data.frame(matrix(ncol = 3, nrow = sum(as.numeric(new_host$n_mosq))))
  colnames(all_mosq) <- c("strain", "host","transmission")
  # assign host host ID for each mosquito
  all_mosq$host <-   unlist(apply(new_host,1,function(x) rep(x["Host_ID"],as.numeric(x["n_mosq"]))))
  # randomly assign new infections to mosquito, from the pool of available mosquito among the hosts:
  all_mosq$strain[sample(1:nrow(all_mosq),min(nrow(all_mosq),nrow(infection_df)),replace=F)] <- infection_df$strain[sample(1:nrow(infection_df),min(nrow(all_mosq),nrow(infection_df)),replace=F)]
  
  # keep only the mosquitos that carry an infection,
  all_mosq$transmission <- 1 
  all_mosq <- as.data.frame(all_mosq[!is.na(all_mosq$strain),] %>% 
                              # change to wide data frame with rows each host to infect, and columns (n=3) the number of mosquito infected with either S_high, Slow, S_mix
                              group_by(strain,host) %>% summarise(transmission = n(), .groups = 'drop')%>% spread(strain,transmission))
  
  
  all_mosq$strain <- "to_det"
  all_mosq$COI <- "clonal"

  # check all strain are as data.frame
  if(!("S_high" %in% colnames(all_mosq))){
    all_mosq <- all_mosq %>% mutate(S_high = NA_real_) 
  }
  if(!("S_low" %in% colnames(all_mosq))){
    all_mosq <- all_mosq %>% mutate(S_low = NA_real_) 
  }
  if(!("S_mix" %in% colnames(all_mosq))){
    all_mosq <- all_mosq %>% mutate(S_mix = NA_real_) 
  }
  
  # assign if the resulting infection in the new host will be monoclonal (S_high or S_low) or polyclonal (S_mix)
  all_mosq$strain[all_mosq$S_high>=1 & all_mosq$S_low>=1] <- "S_mix"
  all_mosq$strain[all_mosq$S_high>=1 & is.na(all_mosq$S_low) & is.na(all_mosq$S_mix)] <- "S_high"
  all_mosq$strain[all_mosq$S_low>=1 & is.na(all_mosq$S_high) & is.na(all_mosq$S_mix)] <- "S_low" 
  all_mosq$strain[all_mosq$S_mix>=1] <- "S_mix"

  # determine if the each host has a single infection (ie clonal only one mosquito) or a super-infection (ie coinfect, at least two mosquito infecting the same host)
  all_mosq$COI[all_mosq$strain=="S_mix"] <- "coinfect"
  all_mosq$COI[all_mosq$strain=="S_mix" & is.na(all_mosq$S_high) & is.na(all_mosq$S_low) & all_mosq$S_mix==1] <- "clonal"
  all_mosq$COI[all_mosq$S_high>1 & all_mosq$strain=="S_high"] <- "coinfect"
  all_mosq$COI[all_mosq$S_low>1 & all_mosq$strain=="S_low"] <- "coinfect"

  return(all_mosq) 
}