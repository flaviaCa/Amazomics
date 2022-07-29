# ---------------------------------------  fct_select_strain.R  ---------------------------------------------------------------------------------------------------------------------

# from the within host data frame, select the strains that correspond to the specified parasite characteristics
# strain characteristics are the sexual conversion rates, the gametocyte development time, the immune cross-reactivity level

fct_select_strain <- function(out,conv_low,tau_low,tau_high,conv_high,cross.react){
  # temp give a name to the strains to be selected
  select_strains <- c(paste(conv_low,tau_low,sep="_"),paste(conv_high,tau_high,sep="_"),"0_0")

  # create a new data-frame with the strains of interests
  strain_df <- out %>% select(time, z_1,z_2,tau_1,tau_2, G1,G2,infectivity_tot
                              , G1__tmt_G.eff_0, G2__tmt_G.eff_0 
                              , infectivity_tot__tmt_G.eff_0, CR.factor, age_treat) %>%
    
    unite("temp_strain1", c(z_1,tau_1), remove = FALSE) %>% unite("temp_strain2", c(z_2,tau_2), remove = FALSE) %>%
    filter(temp_strain1 %in% select_strains) %>%   filter(temp_strain2 %in% select_strains) %>% 
    
    # change time into age of infection, starting at 1 and not 0
    mutate(age_infection = time + 1) %>% select(-time) %>% mutate(age_treatement = age_treat + 1) %>%
    
    # remove the empty rows where infectivity ended
    filter(!is.na(infectivity_tot) | age_infection < 20) %>% # added time so that we keep the first days where infected but not yet infectious
    
    # add a column specifying if it is a co-infection (polyclonal or monoclonal) or single infection
    mutate(COI = ifelse(z_2 == 0,"clonal","coinfect")) %>%
    # add temp column
    mutate(temp_order = case_when(
      z_1 == 0 & z_2 == conv_low & tau_2 == tau_low ~ "not_low",
      z_1 == conv_low & z_2 == 0 & tau_1 == tau_low  ~ "low_not",
      z_1 == 0 & z_2 == conv_high  & tau_2 == tau_high ~ "not_high",
      z_1 == conv_high & z_2 == 0  & tau_1 == tau_high ~ "high_not",
      
      z_1 == conv_low & z_2 == conv_low & tau_1 == tau_low   & tau_2 == tau_low ~ "low_low",
      z_1 == conv_high & z_2 == conv_high & tau_1 == tau_high   & tau_2 == tau_high ~ "high_high",
      
      z_1 == conv_low & z_2 == conv_high & tau_1 == tau_low   & tau_2 == tau_high ~ "low_high",
      z_1 == conv_high & z_2 == conv_low & tau_1 == tau_high   & tau_2 == tau_low ~ "high_low",
      
      TRUE                      ~ NA_character_
    )) %>%
    
    # add strain names
    mutate(strain = case_when(
      temp_order %in% c("low_low","low_not","not_low") ~ "S_low",
      temp_order %in% c("high_high","high_not","not_high") ~ "S_high",
      temp_order %in% c("high_low","low_high") ~ "S_mix",
      TRUE                      ~ NA_character_
    )) %>%
    
    # add gametocyte density for both genotypes
    mutate(gam_dens_s_low = case_when(
      COI == "clonal" & temp_order == "low_not"  ~ G1,
      COI == "clonal" & temp_order == "not_low"  ~ G2,
      COI == "coinfect" & temp_order == "low_high"  ~ G1,
      COI == "coinfect"  & temp_order == "low_low"   ~ G1 +G2,
      COI == "coinfect"  & temp_order == "high_low"  ~ G2,
      TRUE                      ~ NA_real_
    )) %>%
    
    mutate(gam_dens_s_high = case_when(
      COI == "clonal"  & temp_order == "high_not" ~ G1,
      COI == "clonal"  & temp_order == "not_high" ~ G2,
      
      COI == "coinfect"   & temp_order == "high_high" ~ G1 + G2,
      COI == "coinfect"   & temp_order == "low_high"   ~ G2,
      COI == "coinfect"   & temp_order == "high_low"  ~ G1,
      
      TRUE                      ~ NA_real_
    ))  %>%
    
    # add gametocyte densities for both genotypes when treated
    mutate(gam_dens_s_low_treated = case_when(
      COI == "clonal"& temp_order == "low_not"   ~  G1__tmt_G.eff_0,
      COI == "clonal"& temp_order == "not_low"   ~  G2__tmt_G.eff_0,
      
      COI == "coinfect" & temp_order == "low_high"   ~  G1__tmt_G.eff_0,
      COI == "coinfect" & temp_order == "low_low"    ~ G1__tmt_G.eff_0 +  G2__tmt_G.eff_0,
      COI == "coinfect"& temp_order == "high_low"~ G2__tmt_G.eff_0,
      TRUE                      ~ NA_real_
    )) %>%
    
    mutate(gam_dens_s_high_treated = case_when(
      COI == "clonal" & temp_order == "high_not" ~ G1__tmt_G.eff_0,
      COI == "clonal"  & temp_order == "not_high" ~ G2__tmt_G.eff_0,
      
      COI == "coinfect" & temp_order == "high_high" ~  G1__tmt_G.eff_0 + G2__tmt_G.eff_0,
      COI == "coinfect"  & temp_order == "low_high"   ~ G2__tmt_G.eff_0,
      COI == "coinfect"  & temp_order == "high_low"  ~  G1__tmt_G.eff_0,
      TRUE                      ~ NA_real_
    ))  %>%
    
    rename(infectivity_tot_treated = infectivity_tot__tmt_G.eff_0)
  
  strain_df <- strain_df %>% select(-c(G1,G2,G1__tmt_G.eff_0,G2__tmt_G.eff_0,temp_order,temp_strain1,temp_strain2,age_treat))
  return(strain_df)}
