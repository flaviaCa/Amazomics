# functions which assign within host measures 
# Measures include output for treated or untreated for the following:
#   total infectivity;
#   asexual density for higher and lower sexually committing parasites;
#   gametocyte density for higher and lower sexually committing parasites;

match_infectivity_function <- function(indice){
  return(strain_df_long %>% filter(strain == out_df[it,indice,"strain"]
                                   , COI == out_df[it,indice,"COI"]
                                   ,age_infection == out_df[it,indice,"age_inf"]
                                   ,variable == "infectivity_tot" ) %>% pull(measurement))
  
}

match_infectivity_treated_function <- function(indice){
  return(strain_df_long %>% filter(strain == out_df[it,indice,"strain"]
                                   , COI == out_df[it,indice,"COI"]
                                   ,age_infection == out_df[it,indice,"age_inf"]
                                   ,variable == "infectivity_tot_treated" ) %>% pull(measurement))
  
}

match_low_dens_function <- function(indice){
  return(strain_df_long %>% filter(strain == out_df[it,indice,"strain"]
                                   , COI == out_df[it,indice,"COI"]
                                   ,age_infection == out_df[it,indice,"age_inf"]
                                   ,variable == "asex_dens_s_low" ) %>% pull(measurement))
  
}

match_high_dens_function <- function(indice){
  return(strain_df_long  %>% filter(strain == out_df[it,indice,"strain"]
                                   , COI == out_df[it,indice,"COI"]
                                   ,age_infection == out_df[it,indice,"age_inf"]
                                   ,variable == "asex_dens_s_high" ) %>% pull(measurement))
  
}

match_low_gam_function <- function(indice){
  return(strain_df_long %>% filter(strain == out_df[it,indice,"strain"]
                                   , COI == out_df[it,indice,"COI"]
                                   ,age_infection == out_df[it,indice,"age_inf"]
                                   ,variable == "gam_dens_s_low" ) %>% pull(measurement))
  
}

match_low_gam_treated_function <- function(indice){
  return(strain_df_long %>% filter(strain == out_df[it,indice,"strain"]
                                   , COI == out_df[it,indice,"COI"]
                                   ,age_infection == out_df[it,indice,"age_inf"]
                                   ,variable == "gam_dens_s_low_treated" ) %>% pull(measurement))
  
}

match_high_gam_function <- function(indice){
  return(strain_df_long %>% filter(strain == out_df[it,indice,"strain"]
                                   , COI == out_df[it,indice,"COI"]
                                   ,age_infection == out_df[it,indice,"age_inf"]
                                   ,variable == "gam_dens_s_high" ) %>% pull(measurement))
  
}

match_high_gam_treated_function <- function(indice){
  return(strain_df_long %>% filter(strain == out_df[it,indice,"strain"]
                                   , COI == out_df[it,indice,"COI"]
                                   ,age_infection == out_df[it,indice,"age_inf"]
                                   ,variable == "gam_dens_s_high_treated" ) %>% pull(measurement))
  
}
