# ---------------------------------------  fct_transmission_host_to_mosq.R  ---------------------------------------------------------------------------------------------------------------------

# for an infectious host, calculate the number of infections transmitted to mosquitos at a given time step
# inputs are infection characteristics for a given host at current time-step:
#   total infectiousness;
#   gametocyte density for each genotype;
#   number of available mosquitos to bite the infectious host.
#
# other inputs are:
#   the data frame in which the new infection should be stored;
#   normalization factor making sure that the total number of new infections (across all hosts) are kept to a nearly constant level when no treatment is implemented
#   number of days assumed between a mosquito getting infected and a mosquito being infectious.
#
# output: updated data frame storing the new infections, specifying the type of infection (S_high, S_low, or S_mix) and the number of days before it will infect new human hosts


fct_transmission_host_to_mosq <- function(new_inf,infect_tot,G1,G2,lag_days=17,norm_fact,n_mosq){
  # draw the number of infections resulting from this host, in function of infectiousness, number of biting mosquitos, and a normalizing factor
  transmission <- rpois(n=1,lambda=infect_tot*norm_fact*n_mosq)
  if(transmission>0){
    for(t in 1:transmission){
      # determine how many parasites from one or the other genotype are passed to the mosquito
      # there is a shortcut here, as we directly count how many of the (n_spz=10) sporozoites will be passed from the mosquito to the newly infected human 
      trans_s1 <- rbinom(size=n_spz,n=1,p = G1/(G1+G2))
      trans_s2 <- n_spz-trans_s1
      # according to the parasites that will be transmitted, determine if the infection is monoclonal S_high or S_low, or polyclonal S_mix
      if(trans_s1>0&trans_s2==0){
        name_strain <- "S_low"
      }else{
        if(trans_s1==0&trans_s2>0){
          name_strain <- "S_high"
        }else{
          if(trans_s1>0&trans_s2>0){
            name_strain <- "S_mix"
          }
        }
      }
      # store the genotype(s) present in each infected host
      new_inf$strain[new_inf$timer<= -1000][1] <- name_strain
      # set the time to t=lag_days to know when these infections will be passed to to a new human host
      new_inf$timer[new_inf$timer<= -1000 & new_inf$strain == name_strain][1] <- lag_days
       }
  }
  return(new_inf)
}
  
    