# ---------------------------------------  main_Transmission_model.R  ---------------------------------------------------------------------------------------------------------------------
# Runs the transmission model - outputs the the number of infections per genotype at each time step
#
#
#
# ---------------------------------------  Store strains of interest  ---------------------------------------------------------------------------------------------------------------------
# select strain of interest from the within host model output
strain_df <- fct_select_strain(out_all
                               ,conv_low = conv_low,conv_high = conv_high
                               ,tau_low = tau1,tau_high = tau2
                               ,cross.react = cross.react)

# remove duplicates and unnecessary columns
strain_df <- strain_df %>% select(-c(z_1,z_2,tau_1,tau_2)) %>% distinct()

# Add replicate within host dynamics for S_mix result from co-transmission:
temp <- strain_df[strain_df$strain=="S_mix",]
temp$COI <- "clonal"

strain_df <- rbind(strain_df,temp)
rm(temp)

# Cap infections to a maximum to keep them all with same length
max_post_inf  <- strain_df %>% group_by(strain,COI) %>% select(age_infection, strain,COI) %>% summarise(max = max(age_infection),.groups="drop")
keep_min_length <- min(max_post_inf$max) 
max_post_inf$max <- keep_min_length
strain_df <- strain_df %>% filter(age_infection <= keep_min_length)

# convert to long data frame 
strain_df_long <- gather(strain_df, variable, measurement, c(gam_dens_s_low,gam_dens_s_high ,infectivity_tot,age_treatement,gam_dens_s_low_treated,gam_dens_s_high_treated,infectivity_tot_treated)) 

# ---------------------------------------  Determine initial infections  ---------------------------------------------------------------------------------------------------------------------

# determine how many hosts with which genotype
n_per_strain <- fq_strain
n_per_strain$n_per_strain <- as.integer(round(H*n_per_strain$fq_strain))

# If the total number of infections, as the sum of infections of different genotype according to defined frequency, is not exactly equal the number of host: 
if(sum(n_per_strain$n_per_strain)!=H){
  # change the number of infected host with the dominant strain 
  dom_inf <- as.numeric(row.names(n_per_strain)[n_per_strain$n_per_strain==max(n_per_strain$n_per_strain)])
  n_prev_inf <- n_per_strain$n_per_strain[dom_inf]
  n_per_strain$n_per_strain[dom_inf] <- H - sum(n_per_strain$n_per_strain[-dom_inf])
  print(paste("modified number of intitial infections from",n_prev_inf,"to",
              n_per_strain[dom_inf,"n_per_strain"],"for hosts infected with",n_per_strain[dom_inf,"COI"],n_per_strain[dom_inf,"strain"],"to match total number of infected hosts", sep=" "))
}

# to equally distribute the initial infection age, determine how many infections of each age to assign
n_s <- merge(max_post_inf, n_per_strain,by=c("strain","COI")) %>% mutate(n_s_age = round(n_per_strain/max))
# add the at what age of infection the host gets treated
n_s <- merge(n_s,strain_df %>% group_by(strain,COI) %>% select(strain,COI,age_treatement) %>% unique())

# ---------------------------------------  Compute normalizing factor  ---------------------------------------------------------------------------------------------------------------------

# Expected number of cleared infections at each time step (when population is untreated)
E_ended_inf <- round(H/keep_min_length)

# Increase infectivity by assumed fraction super-infection
FOI <- E_ended_inf*(1+p_superInf)

# ---------------------------------------  Initialize infections  ---------------------------------------------------------------------------------------------------------------------
it <- 1

# initialize a matrix which will store all infection measures for each host at each time step
out_df <- array(data = NA,
                dim= c(Time, I+1, 14),
                dimnames = list(c(1:Time), 
                                c("time_step",1:(I)),
                                c("gam_dens_s_low","gam_dens_s_high","infectivity_tot","infectivity_tot_noTmt",
                                  "gam_dens_s_low_treated","gam_dens_s_high_treated","infectivity_tot_treated",
                                  "strain","age_inf","treat","COI","age_treatment","n_mosq","Host_ID")))

# initialize all hosts to be untreated (0)
out_df[it,2:I,"treat"] <- 0 

# assign a host ID to each host
out_df[it,2:(H+1),"Host_ID"] <- paste("ID_",seq(1:H),sep="") 

# add a number of mosquito for each host
out_df[it,2:I,"n_mosq"] <- rnbinom(2:I, k, 1-(A/(A+k)))

# make sure there is at least 1 mosquito per host
out_df[it,2:I,"n_mosq"][out_df[it,2:I,"n_mosq"]=="0"] <- 1

# add all time steps
out_df[,"time_step",] <- c(1:Time) 

# allocate age of infection in each host
out_df[it,2:(H+1),"age_inf"] <-  unlist(apply(n_s,1,function(x) rep_len(1:x["max"],length.out = x["n_per_strain"])))

# allocate at what age of infection the host gets treated
out_df[it,2:(H+1),"age_treatment"] <- unlist(apply(n_s,1,function(x)  rep_len(x["age_treatement"],length.out = x["n_per_strain"])))

# allocate COI status
out_df[it,2:(H+1),"COI"] <- unlist(apply(n_s,1,function(x) rep_len(x["COI"],length.out = x["n_per_strain"])))

# allocate the name of the genotype each host is infected with
out_df[it,2:(H+1),"strain"] <- unlist(apply(n_s,1,function(x) rep_len(x["strain"],length.out = x["n_per_strain"])))

# allocate infectiousness
out_df[it,2:(H+1),"infectivity_tot"] <- sapply(2:(H+1),match_infectivity_function) 
out_df[it,2:(H+1),"infectivity_tot_noTmt"] <- out_df[it,2:(H+1),"infectivity_tot"] 
out_df[it,2:(H+1),"infectivity_tot_treated"] <- sapply(2:(H+1),match_infectivity_treated_function) 

# allocate gametocyte density of each genotype
out_df[it,2:(H+1),"gam_dens_s_low"] <- sapply(2:(H+1),match_low_gam_function) 
out_df[it,2:(H+1),"gam_dens_s_high"] <- sapply(2:(H+1),match_high_gam_function) 

# allocate gametocyte density of each genotype if treated
out_df[it,2:(H+1),"gam_dens_s_low_treated"] <- sapply(2:(H+1),match_low_gam_treated_function) 
out_df[it,2:(H+1),"gam_dens_s_high_treated"] <- sapply(2:(H+1),match_high_gam_treated_function) 

# first iteration is now stored in dataframe

# ---------------------------------------  Initialize infected mosquito before start of simulation  ---------------------------------------------------------------------------------------------------------------------

# account for all infections that happened until 15 days prior the infection, so that we can add them once we start simulation

# data.frame to store new infections:
new_inf <- data.frame(matrix(ncol = 2, nrow = A*I))
colnames(new_inf) <- c("timer", "strain")
new_inf$timer <- -1000
new_inf$strain <- -1000

for(it in -15:1){
  # take all infections from first day of simulation backward
  prev_age_inf <- as.numeric(out_df[1,2:(H+1),"age_inf"]) + (it-1)
  # add NA's manually so that an infection of age 0 doesn't get pushed forward to 1
  prev_age_inf[prev_age_inf<=0] <- NA
  
  # need to take infections backward unless the infection has stopped at previous time step or there was no infection
  infection_backward <- left_join(data.frame("strain" = out_df[1,2:(H+1),"strain"]
                                             ,"COI" = out_df[1,2:(H+1),"COI"]
                                             ,"age_infection" = prev_age_inf
                                             ,"n_mosq" = out_df[1,2:(H+1),"n_mosq"]
                                             ,"Host_ID" = out_df[1,2:(H+1),"Host_ID"]),
                                  strain_df[,c("infectivity_tot","gam_dens_s_low" ,"gam_dens_s_high","strain","age_infection","COI")],by = c("strain", "age_infection","COI"))
  
  # select hosts that are infectious
  n_infectious <- as.data.frame(infection_backward[,c("infectivity_tot","n_mosq","Host_ID")]) %>% 
    filter(
      !is.na(infectivity_tot) & infectivity_tot > 0 &
        !is.na(n_mosq) & n_mosq>0 )  
  who_infectious <- n_infectious$Host_ID
  
  # calculate the total infectious potential, accross all infectious hosts 
  n_infectious$infectious_potential <- apply(n_infectious, 1, function(x) as.numeric(x["infectivity_tot"])*as.numeric(x["n_mosq"]))
  infectious_potential <- sum(n_infectious$infectious_potential,na.rm=T)
  
  # new infections for each of infectious hosts
  for(h in who_infectious){
    # host's infectiousness of host
    infect <-  as.numeric(infection_backward$infectivity_tot[infection_backward$Host_ID==h])
    # host's gametocyte density of each genotype
    G1 <- as.numeric(infection_backward$gam_dens_s_low[infection_backward$Host_ID==h])
    G2 <- as.numeric(infection_backward$gam_dens_s_high[infection_backward$Host_ID==h])
    G1[is.na(G1)] <- 0
    G2[is.na(G2)] <- 0
    
    # host's number of biting mosquitos
    n_mosq <- as.numeric(infection_backward$n_mosq[infection_backward$Host_ID==h])
    
    # calculate the number of mosquitos infected by the host - and specify by which genotype (S_high, S_low or S_mix)
    new_inf <- fct_transmission_host_to_mosq(new_inf,infect_tot=infect,G1=G1,G2=G2,norm_fact = (FOI/infectious_potential), n_mosq = n_mosq)
  }
  
  # remove one time step to the timer of the new infections
  new_inf$timer <- new_inf$timer-1
  
  rm(list=c("infect","G1","G2","infectious_potential","who_infectious",
            "n_infectious","infection_backward","prev_age_inf"))
}

# ---------------------------------------  SIMULATION - run iterations 2 to end of simulation time  ---------------------------------------------------------------------------------------------------------------------

host_fev_th <- numeric() 

# keep track if we need to add hosts
ad_host <- 1

# keep track of the average number of infectious mosquito per infected host, and proportion super-infection in the population
track_mosq <- data.frame(matrix(ncol = 3, nrow = Time))
colnames(track_mosq) <- c("time_step", "p_super_inf","average_mosq_host")

# run simulation
for(it in 2:Time){ 
  # CHECK MINIMUM GENOTYPE FREQUENCY
  # ensure that none of the genotypes go extinct (at least 1 percent of each genotype present as monoclonal infection)
  min_thresh <- round(length(na.omit(out_df[it-1,-1,"strain"]))*0.01,digit = 0)
  n_s1 <- length(na.omit(out_df[it-1,-1,"strain"][out_df[it-1,-1,"strain"]=="S_low"]))
  n_s2 <- length(na.omit(out_df[it-1,-1,"strain"][out_df[it-1,-1,"strain"]=="S_high"]))
  
  if(n_s1 < min_thresh){
    # add enough infections to be above the threshold
    new_inf$strain[new_inf$timer<= -1000][1:(min_thresh-n_s1)] <- rep("S_low",times=min_thresh-n_s1)
    new_inf$timer[new_inf$timer<= -1000][1:(min_thresh-n_s1)] <- rep(0,times=min_thresh-n_s1)
    print(paste("it :",it," ,forced introduction of ", min_thresh-n_s1, " S_low strains",sep=""))
  }
  
  if(n_s2 < min_thresh){
    # add enough infections to be above the threshold
    new_inf$strain[new_inf$timer<= -1000][1:(min_thresh-n_s2)] <- rep("S_high",times=min_thresh-n_s2)
    new_inf$timer[new_inf$timer<= -1000][1:(min_thresh-n_s2)] <- rep(0,times=min_thresh-n_s2)
    print(paste("it :",it," ,forced introduction of ", min_thresh-n_s2, " S_high strains",sep=""))
  }
  
  # SELECT TREATED HOSTS
  # Select hosts which get treated at this iteration:
  if(P_tmt>0){
    # Select hosts who reached infection-age of treatment at previous iteration
      host_fev_th <- as.data.frame(out_df[it-1,-1,c("strain","age_inf","treat","COI","age_treatment","Host_ID")]) %>% filter(age_inf==age_treatment)%>% tibble::rownames_to_column() 
      host_fev_th <- as.vector(host_fev_th$Host_ID)
    if(length(host_fev_th)!=0){ # if infections have reached tmt threshold in previous iteration
      # select a random number of infections based on P_tmt=fraction of treated infections
      to_treat <- sample(host_fev_th,min(length(host_fev_th),rpois(1,length(host_fev_th)*P_tmt))) 
      # tag those hosts as getting treated now (treated=1)
      out_df[it-1,!is.na(out_df[it-1,,"Host_ID"])&out_df[it-1,,"Host_ID"]%in%to_treat,"treat"] <- 1 
      host_fev_th <- vector() # empty for next iteration
      rm(to_treat)
    }
  }
  
  # UPDATE INFECTION STATUS 
  # take all infections from it - 1 one step forward
  prev_age_inf <- as.numeric(out_df[it-1,-1,"age_inf"])
  # add NA's manually so that an infection of age 0 doesn't get pushed forward to 1
  prev_age_inf[prev_age_inf==0] <- NA
  # host ID
  out_df[it,-1,"Host_ID"] <- out_df[it-1,-1,"Host_ID"]
  
  # need to take infections forward unless the infection has stopped at previous time step or there was no infection
  infection_forward <-left_join(data.frame("strain" = out_df[it-1,-1,"strain"]
                                           ,"COI" = out_df[it-1,-1,"COI"]
                                           ,"age_infection" = (prev_age_inf + 1)
                                           ,"age_treatment" = as.numeric(out_df[it-1,-1,"age_treatment"])
                                           ,"tmt_status" = as.numeric(out_df[it-1,-1,"treat"])
                                           ,"n_mosq" = out_df[it-1,-1,"n_mosq"]
                                           ,"Host_ID" = out_df[it,-1,"Host_ID"]),max_post_inf
                                ,by=c("strain","COI"))
  
  infection_forward <- left_join(infection_forward,strain_df[,c("gam_dens_s_low","gam_dens_s_low_treated","gam_dens_s_high","gam_dens_s_high_treated","infectivity_tot","infectivity_tot_treated","strain","age_infection","COI")],by = c("strain", "age_infection","COI"))
  
  infection_forward <- infection_forward %>% 
    # add a column for untreated
    rename(infectivity_tot_noTmt = infectivity_tot) %>%
    # change infectivity and gametocytemia to levels according whether the host is treated or not
    mutate(infectivity_tot = ifelse(tmt_status == 1, infectivity_tot_treated, infectivity_tot_noTmt)) %>%
    mutate(gam_dens_s_low = ifelse(tmt_status == 1, gam_dens_s_low_treated, gam_dens_s_low)) %>%
    mutate(gam_dens_s_high = ifelse(tmt_status == 1, gam_dens_s_high_treated, gam_dens_s_high)) 
  
  # put NA where infection have ended
  infection_forward$age_infection[!is.na(infection_forward$age_infection) & infection_forward$age_infection > infection_forward$max]  <- NA
  infection_forward[is.na(infection_forward$age_infection)
                    ,c("strain",
                       "gam_dens_s_low_treated","gam_dens_s_high_treated","infectivity_tot_treated","gam_dens_s_low","gam_dens_s_high","infectivity_tot","infectivity_tot_noTmt","COI","tmt_status","age_treatment")] <- NA
  
  # add infection characterisitc for iteration it for each host
  # age infection
  out_df[it,-1,"age_inf"] <- infection_forward$age_infection 
  # strain
  out_df[it,-1,"strain"] <-  infection_forward$strain
  # COI
  out_df[it,-1,"COI"] <-  infection_forward$COI
  # treatment status
  out_df[it,-1,"treat"] <-  infection_forward$tmt_status
  # age treatment 
  out_df[it,-1,"age_treatment"] <-  infection_forward$age_treatment
  # infectiousness
  out_df[it,-1,"infectivity_tot"]<- infection_forward$infectivity_tot
  out_df[it,-1,"infectivity_tot_treated"]<- infection_forward$infectivity_tot_treated
  out_df[it,-1,"infectivity_tot_noTmt"]<- infection_forward$infectivity_tot_noTmt
  # gametocyte density
  out_df[it,-1,"gam_dens_s_low_treated"] <- infection_forward$gam_dens_s_low_treated
  out_df[it,-1,"gam_dens_s_high_treated"] <- infection_forward$gam_dens_s_high_treated
  out_df[it,-1,"gam_dens_s_low"] <- infection_forward$gam_dens_s_low
  out_df[it,-1,"gam_dens_s_high"] <- infection_forward$gam_dens_s_high
  # mosquito per host
  out_df[it,-1,"n_mosq"] <- infection_forward$n_mosq
  
  # keep only hosts and not empty cells
  infection_forward <- infection_forward[!is.na(infection_forward$Host_ID),]
  
  # SELECT INFECTIOUS HOSTS
  # total infectious potential across hosts in an untreated population
  n_infectious <- as.data.frame(out_df[it,-1,c("infectivity_tot_noTmt","treat","n_mosq","Host_ID")]) %>% 
    filter(!is.na(treat) & !is.na(infectivity_tot_noTmt) & infectivity_tot_noTmt > 0
           & !is.na(n_mosq) & n_mosq > 0) # here also count for when treat == 1
  n_infectious$infectious_potential <- apply(n_infectious, 1, function(x) as.numeric(x["infectivity_tot_noTmt"])*as.numeric(x["n_mosq"]))
  infectious_potential <- sum(n_infectious$infectious_potential,na.rm=T)
  
  # select hosts that are infectious
  n_transmit <- as.data.frame(out_df[it,-1,c("infectivity_tot","treat","n_mosq","Host_ID")]) %>% 
    filter(!is.na(treat) & !is.na(infectivity_tot) & infectivity_tot > 0
           & !is.na(n_mosq) & n_mosq > 0) # here also count for when treat == 1
  who_infectious <- n_transmit$Host_ID[ n_transmit$treat==0]
  
  # TRANSMISSION HUMAN - MOSQUITO
  # new infections for each of infectious hosts
  for(h in who_infectious){ 
    
    # host's infectiousness of host
    infect <-  as.numeric(infection_forward$infectivity_tot[infection_forward$Host_ID==h])
    
    # host's gametocyte density of each genotype
    G1 <- as.numeric(infection_forward$gam_dens_s_low[infection_forward$Host_ID==h])
    G2 <- as.numeric(infection_forward$gam_dens_s_high[infection_forward$Host_ID==h])
    G1[is.na(G1)] <- 0
    G2[is.na(G2)] <- 0
    
    # host's number of biting mosquitos
    n_mosq <- as.numeric(infection_forward$n_mosq[infection_forward$Host_ID==h])
    
    # calculate the number of mosquitos infected by the host - and specify by which genotype (S_high, S_low or S_mix)
    new_inf  <- fct_transmission_host_to_mosq(new_inf,infect_tot=infect,G1=G1,G2=G2,norm_fact=min(1,(FOI/infectious_potential)),n_mosq = n_mosq) 
    
    rm(list=c("infect","G1","G2"))
  }

# TRANSMISSION MOSQUITO-HUMAN
  # draw random number of new hosts to infect
  new_host <- as.data.frame(out_df[it,2:I+1,c("n_mosq","age_inf","Host_ID")]) %>% filter(is.na(age_inf) & !is.na(Host_ID) & n_mosq>0)
 
  # add hosts if not enough are available
   if(nrow(new_host) < E_ended_inf){
    seq_hosts <- seq(from=ad_host,to= (ad_host+(E_ended_inf-nrow(new_host))-1))
    # update new host numbers
    ad_host <- ad_host+(E_ended_inf-nrow(new_host))
    # add additional hosts to matrix
    out_df[it,is.na(out_df[it,,"Host_ID"]),"Host_ID"][1:(E_ended_inf-nrow(new_host))] <- seq_hosts
    # add number of available mosquitos per new added host
    temp_mosqs <- rnbinom(1:(E_ended_inf-nrow(new_host)), k, 1-(A/(A+k)))
    temp_mosqs[temp_mosqs==0] <- 1 
    out_df[it,out_df[it,,"Host_ID"]%in%seq_hosts,"n_mosq"][-1] <-  temp_mosqs
    
    print(paste("adding : ",E_ended_inf-nrow(new_host), " hosts", sep=""))
    new_host <- as.data.frame(out_df[it,2:I+1,c("n_mosq","age_inf","Host_ID")]) %>% filter(is.na(age_inf) & !is.na(Host_ID))
  }
  
  # assign infected mosquitos ready to infect (timre=0) to available hosts to infect
  new_inf_to_add <- fct_infect_mosq_per_host(new_host = new_host ,infection_df = new_inf[new_inf$timer==0,])
  # total number of infected mosquitos
  new_inf_to_add$n_mosquito <- rowSums(new_inf_to_add[,c("S_high","S_low","S_mix")],na.rm = T)
  # determines if hosts gets a superinfection or a single infection
  new_inf_to_add$super <- ifelse(new_inf_to_add$n_mosquito>1,1,0)
  
  # store the mean number of infectious mosquitos per infected host
  track_mosq$time_step[it] <- it
  track_mosq$average_mosq_host[it] <- mean(new_inf_to_add$n_mosquito)
  # store the fraction of super-infections
  track_mosq$p_super_inf[it] <- sum(new_inf_to_add$super)/nrow(new_inf_to_add)
  
  # if we added new hosts in our dataframe, we need to update
  H <- length(out_df[it,!is.na(out_df[it,,"Host_ID"]),"Host_ID"])-1
  
  # number of hosts to infect
  n_new_inf <- nrow(new_inf_to_add)[1]
  if(n_new_inf != 0){
    # add age treatment
    new_inf_to_add <- left_join(new_inf_to_add,n_s[,c("strain","COI","age_treatement")],by = c("strain","COI"))
    new_inf_to_add$age_treatment <- new_inf_to_add$age_treatement
    new_inf_to_add <- new_inf_to_add[,c("strain","COI","age_treatment","host")]
    
    fill_out <- !is.na(out_df[it,,"Host_ID"])&out_df[it,,"Host_ID"]%in%new_inf_to_add[,"host"]
    out_df[it,fill_out,"age_inf"] <- rep(1,n_new_inf)
    out_df[it,fill_out,"treat"] <- rep("0",n_new_inf)
    which_rows <- as.numeric(names(fill_out[fill_out==T])) + 1 
    
    for(i in which_rows){
      h <- out_df[it,i,"Host_ID"]
      out_df[it,i,"strain"] <- new_inf_to_add[new_inf_to_add$host==h,"strain"]
      out_df[it,i,"COI"] <- new_inf_to_add[new_inf_to_add$host==h,"COI"]
      out_df[it,i,"age_treatment"] <- new_inf_to_add[new_inf_to_add$host==h,"age_treatment"]
    }
    # add infection charactersics of the new infections to the main matrix
    out_df[it,fill_out,"infectivity_tot"] <- sapply(which_rows,match_infectivity_function) 
    out_df[it,fill_out,"infectivity_tot_noTmt"] <- out_df[it,fill_out,"infectivity_tot"]
    out_df[it,fill_out,"infectivity_tot_treated"] <- sapply(which_rows,match_infectivity_treated_function) 
    
    out_df[it,fill_out,"gam_dens_s_low"] <- sapply(which_rows,match_low_gam_function) 
    out_df[it,fill_out,"gam_dens_s_high"] <- sapply(which_rows,match_high_gam_function) 
    
    out_df[it,fill_out,"gam_dens_s_low_treated"] <- sapply(which_rows,match_low_gam_treated_function) 
    out_df[it,fill_out,"gam_dens_s_high_treated"] <- sapply(which_rows,match_high_gam_treated_function) 
    
    # free the space in data frame
    new_inf$strain[new_inf$timer==0] <- -1000 
    new_inf$timer[new_inf$timer==0] <- -1000 
  }
  # remove one time step to the timer of the new infections
  new_inf$timer <- new_inf$timer-1
  rm(list= c("infectious_potential","who_infectious","n_infectious","n_ended_infections","infection_forward","prev_age_inf"))
}

# ---------------------------------------  FINALISATION - format output before saving ---------------------------------------------------------------------------------------------------------------------

rm(new_inf)

# main matrix - with simulation specification and all infections characteristics for each host at each time step
out_n <- as.data.frame.table(out_df)
out_n <- out_n[out_n$Var2 != "time_step",]
colnames(out_n) <- c("time_step","host","measurement","output")
out_n$seed <- seed
out_n$p_superInf <- p_superInf
out_n$cross_react <- cross.react
out_n$tmt_prob <- P_tmt

# all we need of is infected versus not infected 
out_n %>% mutate_if(is.factor, as.character) -> out_n

subset_out <- out_n[out_n$measurement%in%c("age_inf","strain","COI","treat"),]

out_wide <- dcast(subset_out, ... ~ measurement, value.var="output")
out_wide <- out_wide[!is.na(out_wide$treat),]
out_wide <- out_wide[out_wide$treat!=1,]
out_wide$infected <- 0
out_wide$age_inf <- as.numeric(as.character(out_wide$age_inf))
out_wide$age_inf[is.na(out_wide$age_inf)] <- -100
out_wide <- out_wide[out_wide$age_inf!=100,] 

out_wide$infected[out_wide$age_inf>0] <- 1

# summed the number of infected hosts per genotype(s) at each time step
n_inf <- aggregate(infected ~ time_step + strain + seed + p_superInf + cross_react + tmt_prob + COI, data = out_wide, FUN = sum)
n_inf$time_step <- as.numeric(as.character(n_inf$time_step))

