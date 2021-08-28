# With host model as in McKenzie et al,  with added treatment effect

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

# ------ LOAD FUNCTIONS ------------------------------------------------------------------------------------
source(file=paste("~/Functions/fct_wih_DE.R",sep=""))

# ------ SET PARAMETERS ------------------------------------------------------------------------------------

time <- seq(from=0, to=200, by = 1) # :   infection age

# define constant parameters
a_u <- log(16)/2  #    :   daily replication rate of asexual blood forms
y_u <- 1          #    :   survivorship over sequestered period
c_u <- 4.7e-4     #    :   daily removal rate of asexual forms by innate-immunity effectors (asexual forms per effector)
b_uv_fullCR <- matrix( rep(2.3e-4,4), nrow = 2, ncol = 2)
#    :   daily removal rate of asexual forms by acquired-immunity effectors (asexual forms per effector) . diag = 0 if no cross reactivity
s_v <- 0.168      #    :   rate of stimulation of the innate immune response (effectors per asexual form)
r_uv_fullCR <- matrix( rep(0.067,4), nrow = 2, ncol = 2)
#    :   rate of stimulation of the acquired immune response (effectors per asexual form) diag = 0 if no cross reactivity
k_v <- 2.8e-5     #    :   removal rate of innate-immunity effectors (effectors per effector per asexual form)
w_u <- 0          #    :   removal rate of acquired-immunity effectors (effectors per effector per asexual form)
q <- 0.285        #    :   decay rate of innate-immunity effectors
p <- 0.50         #    :   decay rate of mature gametocytes

fact_G <- 0       #   :   effect of treatment on gametocytes
P_f <- 10^3.5     #   :   threshold of fever, and thus of treatment - asexual parasite density - random value!

# Define varying parameters:
factor_CR <- seq(0,1,by=.25)
#    :   Cross reactivity, i.e how much does the immune response of one strain affect the other strain (factor between 0 and 1)
factor_activation <- 0.5*factor_CR
#    :   Cross reactivity of immune stimulation, i.e how much does one strain stimulate the immune response to the other strain  (factor between 0 and 1)
tau_all <- 9 #c(6,9) #    :   development time of gametocytes, ie number of days seuqestered before released in the blood stream
z_u_all <- seq(0.01,0.3,by=.01)
#    :   daily conversion rate of asexual forms to gametocytes
fact_G_all <- 0 #c(0,0.25,0.5,1)
#    :   effect of treatment on gametocytes (between 0 and 1)

# initial conditions:
state <- c(M = c(0.01,0.01), G = c(0,0) , J = c(0,0) , I = 0 ) # initial for co-infections
state_one_inf <- c(M = c(0.01,0), G = c(0,0) , J = c(0,0) , I = 0 ) # initial for single infections

# Variables
# M : merozoites
# G : gametocytes
# J : acquired immune response
# I : innate immune response

# ------ RUN MODEL WITHOUT TREATMENT ------------------------------------------------------------------------------------

sim <- 1
switch <- 0
t_treat <- 250 # which day is the day the infecrtion reaches the given threshold, so here no treatment

# create empty dataframe to store outputs
out <- data.frame(sim = integer(), tau_1 = numeric(), tau_2 = numeric(), z_1 = numeric(), z_2 = numeric(),CR.factor = numeric(),  timing_second_inf = numeric(),
                  time = integer(), M1= numeric(), M2= numeric(), G1= numeric(), G2= numeric(),J1 = numeric(),J2= numeric())


# Generate outputs for single infection

b_uv <- b_uv_fullCR # only one infection, thus cross-reactivity doesn't matter
r_uv <- r_uv_fullCR # only one infection, thus cross-reactivity doesn't matter

for(lag in 1:length(tau_all)){ # gametocyte developement time
  for(z in 1:length(z_u_all)){ # sexual conversion rate
    tau1 <- tau_all[lag]
    tau2 <- tau_all[lag]
    z_u <- z_u_all[z]
    time_intro <- c(0,tail(time,1)+1) # second infection doesn't occur
    
    # run the diff. equation solver
    yout <- dede(y = state_one_inf, times = time, func = wih_tau_kill, parms = NULL)
    
    # store the outputs, including parameter values
    temp <- as.data.frame(yout)
    temp$sim <- sim
    temp$tau_1 <- tau1
    temp$tau_2 <- 0
    temp$timing_second_inf <- "none" 
    temp$z_1 <- z_u
    temp$z_2 <- 0
    temp$CR.factor <- "none"
    out <- rbind(out,temp)
    sim <- sim+1
  } # sexual converstion rate 
} # gametocyte developement time


# Generate outputs for co-infection

for(cr in 1:length(factor_CR)){ # cross reactivivity
  for(lag1 in 1:length(tau_all)){ # gametocyte development time for strain 1
    for(lag2 in 1:length(tau_all)){ # gametocyte development time for strain 2
      for(z1 in 1:length(z_u_all)){ # sexual conversion rate of strain 1
        for(z2 in 1:length(z_u_all)){ # sexual conversion rate of strain 2
          b_uv <- b_uv_fullCR * c(1,factor_CR[cr],factor_CR[cr],1) # decrease Cross reactive immunity effect by given factor
          r_uv <- r_uv_fullCR * c(1,factor_activation[cr],factor_activation[cr],1) # decrease Cross reactive immunity activation by given factor
          z_u <- c(z_u_all[z1],z_u_all[z2])
          tau1 <- tau_all[lag1]
          tau2 <- tau_all[lag2]
          time_intro <- c(0,0) # introduction time of both strains
          # run the diff. equation solver
          yout <- dede(y = state, times = time, func = wih_tau_kill, parms = NULL)
          
          # store the outputs, including parameter values
          temp <- as.data.frame(yout)
          temp$sim <- sim
          temp$tau_1 <- tau1
          temp$tau_2 <- tau2
          temp$timing_second_inf <- time_intro[2] 
          temp$z_1 <- z_u[1]
          temp$z_2 <- z_u[2]
          temp$timing_second_inf <- time_intro[2] 
          temp$CR.factor <- factor_CR[cr]
          out <- rbind(out,temp)
          sim <- sim+1
        }# sexual conversion rate of strain 2
      }# sexual conversion rate of strain 1
    } # gametocyte development time for strain 2
  }# gametocyte development time for strain 1
} # cross reactivivity

out$treat <-0
out$fact_G_treat <- "no.tmt"

# ------ RUN MODEL WITH TREATMENT ------------------------------------------------------------------------------------

# create empty dataframe to store outputs
out_treat <- data.frame(sim = integer(), tau_1 = numeric(), tau_2 = numeric(), z_1 = numeric(), z_2 = numeric(),CR.factor = numeric(),  timing_second_inf = numeric(),
                        time = integer(), M1= numeric(), M2= numeric(), G1= numeric(), G2= numeric(),J1 = numeric(),J2= numeric(),age_treat = integer(),fact_G_treat = numeric())


# first check when threshold for treatment is reached
treat_time <- out %>% group_by(sim,time)  %>% mutate(M_tot = sum(M1,M2,na.rm=T)) %>% filter(M_tot>=P_f)%>% group_by(sim) %>% filter(time == min(time))%>% mutate(time_treat = time) %>% select(time_treat,sim,tau_1,tau_2,z_1,z_2,CR.factor) %>% ungroup()

# remove dupplicate
treat_time <- treat_time %>% select(-sim) %>%  distinct()

# Generate outputs for single infection

b_uv <- b_uv_fullCR # only one infection, thus cross-reactivity doesn't matter
r_uv <- r_uv_fullCR # only one infection, thus cross-reactivity doesn't matter

for(lag in 1:length(tau_all)){ # gametocyte developement time
  for(z in 1:length(z_u_all)){ # sexual conversion rate 
    for(fact_G in fact_G_all){ # effect of treatment on gametocytes
      tau1 <- tau_all[lag]
      tau2 <- tau_all[lag]
      z_u <- z_u_all[z]
      time_intro <- c(0,tail(time,1)+1) # second infection doesn't occur
      state_one_inf <- c(M = c(0.01,0), G = c(0,0) , J = c(0,0) , I = 0 )
      # when is the treatment threshold reached?
      t_treat <- as.numeric(treat_time %>% filter(tau_1==tau1 & tau_2==0 & z_1==z_u & z_2==0) %>% select(time_treat))
      # run the diff. equation solver
      yout <- dede(y = state_one_inf, times = time, func = fct_wih_DE, parms = NULL)
      # store the outputs, including parameter values
      temp <- as.data.frame(yout)
      temp$sim <- sim
      temp$tau_1 <- tau1
      temp$tau_2 <- 0
      temp$timing_second_inf <- "none" 
      temp$z_1 <- z_u
      temp$z_2 <- 0
      temp$CR.factor <- "none"
      temp$age_treat <- t_treat
      temp$fact_G_treat <- fact_G
      out_treat <- rbind(out_treat,temp)
      sim <- sim+1
    }# effect of treatment on gametocytes
  } # sexual conversion rate 
} # gametocyte development time

# Generate outputs for co-infection

for(cr in 1:length(factor_CR)){ # cross reactivivity
  for(lag1 in 1:length(tau_all)){ # gametocyte development time for strain 1
    for(lag2 in 1:length(tau_all)){ # gametocyte development time for strain 2
      for(z1 in 1:length(z_u_all)){ # sexual conversion rate of strain 1
        for(z2 in 1:length(z_u_all)){ # sexual conversion rate of strain 2
          for(fact_G in fact_G_all){ # effect of treatment on gametocytes
            b_uv <- b_uv_fullCR * c(1,factor_CR[cr],factor_CR[cr],1) # decrease Cross reactive immunity effect by given factor
            r_uv <- r_uv_fullCR * c(1,factor_activation[cr],factor_activation[cr],1) # decrease Cross reactive immunity activation by given factor
            z_u <- c(z_u_all[z1],z_u_all[z2])
            tau1 <- tau_all[lag1]
            tau2 <- tau_all[lag2]
            time_intro <- c(0,0) # time of introduction of each strain
            # when is the treatment threshold reahced?
            t_treat <- as.numeric(treat_time %>% filter(tau_1==tau1 & tau_2==tau2 & z_1==z_u[1] & z_2==z_u[2] & CR.factor == factor_CR[cr]) %>% select(time_treat))
            # run the diff. equation solver
            yout <- dede(y = state, times = time, func = wih_tau_kill, parms = NULL)
            # store the outputs, including parameter values
            temp <- as.data.frame(yout)
            temp$sim <- sim
            temp$tau_1 <- tau1
            temp$tau_2 <- tau2
            temp$timing_second_inf <- time_intro[2] 
            temp$z_1 <- z_u[1]
            temp$z_2 <- z_u[2]
            temp$timing_second_inf <- time_intro[2] 
            temp$age_treat <- t_treat
            temp$fact_G_treat <- fact_G
            temp$CR.factor <- factor_CR[cr]
            out_treat <- rbind(out_treat,temp)
            sim <- sim+1
            
          } # effect of treatment on gametocytes
        } # sexual conversion rate of strain 2
      } # sexual conversion rate of strain 1
    } # gametocyte development time for strain 2
  } # gametocyte development time for strain 1
} # cross reactivivity

out_treat$treat <- 1

# ------ MERGE DATA FRAMES AND ADD INFECTIVITY ------------------------------------------------------------------------------------

# output without treatment:

# put a threshold on asexual concentration:
out$M1[out$M1<10^-5] <- 0
out$M2[out$M2<10^-5] <- 0
# add total asexual concentration:
out$M_tot <- out$M1 + out$M2
out$G_tot <- out$G1 + out$G2

# add infectivity
out$infectivity_1 <- NA
out$infectivity_2 <- NA
out$infectivity_tot <- NA
out$infectivity_1[out$G1>2/3] <- 1.08*exp(-exp(-0.86*(log10(out$G1[out$G1>2/3] )-1.48))) 
out$infectivity_2[out$G2>2/3] <- 1.08*exp(-exp(-0.86*(log10(out$G2[out$G2>2/3] )-1.48)))
out$infectivity_tot[out$G_tot>2/3] <- 1.08*exp(-exp(-0.86*(log10(out$G_tot[out$G_tot>2/3] )-1.48)))

# output with treatment

# put a threshold on asexual concentration:
out_treat$M1[out_treat$M1<10^-5] <- 0
out_treat$M2[out_treat$M2<10^-5] <- 0
# add total asexual concentration:
out_treat$M_tot <- out_treat$M1 + out_treat$M2
out_treat$G_tot <- out_treat$G1 + out_treat$G2

# add infectivity
out_treat$infectivity_1 <- NA
out_treat$infectivity_2 <- NA
out_treat$infectivity_tot <- NA
out_treat$infectivity_1[out_treat$G1>2/3] <- 1.08*exp(-exp(-0.86*(log10(out_treat$G1[out_treat$G1>2/3] )-1.48))) 
out_treat$infectivity_2[out_treat$G2>2/3] <- 1.08*exp(-exp(-0.86*(log10(out_treat$G2[out_treat$G2>2/3] )-1.48)))
out_treat$infectivity_tot[out_treat$G_tot>2/3] <- 1.08*exp(-exp(-0.86*(log10(out_treat$G_tot[out_treat$G_tot>2/3] )-1.48)))

# merge both data frames
out <- out  %>% select(-sim) %>%  distinct()
out_all <- merge(out, out_treat,by=c("time","tau_1","tau_2","z_1","z_2","timing_second_inf","CR.factor"))

# long data
long_out <- out  %>% gather(variable, value,c(M1,M2,M_tot,G1,G2, G_tot,infectivity_1,infectivity_2, infectivity_tot,J1,J2,I))
long_out_treat <- out_treat %>% select(-c(sim)) %>% gather(variable, value,c(M1,M2,M_tot,G1,G2, G_tot,infectivity_1,infectivity_2, infectivity_tot,age_treat,J1,J2,I))
long_out <- rbind(long_out,long_out_treat)

# ------ SAVE OUTPUT ------------------------------------------------------------------------------------

save(out_all,long_out,file=paste("~/Within-host-dynamics/out_within-hosts_",Sys.Date(),".Rdata",sep=""))