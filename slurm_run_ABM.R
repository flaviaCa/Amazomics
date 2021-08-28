# ---------------------------------------  Slurm_run_ABM.R ---------------------------------------------------------------------------------------------------------------------
# Run the transmission model
# Script defining input parameters, calls the main script and stores the output
#
# created by fcamponovo@hsph.harvard.edu

# ---------------------------------------  Libraries ---------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(reshape2)

# ---------------------------------------  Define paths  ---------------------------------------------------------------------------------------------------------------------

main_dir <- "/n/home09/fcamponovo/Conversion-selection/"
main_dir_code <- paste(main_dir,"code/",sep="")
main_dir_analyis <- paste(main_dir,"Analysis/",sep="")

# path to the individual within host dynamics:
within_host_df <-  paste(main_dir_analyis,"Within-host-dynamics/out_within-hosts.Rdata",sep="")

# baseline fraction super-infection for different average infectious mosquitos per infected host.
# This table was computed running the simulation assuming only one genotype in the population
base_co.inf <- read.csv(file=paste(main_dir,"/Analysis/Inputs/baseline_one_genotype.csv",sep=""))

# define where to store output
main_output_dir <- paste(main_dir_analyis,"Transmission_model/",sep="")
secondary_folder <- paste(main_output_dir,"/aggregated",sep="")

# create the output folder if non-existent
if (!(file.exists(main_output_dir))){
  dir.create(file.path(main_output_dir))
} 
if (!(file.exists(secondary_folder))){
  dir.create(file.path(secondary_folder))
}
# ---------------------------------------  Load within host dynamics  ---------------------------------------------------------------------------------------------------------------------

# load the within host dynamics data frame - include infectivity, aseuxal and gametocyte densities, for treated and untreated infections
# inludes dynamics for different genotypes
source(file=paste(within_host_df))
# remove unnecessary df
rm("long_out")

# ---------------------------------------  Load functions  ---------------------------------------------------------------------------------------------------------------------

source(file=paste(main_dir_code,"/Functions/fct_assign_wih_mst.R",sep=""))
source(file=paste(main_dir_code,"/Functions/fct_infect_mosq_per_host.R",sep=""))
source(file=paste(main_dir_code,"/Functions/fct_transmission_host_to_mosq.R",sep=""))
source(file=paste(main_dir_code,"/Functions/fct_select_strain.R",sep=""))

# ---------------------------------------  Define Parameters  ---------------------------------------------------------------------------------------------------------------------

# simulation parameters
cmd_args <- commandArgs(trailingOnly = TRUE) #  : parameter input from submission file
Time <- 3000  #   : time of simulation, in a 1-day time step
seed <- cmd_args[4]

# parasite parameters:
conv_low <- 0.05      #   : conversion rate for S_low
conv_high <- 0.2      #   : conversion rate for S_high
tau1 <- 9             #   : gametocyte development time S_low
tau2 <- 9             #   : gametocyte development time S_high
cross.react <- 0.75   #   : assumed immune cross reactivity level between the two genotypes
asex.mult <- 16       #   : select the multiplication factor of asexual replication
n_spz <- 10           #   : specify how many sporozoites are infecting humans

# vector parameter
k <- 2.5           #   : dispersion parameter of the negative binomial distribution of the number of mosquitos per host
A <- 10            #   : mean number of mosquitos per host

# population parameters:  
I <- 5500          #   : number of total hosts (this allows for some room if number of infections increase slightly during simulation)
H <- 5000          #   : number of initial infected hosts 
prop_low_strain <- cmd_args[2]  # : initial proportion of S_low genotype infections among the monoclonal infections
P_tmt <- as.numeric(cmd_args[3])  # : Fraction of treated infections

expected_superInf <-  as.numeric(cmd_args[1])  # : expected fraction super-infections
# initial fraction polyclonal infections, assuming a given super-infection level
p_superInf <- base_co.inf$Baseline.coInf[base_co.inf$p_superInf == expected_superInf] # observed fraction super-infections in baseline simulations with only one genotype

# define the frequency of monoclonal and polyclonal infections
fq_strain <- data.frame("strain"= c(rep(c("S_low","S_high","S_mix"),times=2))
                        ,"COI" = c(rep("clonal",times=3),rep("coinfect",times=3))
                        ,"fq_strain" = c((1-p_superInf)*prop_low_strain,(1-p_superInf)*(1-prop_low_strain)
                            ,0,0,0
                            ,p_superInf))

# ---------------------------------------  Run simulations  ---------------------------------------------------------------------------------------------------------------------

source(paste(main_dir_code,"/main.R",sep=""))

# ---------------------------------------  Store output  ---------------------------------------------------------------------------------------------------------------------

#out_n %>% mutate_if(is.factor, as.character) -> out_n
#save(out_n,strain_df,n_s,treated_cleared, file = paste(main_output_dir,"/out_ABM_CoInf_",p_superInf,"pLow",prop_low_strain,"_Tmt_",P_tmt,"seed",seed,"_",Sys.Date(),".Rdata",sep=""))
save(n_inf,strain_df,n_s,treated_cleared,FOI,track_mosq, file = paste(secondary_folder,"/aggregated_out_ABM_CoInf_",p_superInf,"pLow",prop_low_strain,"_Tmt_",P_tmt,"seed",seed,"_",Sys.Date(),".Rdata",sep=""))
