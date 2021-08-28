# differential equations of the within host dynamics
fct_wih_DE <- function(t, y, parameters){
  tlag1 <- t - tau1
  tlag2 <- t - tau2
  
  bool <- as.numeric(t >= time_intro) # let's see if we can do something here
  if(tlag1 < 0)
    Mlag1 <- bool[1]*c(state[["M1"]],state[["M2"]]) # takes the initial condition, or leaves it at 0 if time of introduction is not yet reached (bool = 0 for the second infection)
  else
    Mlag1 <- lagvalue(tlag1)
  
  if(tlag2 < 0)
    Mlag2 <- bool[2]*c(state[["M1"]],state[["M2"]]) # takes the initial condition, or leaves it at 0 if time of introduction is not yet reached (bool = 0 for the second infection)
  else
    Mlag2 <- lagvalue(tlag2)
  
  if(switch == 0){
    e_kill <- 0
    if(t >= t_treat){
      e_kill <- 10
      switch <- 1
    }
  }
  else{
    e_kill <- 10
    switch <- 1
  }
  dy1 <- (a_u - z_u - c_u*y[["I"]] - b_uv %*% c(y[["J1"]],y[["J2"]])-e_kill) * bool*c(y[["M1"]], y[["M2"]]) # one equation per genotype
  dy3 <- z_u * y_u * c(Mlag1[1],Mlag2[2]) - p * c(y[["G1"]],y[["G2"]]) -fact_G*e_kill * c(y[["G1"]],y[["G2"]])  # one equation per genotype
  dy5 <- bool * r_uv %*% c(y[["M1"]],y[["M2"]]) - w_u * c(y[["J1"]],y[["J1"]]) # one equation per genotype
  dy7 <- sum(s_v * bool * c(y[["M1"]],y[["M2"]])) - y[["I"]]*sum(k_v *  bool*c(y[["M1"]],y[["M2"]])) - q* y[["I"]] # one for all genotype
  list(c(dy1,dy3,dy5,dy7))
}
