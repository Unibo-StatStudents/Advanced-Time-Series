################
##### Simulation
################
uDCS_t_model_simulator <- function(T, omega, phi, k, varsigma, nu){
  
  ###Define the Processes
  y   <- array(data = NA, dim = c(T) ) 
  
  ###Define Dynamic Location and Innovations
  mu_t <- array(data = NA, dim = c(T))
  u_t  <- array(data = NA, dim = c(T-1))
  
  ###Initial value for the recursion 
  mu_t[1]   <- omega
  
  ###Generate the first observations of the process
  y[1]   <- uSTDT_rnd(1, mu_t[1], varsigma, nu)
  
    ###Dynamics 
    for (t in 2:T) {
      
      ###Factor Innovations
      u_t[t-1] <- martingale_diff_u_t(y[t-1], mu_t[t-1], varsigma, nu)
        
        ###Updating Filters                    
        mu_t[t]   <- omega + phi * (mu_t[t-1] - omega) + k * u_t[t-1]
        
          ###Generate the observations of the processes
          y[t] <- uSTDT_rnd(1, mu_t[t], varsigma, nu)
    }
  ######################
  ####### OUTPUT #######
  ######################
  
  ###Make List
  out <- list(y_t_gen          = as.ts(y),
              Dynamic_Location = as.ts(mu_t),
              Innovation_u_t   = as.ts(u_t))
  
  return(out)
}

############################################################################
################# Univariate Student's t Random Generator ##################
############################################################################

uSTDT_rnd <- function(n, mu, varsigma, nu) {
  
  z <- rt(n, df = nu) 
  y <- numeric()
  for(i in 1:n){
    y[i] <- c(mu + z[i]* sqrt(varsigma) ) 
  }
  
  return(y)
}
