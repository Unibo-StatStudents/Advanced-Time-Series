###############
##### Filtering
###############
uDCS_t_model_filter <- function(y, theta){
  
  ###Take T
  T <- length(y)
  
  ###Define LogLikelihoods
  dloglik <- array(data = NA, dim = c(T))
  loglik  <- numeric()
  
  ###Parameter Selections Dynamic Location
  omega <- theta[1]
  phi   <- theta[2]
  k     <- theta[3]
  
  varsigma <- theta[4]
  nu       <- theta[5]
  
  ###Define Dynamic Location and Innovations
  mu_t <- array(data = NA, dim = c(T+1))
  u_t  <- array(data = NA, dim = c(T))
  
  ###Initialize Dynamic Location
  mu_t[1]   <- (omega)
    
    ###Initialize Likelihood
    dloglik[1] <- uSTDT_uDCS_t(y[1], mu_t[1], varsigma = varsigma, nu = nu, log = TRUE)
    loglik     <- dloglik[1]
    
    for(t in 2:(T+1)) {
      ###Dynamic Location Innovations
      u_t[t-1] <- martingale_diff_u_t(y[t-1], mu_t[t-1], varsigma, nu)
      ###Updating Filter                    
      mu_t[t]   <- omega + phi * (mu_t[t-1] - omega) + k * u_t[t-1]
    
      if(t < (T+1)){
        ###Updating Likelihoods
        dloglik[t] <- uSTDT_uDCS_t(y[t], mu_t = mu_t[t], varsigma = varsigma, nu = nu, log = TRUE)
        loglik     <- loglik + dloglik[t]
      }
    }

  ######################
  ####### OUTPUT #######
  ######################
  mu_t <- ts(mu_t, start = start(y), frequency = frequency(y))
  u_t  <- ts(u_t, start = start(y), frequency = frequency(y))
  
  ###Make List
  out <- list(Dynamic_Location = mu_t,
              Innovation_u_t   = u_t,
              Log_Densities_i  = dloglik,
              Log_Likelihood   = loglik)
  
  return(out)
}

############################################
####### ADDITIONAL FUNCTIONS #######
############################################

########################################################
martingale_diff_u_t <- function(y, mu_t, varsigma, nu){
  
  u_t <- c((1 / (1 + (y - mu_t)^2/(nu*varsigma)) * (y - mu_t)))
  
  return(u_t)
}

########################################################
uSTDT_uDCS_t <- function(y, mu_t, varsigma, nu, log = TRUE){
  
  ulpdf <- (lgamma((nu + 1) / 2) - lgamma(nu / 2) - (1/2) * log(varsigma) -
              (1/2)  * log(pi * nu) - ((nu + 1) / 2) * log(1 + (y - mu_t)^2 / (nu*varsigma) ))
  
  if(log != TRUE){
    ulpdf <- exp(ulpdf)
  } 
  
  return(ulpdf)
}