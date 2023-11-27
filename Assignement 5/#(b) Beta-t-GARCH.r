    #(b) Beta-t-GARCH


    Beta_martingale_diff_u_t <- function(y, varsigma, nu){
    
    u_t <- c(((nu +1)*y^2)/((nu-2)*varsigma +y^2)-1)
    
    return(u_t)
    }

    uSTDT_rnd <- function(n, mu, varsigma = 1, nu = 3) {
    
    z <- rt(n, df = nu) 
    y <- numeric()
    for(i in 1:n){
        y[i] <- c(mu + z[i]* sqrt(varsigma) ) 
    }
    
    return(y)
    }



    Beta_model_simulator <- function(T, omega, phi, k, varsigma, nu){
    
    ###Define the Processes
    y   <- array(data = NA, dim = c(T) ) 
    
    ###Define Dynamic Scale and Innovations
    sigma_t <- array(data = NA, dim = c(T))
    u_t  <- array(data = NA, dim = c(T-1))
    
    ###Initial value for the recursion 
    sigma_t[1]  <- omega
    
    ###Generate the first observations of the process
    y[1]   <- uSTDT_rnd(1, varsigma)
    
        ###Dynamics 
        for (t in 2:T) {
        
        ###Factor Innovations
        u_t[t-1] <- Beta_martingale_diff_u_t(y[t-1], sigma_t[t-1], nu)
            
            ###Updating Filters                    
            sigma_t[t]  <- omega + phi * sigma_t[t-1]^2 + k*sigma_t[t-1]^2*u_t[t-1]
            
            ###Generate the observations of the processes
            y[t] <- uSTDT_rnd(1, varsigma)
        }
    ######################
    ####### OUTPUT #######
    ######################
    
    ###Make List
    out <- list(y_t_gen          = as.ts(y),
                Dynamic_Scale    = as.ts(sigma_t),
                Innovation_u_t   = as.ts(u_t))
    
    return(out)
    }


a <- Beta_model_simulator(T = 20, omega = 1, phi = 1, k = 1, varsigma = 1, nu = 10)
a





Beta_t_model_filter <- function(y, theta){
  
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
  sigma_t <- array(data = NA, dim = c(T+1))
  u_t  <- array(data = NA, dim = c(T))
  
  ###Initialize Dynamic Location
  sigma_t[1]   <- (omega)
    
    ###Initialize Likelihood
    dloglik[1] <- uSTDT_t(y[1], sigma_t[1], varsigma = varsigma, nu = nu, log = TRUE)
    loglik     <- dloglik[1]
    
    for(t in 2:(T+1)) {
      ###Dynamic Location Innovations
      u_t[t-1] <- Beta_martingale_diff_u_t(y[t-1], sigma_t[t-1], varsigma)
      ###Updating Filter                    
      sigma_t[t]   <- omega + phi * (sigma_t[t-1] - omega) + k * u_t[t-1]
    
      if(t < (T+1)){
        ###Updating Likelihoods
        dloglik[t] <- uSTDT_t(y[t], sigma_t = sigma_t[t], varsigma = varsigma, nu = nu, log = TRUE)
        loglik     <- loglik + dloglik[t]
      }
    }

  ######################
  ####### OUTPUT #######
  ######################
  sigma_t <- ts(sigma_t, start = start(y), frequency = frequency(y))
  u_t  <- ts(u_t, start = start(y), frequency = frequency(y))
  
  ###Make List
  out <- list(Dynamic_Location = sigma_t,
              Innovation_u_t   = u_t,
              Log_Densities_i  = dloglik,
              Log_Likelihood   = loglik)
  
  return(out)
}

Beta_t_model_filter(y=c(1,2,3,4,5,6,4,3,2,1), theta = c(0.8,0.8,1,0.7,10))

uSTDT_t <- function(y, sigma_t, varsigma, nu, log = TRUE){
  
  ulpdf <- (lgamma((nu + 1) / 2) - lgamma(nu / 2) - (1/2) * log(varsigma) -
              (1/2)  * log(pi * nu) - ((nu + 1) / 2) * log(1 + (y - sigma_t)^2 / (nu*varsigma) ))
  
  if(log != TRUE){
    ulpdf <- exp(ulpdf)
  } 
  
  return(ulpdf)
}

Beta_t_model_estimator <- function(dati, param){
  
  Start <- Sys.time()
  ###Take T
  T <- length(dati)
  
  ###Parameter Selections Dynamic Location
  omega <- param[1]
  phi   <- param[2]
  k     <- param[3]
  
  varsigma <- param[4]
  nu       <- param[5]

  ###Create a vector with the parameters
  theta_st <- c(omega, phi, k, varsigma, nu)
  
  ###Take Bounds
  lower <- c(-Inf, -0.999, -2, 1e-05, 2.099)
  upper <- c( Inf,  0.999,  2, Inf, 300)
  
  #------> Optimize uDCS_t_model Filters w/L-BFGS-B
  #optimizer <- suppressWarnings(optim(par = theta_st, fn = interprete_uDCS_t_model, 
  #                                    dati = dati, method = "L-BFGS-B", 
  #                                    control = list(trace = 1), hessian = FALSE,
  #                                    lower = lower, upper = upper))
  
  #------> Optimize every uDCS_t_model Filters w/solnp := WAY faster than optim
  #optimizer <- suppressWarnings(Rsolnp::solnp(pars = theta_st, fun = interprete_uDCS_t_model, 
  #                                            dati  = dati, control = list(trace = 1), 
  #                                            LB = lower, UB = upper)) 
  
  #------> Optimize every uDCS_t_model Filters w/nlminb := WAY faster than optim
  optimizer <- suppressWarnings(nlminb(start = theta_st, objective = interprete_Beta_t_model, 
                                       dati  = dati, gradient = NULL, 
                                       control = list(trace = 0), hessian = NULL,
                                       lower = lower, upper = upper))
  
  #------> Save the optimized parameters Dynamic Location
  omega_opt <- optimizer$par[1]  
  phi_opt <- optimizer$par[2]
  k_opt  <- (optimizer$par[3])
  
  varsigma_opt <- optimizer$par[4]
  nu_opt       <- optimizer$par[5]
  
  ###Create a vector with ALL the optimized parameters
  theta_opt <- c(omega_opt, phi_opt, k_opt, varsigma_opt, nu_opt)
  
  ###Create a list with ALL the optimized parameters
  theta_list <- list(omega = omega_opt,
                     phi = phi_opt,
                     k  = k_opt,
                     varsigma = varsigma_opt,
                     nu    = nu_opt)
  
  ######################
  ####### OUTPUT #######
  ######################
  
  #------> Some detail
  Elapsed_Time <- Sys.time() - Start
  print(paste("Elapsed Time: ", toString(Elapsed_Time)))
  
  ###Make List
  out <- list(theta_list = theta_list,
              theta      = theta_opt,
              optimizer  = optimizer)
  
  return(out) 
}

interprete_Beta_t_model <- function(dati, param){
    
  ###Take T
  T <- length(dati)
  
  ###Parameter Selections Dynamic Location
  omega <- param[1]
  phi   <- param[2]
  k     <- param[3]
  
  varsigma <- param[4]
  nu       <- param[5]
    
   ###Create a new vector with the parameters
   theta_new <- c(omega, phi, k, varsigma, nu)
    
    #------> Fitness Functions
    fitness <- Beta_t_model_filter(dati, theta_new)$Log_Likelihood
    
    if(is.na(fitness) | !is.finite(fitness)) fitness <- -1e10
    if(fitness != fitness) fitness <- -1e10
    
  return(-fitness)
} 



set.seed(1234)
T <- 200
omega <- 0
phi <- 0.8
k <- 0.8
varsigma <- 1

nu_values <- c(3, 10, 200)

for (nu in nu_values){

    theta <- c(omega, phi, k, varsigma, nu)

    simu <- Beta_model_simulator(T, omega, phi, k, varsigma, nu)
    y = simu$y_t_gen
    ts.plot(y)
    
    est  <- Beta_t_model_estimator(y, theta)
    est$theta_list

    filter <- Beta_t_model_filter(y, est$theta)
    
    plot_title <- paste("nu = ", nu)
    plot(y, main = plot_title)
    lines(filter$Dynamic_Location,col = "red")
}



