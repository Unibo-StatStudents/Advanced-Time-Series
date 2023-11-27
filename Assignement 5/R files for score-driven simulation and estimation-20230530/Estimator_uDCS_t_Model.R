################
##### Estimation
################
uDCS_t_model_estimator <- function(dati, param){
  
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
  optimizer <- suppressWarnings(nlminb(start = theta_st, objective = interprete_uDCS_t_model, 
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
  
################
##### Interprete
################
interprete_uDCS_t_model <- function(dati, param){
    
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
    fitness <- uDCS_t_model_filter(dati, theta_new)$Log_Likelihood
    
    if(is.na(fitness) | !is.finite(fitness)) fitness <- -1e10
    if(fitness != fitness) fitness <- -1e10
    
  return(-fitness)
} 

