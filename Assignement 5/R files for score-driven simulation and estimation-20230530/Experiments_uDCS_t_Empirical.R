##################################################################################
# this file does fit a DCS-t location model to data that have to be uploaded 
# it uses the functions defined in Estimator, Filter
##################################################################################
# initial values of the parameters that will enter into theta (static parameters)
# better call it theta_0 
omega <- 0
phi <- 0.8
k <- 0.8
varsigma <- 1
nu <- 8

theta <- c(omega, phi, k, varsigma, nu)

# estimate the parameters (uses functions in Estimator)
# input: data, y, and INITIAL parameter values, theta
est  <- uDCS_t_model_estimator(y, theta)
est$theta_list

# recover the estimated time varying signal (uses functions in Filter)
# input: data, y, and ESTIMATED parameter values, theta
filter <- uDCS_t_model_filter(y, est$theta)
ts.plot(y)
lines(filter$Dynamic_Location,col = "red")

# compare the estimated one-step-ahead prediction error with the score 
# both are mds, the score is nonlinear and robust if nu is small, 
# if nu->infinity, then the score u_t converges to the one-step-ahead prediction error
# which is then linear and Gaussian

T = NROW(y)
v = 0 
for(t in 2:(T+1)){
  v[t-1]=y[t-1]-filter$Dynamic_Location[t-1]
}

# plot of the innnovation error v_t against the robust score u_t
ts.plot(v)
lines(filter$Innovation_u_t,col="green")
