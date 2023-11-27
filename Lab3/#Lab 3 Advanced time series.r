#Lab 3 Advanced time series

################################################################################
# we shall create our own functions to simulate and estimate the parameters as 
# well as recovering the predicted state by the KF for an AR(1) 
# signal plus noise model with our own created function 
# 
### warning: this file will be quite rough, you will then make it both personal and elegant 
################################################################################

# dgp stands for data generating process, i.e we are simulating an AR(1) plus noise 
# process: the input is made of the length of the series and the static parameters 

dgp <-function(n, phi, sigma_e, sigma_eta){
         "n: length of the time series"
         "phi: coefficient of the autoregressive term"
         "sigma_e: variance of the noise term"
         "sigma_eta: variance of the AR process innovation term"

eta = sigma_eta * rnorm(n)
        
e = sigma_e * rnorm(n)

mu = 0
mu[1] = 0  # unconditional mean of an AR(1) process with no intercept
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}

ts.plot(y)
# lines(mu,col='red')

out <- list(y = as.ts(y), mu = as.ts(mu))

return(out)
        "y: simulated time series"
        "mu: the underlying AR process whitout the noise term"
}


################################################################################
# KALMAN filter recursions with our own function
################################################################################

KF <-function(y, phi, sigma_e, sigma_eta){
    "input: same as dgp function"

# we here separate the arguments of the function that will be estimated - phi, sigma_e,
# sigma_eta by the arguments of the function that will not - mu_1|0 and P_1|0 
  
m10 = 0
P10 = 1
  
n = NROW(y)
  
# allocate space 
mu_pred = array(data = NA, dim = c(n))  # this is mu_{t|t-1}
P = 0                                   # this is P_{t|t-1} 
v = array(data = NA, dim = c(n))        # this  will be the innovation error 
K = 0                                   # the Kalman gain 
F = 0                                   # the conditional variance of v_t 
dllk = array(data = NA, dim = c(n))    #the log-likelihood value 
llk = 0 

# initialise the recursion 

mu_pred[1] = m10
P[1] = P10

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  dllk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
  llk  = llk + dllk[t]
}

#llk = sum(dllk)

out <- list(mu_pred = as.ts(mu_pred), llk = llk)
#out <- list(v, mu_pred)

return(out)
        "mu_pred: the filtered estimates of the hidden state at each time step"
        "llk: the log-likelihood of the observed time series"
}


################################################################################
# KALMAN filter likelihood 
################################################################################

loglikelihood <- function(par,y){
    "par: A vector of model parameters (phi, sigma_e, sigma_eta)"
    "y is the observed time series"

phi       <- par[1]
sigma_e   <- par[2]
sigma_eta <- par[3]

theta_new <- c(phi, sigma_e,sigma_eta)

obj = KF(y, par[1],par[2],par[3])$llk

return(-obj)
        "-obj: the negative log-likelihood value"
}


################################################################################
# KALMAN filter parameter estimation
################################################################################
# input: the data and the initial values of the parameters that we aim to estimate
# objective function: the KF Gaussian likelihood given from the KF 
# (prediction error decomposition), see the function above 
# output: the estimated parameters
# the filtered mu will be recovered from the KF recursions 

estimator <- function(y,par){

n = NROW(y)

phi       <- par[1]
sigma_e   <- par[2]
sigma_eta <- par[3]

theta_0 <- c(phi, sigma_e, sigma_eta)

#hat = nlminb(theta_0, loglikelihood)
hat = nlminb(start = theta_0, objective = loglikelihood, y = y)

hat_phi = hat$par[1]
hat_sigma_e   <- hat$par[2]
hat_sigma_eta <- hat$par[3]

theta_hat <- c(hat_phi,hat_sigma_e,hat_sigma_eta)


###Create a list with ALL the optimized parameters
theta_list <- list(hat_phi = hat_phi,
                   hat_sigma_e = hat_sigma_e,
                   hat_sigma_eta = hat_sigma_eta)

out <- list(theta_list = theta_list)

return(out) 
}          

#optimizer <- suppressWarnings(nlminb(start = theta_0, objective = llk, 
 #                                    dati  = y, gradient = NULL, 
  #                                   control = list(trace = 0), hessian = NULL,
   #                                  lower = lower, upper = upper))



library(TSA)
data(Nile)

y = Nile

sigma_e = 15099
sigma_eta = 1469.1 
phi = 1 

filter <- KF(y, phi, sigma_e, sigma_eta)
est_res <- estimator(y, c(phi, sigma_e, sigma_eta))
plot(y)
lines(filter$mu_pred,col = "green")




################################################################################
# simulate an AR(1) signal plus noise model and estimate the static parametrs 
# and recover the filtered state by the KF 
################################################################################

############################
# set some parameter values 
############################
n = 300
phi = 0.8
sigma_e = 1
sigma_eta = 1
##############

dgp = dgp(n,phi, sigma_e, sigma_eta)

theta = c(phi,  sigma_e, sigma_eta)

y_sim = dgp$y
true_mu = dgp$mu

ts.plot(y_sim)
lines(true_mu, col='red')

y = y_sim
theta_0 = c(0.7,.8, .7)

est = estimator(y,theta_0)
est$theta_list

theta_hat = unlist(est$theta_list)
theta_hat = as.numeric(theta_hat)

filter <- KF(y,theta_hat[1],theta_hat[2],theta_hat[3])
#ts.plot(y)
lines(filter$mu_pred,col = "green")
lines(true_mu, col='blue')






# Load your dataset or create a sample dataset
# For example, let's assume your dataset is in a CSV file named "Nile.csv"
dataset <- read.csv("Nile.csv")

# Extract the time series data from your dataset
# Assuming the time series is in a column named "value"
time_series <- dataset$value

# Specify the initial parameter values
phi_initial <- 0.5
sigma_e_initial <- 1
sigma_eta_initial <- 1

