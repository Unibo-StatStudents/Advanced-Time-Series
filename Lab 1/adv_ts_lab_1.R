################################################################################
# Simulate a LLM with varying values of q, arising from the pairs: 
# (sigma_e, sigma_eta) = (0.1,1);(1,1);(1,0.1)
################################################################################

# fix the length of the series (the sample size in real data)

n = 300

# fix the sqrt of the variance of eta (notation of the slides)
# and simulate the error process eta 

sigma_eta = 0.5
eta = sigma_eta * rnorm(n)
mean(eta)
var(eta)

# check: starts
eta
ts.plot(eta)

acf(eta, drop.lag.0 = FALSE, lag.max = 120) #autocorrelation function with the maximum of lags of 120
# check: ends

# simulate epsilon_t 

sigma_e = 1
e = sigma_e * rnorm(n)

# check: starts
e
ts.plot(e)
acf(e, lag.max = 80, drop.lag.0 = FALSE) #autocorrelation function with the maximum of lags of 120
# check: ends

# q, signal-to-noise ratio (exact from simulations)

q = sigma_eta^2/sigma_e^2

# sample signal-to-noise ratio 
# hat_q = var(eta)/var(e)

# start the recursion 

mu = 0
mu[1] <- 0  # initial condition  
y = 0 

####Local level model iteration

for(t in 1:(n-1)){
  mu[t+1] = mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}

ts.plot(y)
lines(mu,col="red")

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)
ts.plot(eta)
lines(e, col = "red")

# the stochastic properties of: 
# y_t 
# mu_t
# (1-L)y_t 
# remind that they essentially depend on the signal-to-noise ratio 

# let's start with y_t (that from theory is ARIMA(0,1,1) with constraints)
ts.plot(y)
acf(y,lag.max = 100, drop.lag.0 = FALSE)
pacf(y,lag.max = 100, drop.lag.0 = FALSE) #partial autocorrelation function

ts.plot(mu)
acf(mu,lag.max = 100, drop.lag.0 = FALSE)
pacf(mu,lag.max = 100, drop.lag.0 = FALSE) #partial autocorrelation function

# as expected, the ACF of y and mu are slowly linearly decreasing as the ACFs of 
# nonstationary processes 

dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

# remind that (1-L)y_t is a MA(1) with rho_1 = -1/(2+q)
# when q goes to infinity then rho_1 -> 0

# while the reduced form of (1-L)y_t is MA(1),  (1-L)mu_t 
# is a white noise process, actually, BY CONSTRUCTION, is eta_t

dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

################################################################################
# LLT: local linear trend
################################################################################
n = 400

sigma_eta = 0.2
eta = sigma_eta * rnorm(n)
sigma_e = 0.8
e = sigma_e * rnorm(n)
sigma_xi = 2.1
xi = sigma_xi * rnorm(n)

mu = 0
mu[1] <- 0   
beta = 0
beta[1] <- 0   
y = 0 

#local linear trend iteration
for(t in 1:(n-1)){
  beta[t+1] = beta[t] + xi[t]
  mu[t+1] = beta[t] + mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}

ts.plot(y)
lines(mu,col="red")
lines(beta,col="blue")

ts.plot(y)
ts.plot(mu)
ts.plot(beta)


###############################################################################







