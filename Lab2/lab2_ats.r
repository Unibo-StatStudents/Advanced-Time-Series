################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
################################################################################

# write detailed comments each line of the following code and improve it whenever you
# you feel able to do that 

n = 300

sigma_eta = 1 #initialize the variance
eta = sigma_eta * rnorm(n) #simulate the eta noise


sigma_e = 1 #initalize the variance
e = sigma_e * rnorm(n) #simulate the epsilin noise



q = sigma_eta^2/sigma_e^2 #signal-to-noise ratio


#if we want to use a t-student distribution for the noise part
#  nu = 3   degree of freedom
#  e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)    standardization of e
#  ts.plot(e)


ts.plot(e)
lines(eta, col = "red")


#set the initial parameters
phi <- 0.8 
mu = 0
mu[1] <- 0   
y = 0 

#iterate the AR(1) signal plus noise model
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

ts.plot(y)
lines(mu,col="red")


acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

#' initialize the parameters of the kalman filter
mu_pred <- 0   #prediction of the mean
P <- 0         #variance
v <-0          #error measurament betwwen the observed value and the predicted
K <- 0         #kalman gain
F <- 0         #
llk<-0         #likelihood


mu_pred[1] = 0  #initialize the mu_pred parameter
P[1] = (sigma_eta^2)/(1-phi^2)  #autocovariance at the first lag
llk[1]=0   #initialize the likelihood


for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

