#Lab2_ats
################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
################################################################################

# write detailed comments each line of the following code and improve it whenever you
# you feel able to do that 

n = 300

sigma_eta = 0.1
eta = sigma_eta * rnorm(n)


sigma_e = 1

q = sigma_eta^2/sigma_e^2 #why we need to define q?

e = sigma_e * rnorm(n)

#nu = 3
#e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
#ts.plot(e)


ts.plot(e)
lines(eta, col = "red")

phi <- 0.8
mu = 0
mu[1] <- 0   
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)
mean(y)
# y = mu + e

ts.plot(y)
lines(mu,col="red")


acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)


n = 300

sigma_eta = 1
eta = sigma_eta * rnorm(n)



#sigma_e = 1

q = sigma_eta^2/sigma_e^2 #why we need to define q?

#e = sigma_e * rnorm(n)

nu = 3


ts.plot(e)
lines(eta, col = "red")

phi <- 0.8
mu = 0
mu[1] <- 0   
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e
mean(y)
ts.plot(y)
lines(mu,col="red")


acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################


mu_pred <- 0   
P <- 0         
v <-0          
K <- 0         
F <- 0           
llk<-0       


mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)  #\sigma^2_{epsilon} 
                                #why \phi and not \theta?
                                #where P[1] come from?
llk[1]=0


for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  #llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')
