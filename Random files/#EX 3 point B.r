#EX 3 point B

#simulate the following AR(1) plus noise process for \phi = 0.8
#initialize the population
n <- 3000
#simulate epsilon
sigma_epsilon = 0.6
epsilon = sigma_epsilon*rnorm(n)

#simulate eta
sigma_eta = 0.4
eta = sigma_eta * rnorm(n)

mu = 0
mu[1] <- 0
phi <- 0.8
y = 0
y[1]<-0

for (t in 1:(n-1)){
  mu[t+1]= phi*mu[t] + eta[t]
  y[t] = mu[t] + epsilon[t]
}

ts.plot(y)
mean(y)
var(y)
pacf(y, lag.max = 100, drop.lag.0 = FALSE)








