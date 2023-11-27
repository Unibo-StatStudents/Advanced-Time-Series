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


a <- read.table("daxdow.txt")

uk <- a[, 1]
deu <- a[,2]

uk_log <- log(uk)
deu_log <- log(deu)

plot(uk_log)
acf(deu_log)

acf(diff(uk_log))
acf(diff(deu_log))


pacf(diff(uk_log))
pacf(diff(deu_log))


y <- 4*3
y

