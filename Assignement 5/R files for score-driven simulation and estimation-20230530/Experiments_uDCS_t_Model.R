##################################################################################
# this file does fit a DCS-t location model to data that ARE SIMULATED
# it uses the functions defined in Simulator, Estimator, Filter
##################################################################################

T <- 600
omega <- 0
phi <- 0.8
k <- 0.8
varsigma <- 1
nu <- 8

theta <- c(omega, phi, k, varsigma, nu)

simu <- uDCS_t_model_simulator(T, omega, phi, k, varsigma, nu)
y = simu$y_t_gen
ts.plot(y)
# forecast::autoplot(simu$y_t_gen)

est  <- uDCS_t_model_estimator(y, theta)
est$theta_list

filter <- uDCS_t_model_filter(y, est$theta)
ts.plot(y)
lines(filter$Dynamic_Location,col = "red")
#autoplot(y) + autolayer(filter$Dynamic_Location)
