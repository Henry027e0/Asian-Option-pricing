set.seed(456)
library(quantmod)
library(RQuantLib)
library(ggplot2)


getSymbols("CL=F", from = Sys.Date() - 365*20, to = Sys.Date())
wti_df <- data.frame(Date = index(`CL=F`), Price = as.numeric(Ad(`CL=F`)))
wti_df$LogReturn <- c(NA, diff(log(wti_df$Price)))
log_returns <- na.omit(c(NA, diff(log(wti_df$Price))))

#VOLATILITY AND DRIFT
S0 <- tail(wti_df$Price, 1)
volatility <- sd(log_returns)
drift <- mean(log_returns) + 0.5 * volatility^2

#RISK-FREE RATE
getSymbols("DGS10", src = "FRED") 
long_rate <- as.numeric(tail(DGS10, 1)) / 100
getSymbols("SOFR", src = "FRED")
current_sofr <- as.numeric(tail(SOFR, 1)) / 100


#PARAMETERS
K <- 75
expiry_date <- as.Date("2025-10-31")
current_date <- Sys.Date()
days_to_maturity <- as.numeric(expiry_date - current_date)
trading_days <- 252
T <- days_to_maturity / trading_days
n_steps <- max(1, floor(T * trading_days))
dt <- T / n_steps  

#RISK-FREE DECISION
r_T <- ifelse(T <= 0.25, current_sofr,
              ifelse(T >= 10, long_rate,
                     current_sofr + (long_rate - current_sofr) * (T / 10)))



#PARAMETERS
sigma_annual <- volatility * sqrt(trading_days)  #annual vol
r_annual_adj <- drift
total_sim <- 100000
n_pilot <- 10000
n_main <- total_sim - n_pilot

#PILOT SIMULATION
pilot_A <- numeric(n_pilot)
pilot_G <- numeric(n_pilot)
pilot_Y <- numeric(n_pilot)

for (i in 1:n_pilot) {
  epsilon <- rnorm(n_steps)
  drift_daily <- (r_T - 0.5 * sigma_annual^2) * dt
  noise_daily <- sigma_annual * sqrt(dt) * epsilon
  log_returns <- drift_daily + noise_daily
  
  S_path <- S0 * exp(cumsum(log_returns))
  S_path <- c(S0, S_path)
  
  pilot_A[i] <- mean(S_path)
  pilot_G[i] <- exp(mean(log(S_path + 1e-8))) #1e-8 to avoid null arg in log
  pilot_Y[i] <- exp(-r_T * T) * max(pilot_A[i] - K, 0)
}

#KEMNA VORST CONTROL VARIATE
n <- n_steps
sigma_adj <- sigma_annual * sqrt((2*n + 1)/(6*(n + 1)))
mu_adj <- (r_T - 0.5*sigma_annual^2) * (n + 1)/(2*n) + 0.5*sigma_adj^2
d1 <- (log(S0/K) + (mu_adj + 0.5*sigma_adj^2)*T) / (sigma_adj * sqrt(T))
d2 <- d1 - sigma_adj * sqrt(T)
E_geo <- exp(-r_T * T) * (S0 * exp(mu_adj * T) * pnorm(d1) - K * pnorm(d2))

#OPTIMAL COEFFICIENT
cov_YG <- cov(pilot_Y, pilot_G)
var_G <- var(pilot_G)
c_star <- -cov_YG / var_G

#MAIN SIMULATION
main_Z <- numeric(n_main)

for (i in 1:n_main) {
  epsilon <- rnorm(n_steps)
  drift_daily <- (r_T - 0.5 * sigma_annual^2) * dt
  noise_daily <- sigma_annual * sqrt(dt) * epsilon
  log_returns <- drift_daily + noise_daily
  
  S_path <- S0 * exp(cumsum(log_returns))
  S_path <- c(S0, S_path)
  
  A_i <- mean(S_path)
  G_i <- exp(mean(log(S_path + 1e-10)))
  
  Y_i <- exp(-r_T * T) * max(A_i - K, 0)
  Y_geo_i <- exp(-r_T * T) * max(G_i - K, 0)
  
  main_Z[i] <- Y_i + c_star * (Y_geo_i - E_geo)
}

option_price <- mean(main_Z)
option_price

#PLOT

S0 <- 64
mu <- drift
sigma <- volatility   
T <- 50/252       
N <- 50            
dt <- T/N          
n_paths <- 50      


paths <- matrix(nrow = N+1, ncol = n_paths)
paths[1,] <- S0    


for (i in 1:n_paths) {
  for (t in 2:(N+1)) {
    epsilon <- rnorm(1)
    paths[t,i] <- paths[t-1,i] * exp((mu - 0.5*sigma^2)*dt + sigma*sqrt(dt)*epsilon)
  }
}


time <- seq(0, T, length.out = N+1)


plot(time, paths[,1], type = "l", col = rainbow(n_paths)[1],
     xlab = "Time", ylab = " WTI Price", main = "GBM Paths", ylim = c(0.95*S0, 1.05*S0))
for (i in 2:n_paths) {
  lines(time, paths[,i], col = rainbow(n_paths)[i])
}



