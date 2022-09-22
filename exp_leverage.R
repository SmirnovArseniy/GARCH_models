library("pROC")
library("margins")
library("boot")
library("lmtest")
library("numDeriv")
library("rJava")
library("xlsx")
library("forecast")
library("rugarch")
library('chron')
require(xts)
rm(list = ls())
options(scipen = 999)
n <- 10000


mu <- 1
omega <- 0.1
alpha <- 0.15
beta <- 0.5
xi <- 0.2
delta <- 0.6
nu1 <- 0.1
nu2 <- -0.5
lambda <- 1

ln.sigma2_0 <- log((omega + xi * alpha) / (1 - beta - delta * alpha))

z <- rnorm(n)                   
eps <- rep(NA, n)               
ln.sigma2 <- rep(NA, n)         
ln.x <- rep(NA, n)           
u <- rnorm(n, sd = lambda)                    
ln.sigma2[1] <- ln.sigma2_0                 
eps[1] <- sqrt(exp(ln.sigma2[1])) * z[1]      
ln.x[1] <- xi + delta * ln.sigma2[1] + 
  nu1 * z[1] + nu2 * (abs(z[1])-sqrt(2/pi)) + u[1] 

for(t in 2:n)
{
  
  ln.sigma2[t] <- omega + alpha * ln.x[t - 1] + beta * ln.sigma2[t - 1]
  ln.x[t] <- xi + delta * ln.sigma2[t] + nu1 * z[t] + nu2 * (abs(z[t])-sqrt(2/pi)) + u[t]
  eps[t] <- sqrt(exp(ln.sigma2[t])) * z[t]
}


y <- mu + eps
RV <- exp(ln.x)
data <- data.frame("y" = y, 'RV' = RV)

lnL <- function (x,                       
                 data,                    
                 is_aggregate = TRUE)      
{ 
  mu <- x[1]                              
  omega <- x[2]
  alpha <- x[3]
  beta <- x[4]
  xi <- x[5]
  delta <- x[6]
  nu1 <- x[7]                            
  nu2 <- x[8]
  lambda <- x[9]
  
  data <- as.matrix(data)
  
  n = length(data[, 1])                       
  y <- data[, 1]                             
  x <- data[, 2]                              
  ln.x <- log(x)
  ln.sigma2 <- rep(NA, n)                 
  z <- rep(NA, n)                          
  
  u <- rep(NA, n)                          
  ln.sigma2_0 <- log((omega + xi * alpha) / (1 - beta - delta * alpha)) 
  ln.sigma2[1] <- ln.sigma2_0                                           
  eps <- y - mu                               
  z_0 <- eps[1] / sqrt(exp(ln.sigma2[1]))    
  z[1] <- z_0
  u[1] <- ln.x[1] - xi - delta * ln.sigma2[1] -nu1 * z[1] - nu2 * (abs(z[1])-sqrt(2/pi)) 
  
  for(t in 2:n)
  {
    ln.sigma2[t] <- omega + alpha * ln.x[t - 1] + beta * ln.sigma2[t - 1]
    z[t] <-  eps[t] / sqrt(exp(ln.sigma2[t]))
    u[t] <- ln.x[t] - xi - delta * ln.sigma2[t] - nu1 * z[t] - nu2 * (abs(z[t])-sqrt(2/pi))
    
  }
  
  sigma <- sqrt(exp(ln.sigma2))
  L_vector <- dnorm(eps, mean = 0, sd = sigma) * dnorm(u, mean = 0, sd = lambda)
  
  
  lnL_value <- log(L_vector)              
  
  
  if(!is_aggregate)                         
  {                                               
    return(lnL_value)                      
  }
  
  return(sum(lnL_value))                   
}                   
x0 <- c(mu, omega, alpha, beta, xi, delta, nu1, nu2, lambda)

lnL(x0, data)

result <- optim(par = x0,                              # Р Р…Р В°РЎвЂЎР В°Р В»РЎРЉР Р…Р В°РЎРЏ РЎвЂљР С•РЎвЂЎР С”Р В° (Р В»РЎС“РЎвЂЎРЎв‚¬Р Вµ РЎвЂљР С•Р В¶Р Вµ Р В°РЎР‚Р С–РЎС“Р СР ВµР Р…РЎвЂљР С•Р С Р С—Р С•Р Т‘Р В°Р Р†Р В°РЎвЂљРЎРЉ, Р Р…Р С• Р Т‘Р В»РЎРЏ Р С—РЎР‚Р С•РЎРѓРЎвЂљР С•РЎвЂљРЎвЂ№ Р С—Р С•Р С”Р В° Р С—РЎР‚РЎРЏР СР С• Р Р† РЎвЂљР ВµР В»Р Вµ РЎвЂћРЎС“Р Р…Р С”РЎвЂ Р С‘Р С‘)
                method = "BFGS",                       # РЎвЂЎР С‘РЎРѓР В»Р ВµР Р…Р Р…РЎвЂ№Р в„– Р СР ВµРЎвЂљР С•Р Т‘ Р С•Р С—РЎвЂљР С‘Р СР С‘Р В·Р В°РЎвЂ Р С‘Р С‘
                fn = lnL,                              # Р СР В°Р С”РЎРѓР С‘Р СР С‘Р В·Р С‘РЎР‚РЎС“Р ВµР СР В°РЎРЏ РЎвЂћРЎС“Р Р…Р С”РЎвЂ Р С‘РЎРЏ Р С—РЎР‚Р В°Р Р†Р Т‘Р С•Р С—Р С•Р Т‘Р С•Р В±Р С‘РЎРЏ
                control = list(maxit = 10000,          # РЎвЂЎРЎвЂљР С•Р В±РЎвЂ№ Р СР С‘Р Р…Р С‘Р СР С‘Р В·Р В°РЎвЂ Р С‘Р С•Р Р…Р Р…РЎС“РЎР‹ Р В·Р В°Р Т‘Р В°РЎвЂЎРЎС“ Р С—РЎР‚Р ВµР Р†РЎР‚Р В°РЎвЂљР С‘РЎвЂљРЎРЉ
                               fnscale = -1,           # Р Р† Р СР В°Р С”РЎРѓР С‘Р СР С‘Р В·Р В°РЎвЂ Р С‘Р С•Р Р…Р Р…РЎС“РЎР‹ РЎС“Р СР Р…Р С•Р В¶Р В°Р ВµР С РЎвЂћРЎС“Р Р…Р С”РЎвЂ Р С‘РЎР‹ Р Р…Р В° -1
                               reltol = 1e-10),        # РЎС“РЎРѓРЎвЂљР В°Р Р…Р С•Р Р†Р С‘Р С Р Т‘Р С•РЎРѓРЎвЂљР В°РЎвЂљР С•РЎвЂЎР Р…Р С• Р Р†РЎвЂ№РЎРѓР С•Р С”РЎС“РЎР‹ РЎвЂљР С•РЎвЂЎР Р…Р С•РЎРѓРЎвЂљРЎРЉ          
                hessian = FALSE,                       # Р Р†Р ВµРЎР‚Р Р…Р ВµР С Р вЂњР ВµРЎРѓРЎРѓР С‘Р В°Р Р… РЎвЂћРЎС“Р Р…Р С”РЎвЂ Р С‘Р С‘
                data = data)   
result
