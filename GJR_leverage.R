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
  nu1 * z[1] ^ 2 * exp(ln.sigma2[1]) * ifelse(z[1] < 0, z[1], 0) + u[1] 

for(t in 2:n)
{
  
  ln.sigma2[t] <- omega + alpha * ln.x[t - 1] + beta * ln.sigma2[t - 1]
  ln.x[t] <- xi + delta * ln.sigma2[t] + nu1 * z[t] ^ 2 * exp(ln.sigma2[t]) * ifelse(z[t] < 0, z[t], 0) + u[t]
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
  u[1] <- ln.x[1] - xi - delta * ln.sigma2[1] -nu1 * z[1] ^ 2 * exp(ln.sigma2[1]) * ifelse(z[1] < 0, z[1], 0)
  
  for(t in 2:n)
  {
    ln.sigma2[t] <- omega + alpha * ln.x[t - 1] + beta * ln.sigma2[t - 1]
    z[t] <-  eps[t] / sqrt(exp(ln.sigma2[t]))
    u[t] <- ln.x[t] - xi - delta * ln.sigma2[t] - nu1 * z[t] ^ 2 * exp(ln.sigma2[t]) * ifelse(z[t] < 0, z[t], 0)
    
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

result <- optim(par = x0,                             
                method = "BFGS",                      
                fn = lnL,                              
                control = list(maxit = 10000,         
                               fnscale = -1,          
                               reltol = 1e-10),       
                hessian = FALSE,                      
                data = data)   
result
