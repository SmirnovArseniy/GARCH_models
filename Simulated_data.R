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

########################## ФУНКЦИЯ РЫЧАГА ########################################
func_value <- function(type, tau1, xi, tau2=NULL, ln_sigma2=NULL) {
  if (type == 'polynomial') {
    return(tau1 * xi + tau2 * (xi ^ 2 - 1)) 
  } else if (type == 'egarch') {
    return(tau1 * xi + tau2 * (abs(xi) - sqrt(2 / pi)))
  } else if (type == 'gjr') {
    return(tau1 * xi ^ 2 * exp(ln_sigma2) * ifelse(xi < 0, xi, 0))
  }
}
#######################################################################################
n <- 10000


params = list(
  mu=0.5, 
  omega = 0.001,
  beta=0.02,
  gamma=0.01,
  dz=0.005,
  phi=-0.15,
  tau1=0.009,
  tau2=-0.0036,
  sigma_u=0.04,
  sigma0= 0.1
)

####################### ГЕНЕРАЦИЯ ДАННЫХ #################################################

generate <- function(type, specification, params, n_observation = 10 ^ 5, random_state=2001) {
  mu <- params$mu
  omega <- params$omega
  beta <- params$beta
  gamma <- params$gamma
  dz <- params$dz
  phi <- params$phi
  tau1 <- params$tau1
  tau2 <- params$tau2
  sigma_u <- params$sigma_u
  ln_sigma2_0 <- log(params$sigma0^2)
  set.seed(random_state)
  xi <- rnorm(n_observation, mean=0, sd=1)
  u <- rnorm(n_observation, mean=0, sd=sigma_u) 
  
  ln_sigma2 <- rep(NaN, times = n_observation)
  eps <- rep(NaN, times = n_observation)
  ln_x <- rep(NaN, times = n_observation)
  
  if (is.null(ln_sigma2_0)) {
    ln_sigma2_0 <- log((omega + gamma * dz)/(1 - (beta + phi * gamma)))
  }
  if (specification == 'log') {
    ln_sigma2[1] <- ln_sigma2_0
    eps[1] <- sqrt(exp(ln_sigma2[1])) * xi[1]
    ln_x[1] <- dz + phi * ln_sigma2[1] + func_value(type=type, tau1=tau1, tau2=tau2, ln_sigma2=ln_sigma2[1], xi=xi[1]) + u[1]
    
    for (t in 2:n_observation) {
      ln_sigma2[t] <-omega + beta * ln_sigma2[t - 1] + gamma * ln_x[t - 1]
      ln_x[t] <- dz + phi * ln_sigma2[t] + func_value(type=type, tau1=tau1, tau2=tau2, ln_sigma2=ln_sigma2[t], xi=xi[t]) + u[t]
      eps[t] <- sqrt(exp(ln_sigma2[t])) * xi[t]
    }
    return(list(y= eps+mu,
                x = exp(ln_x)))
  }
}

gen <- generate(type = 'polynomial',specification = 'log', params = params, n_observation = 10^5)
########################################################################################################################

###################################### ЗАДАДИМ ПАРАМЕТРЫ В ВЕКТОРЕ ######################################################
params_c <- c(0.3, 0.02, 0.05, 0.03, 0.01, -0.001, 0.002, 0.00004, 0.1, 1.0)

#########################################################################################################################

####################################################### СДЕЛАЕМ ТАБЛИЦУ СДАННЫМИ ########################################
data <- data.frame('y'=gen$y,'ln_x'=gen$x)
data$y <- as.numeric(data$y)
data$ln_x <- as.numeric(data$ln_x)
###########################################################################################################################

################################################# ФУНКЦИЯ ПРАДОПОДОБИЯ ####################################################
log_likelihood_pol <- function(params, data, is_aggregate = TRUE)
{

  mu <- params_c[1]
  omega <- params_c[2]
  beta <- params_c[3]
  gamma <- params_c[4]
  dz <- params_c[5]
  phi <- params_c[6]
  tau1 <- params_c[7]
  tau2 <- params_c[8]
  sigma_u <- params_c[9]
  ln_sigma2_0 <- log(params_c[10]^2)
  
  data <- as.matrix(data)
  
  n = length(data[, 1])                       # Р С”Р С•Р В»Р С‘РЎвЂЎР ВµРЎРѓРЎвЂљР Р†Р С• Р Р…Р В°Р В±Р В»РЎР‹Р Т‘Р ВµР Р…Р С‘Р в„–
  y <- data[, 1]                              # Р Т‘Р С•РЎвЂ¦Р С•Р Т‘Р Р…Р С•РЎРѓРЎвЂљР С‘
  ln_x <- data[, 2]   
  eps <- y - mu
  ln_sigma2 <- rep(NaN, times = n)
  
  ln_sigma2_0 <- log((omega + gamma * dz) / (1 - (beta + phi * gamma)))
  
  ln_sigma2[1] = ln_sigma2_0
  
  xi <- rep(NaN, times = n)
  u <- rep(NaN, times = n)

  xi_0 <- eps[1] /sqrt(exp(ln_sigma2_0))
  
  xi[1] <- xi_0

  u[1] <- ln_x[1] - dz - phi * ln_sigma2[1] - (tau1 * xi[1] + tau2 * (xi[1] ^ 2 - 1))

  for (t in 2:n) 
    {
      ln_sigma2[t] <- omega + beta * ln_sigma2[t - 1] + gamma * ln_x[t - 1]
      xi[t] <- eps[t] / sqrt(exp(ln_sigma2[t]))
      u[t] <- ln_x[t] - dz - phi * ln_sigma2[t] - tau1 * xi[t] - tau2 * (xi[t] ^ 2 - 1)
    }
  sigma <- sqrt(exp(ln_sigma2))
  L <- dnorm(eps, mean = 0, sd = sigma) * dnorm(u, mean = 0, sd = sigma_u)
  return(sum(log(L)))
}
##################################################################################################
log_likelihood_pol(params_c, data)
############################# МАКСИМИЗИРУЕМ ФУНКЦИЮ ПРАВДОПОДОБИЯ ###############################
data <- as.matrix(data)
result <- optim(par = params_c,       
                method = "BFGS",                
                fn = log_likelihood_pol,                              
                control = list(maxit = 10000,          
                               fnscale = -1,           
                               reltol = 1e-10),                
                hessian = FALSE,                      
                data = data) 
result

##################### ЮРААААААААААААААААААААААААААААААААААААААА#############################
