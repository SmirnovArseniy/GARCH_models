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
options(scipen = 999)

              #####################################################
              # Процесс симуляции данных согласно модели RV-GARCH #
              #####################################################

# число наблюдений
n <- 10000

# истинные параметры GARCH модели
mu <- 1
omega <- 0.1
alpha <- 0.15
beta <- 0.5
xi <- 0.2
delta <- 0.6
nu1 <- 0.1
nu2 <- -0.5
lambda <- 2

# волатильность в начальный момент времени (безусловная дисперсия)
ln.sigma2_0 <- log((omega + xi * alpha) / (1 - beta - delta * alpha))

z <- rnorm(n)                    # независимые стандартные нормальные (шоки в доходности стандартизированные)
eps <- rep(NA, n)                # случайная ошибка
ln.sigma2 <- rep(NA, n)          # волатильность
ln.x <- rep(NA, n)               # реализованная волатильность (вектор rt)
u <- rnorm(n, sd = lambda)                    # случайная ошибка в уравнении RV
ln.sigma2[1] <- ln.sigma2_0                   # волатильность для первого момента времени
eps[1] <- sqrt(exp(ln.sigma2[1])) * z[1]      # случайная ошибка для первого момента времени
ln.x[1] <- xi + delta * ln.sigma2[1] + 
           nu1 * z[1] + nu2 * (z[1] ^ 2 - 1) + u[1]  # значение RV для первого момента времени

# Считаем для всех остальных моментов времени
for(t in 2:n)
{
  
  ln.sigma2[t] <- omega + alpha * ln.x[t - 1] + beta * ln.sigma2[t - 1]
  ln.x[t] <- xi + delta * ln.sigma2[t] + nu1 * z[t] + nu2 * (z[t] ^ 2 - 1) + u[t]
  eps[t] <- sqrt(exp(ln.sigma2[t])) * z[t]
}

# создадим доходности 
y <- mu + eps
RV <- exp(ln.x)

data <- data.frame("y" = y, 'RV' = RV) # создадим выборку с доходностями и значениями реализованной волатильности 

# переведем данные в формат xts (рандромно даты расставим для симуляций)
date <- seq.Date(from = as.Date('2017-01-01'), to = as.Date('2019-09-27'), by = 'days') # base
data_time <- data.frame(as.POSIXct(date), data)
data_xts <- xts(data_time[-1], order.by=data_time$as.POSIXct.date.)

                  ###################################################
                               # Оценим пакетом модель #
                  ###################################################

spec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = TRUE), variance.model = list(model = 'realGARCH', garchOrder = c(1, 1)))
fit = ugarchfit(spec, data_xts$y, solver = 'hybrid', realizedVol = data_xts$RV)
cbind(est = coef(fit), true = c(mu, omega, alpha, beta, nu1, nu2, delta, lambda, xi))
package_est_vector <- c(coef(fit)[1:4],coef(fit)[9], coef(fit)[7], coef(fit)[5:6], coef(fit)[8])     # сохраним оценки из пакета в нужном нам порядке

                  ###################################################
                       # Запрограммируем функцию правдоподобия #
                  ###################################################

lnL <- function (x,                        # вектор параметров
                 data,                     # реализация выборки
                 is_aggregate = TRUE)      # возвращаем сумму вкладов наблюдений в функцию  
                                           # правдоподобия (TRUE) или по отдельности (FALSE)
{ 
  mu <- x[1]                               # для удобства создадим отдельные переменные под значения оцениваемых параметров
  omega <- x[2]
  alpha <- x[3]
  beta <- x[4]
  xi <- x[5]
  delta <- x[6]
  nu1 <- x[7]                             # параметр, отвечающий за эффект рычага
  nu2 <- x[8]
  lambda <- x[9]
  
  data <- as.matrix(data)
  
  n = length(data[, 1])                       # количество наблюдений
  y <- data[, 1]                              # доходности
  x <- data[, 2]                              # реализованная волатильность
  ln.x <- log(x)
  ln.sigma2 <- rep(NA, n)                  # волатильность
  z <- rep(NA, n)                          # стандартные нормальные случайные величины, определяющие эпсилон
  
  u <- rep(NA, n)                          # случайные ошибки уравнения RV
  ln.sigma2_0 <- log((omega + xi * alpha) / (1 - beta - delta * alpha)) # моделируем логарифм сигмы в квадрате для первого момента
  ln.sigma2[1] <- ln.sigma2_0                                           
  eps <- y - mu                               # определяем случайные шоки epsilon
  z_0 <- eps[1] / sqrt(exp(ln.sigma2[1]))     # определяем значение стандартной нормальной случайной величины в первый момент времени
  z[1] <- z_0
  u[1] <- ln.x[1] - xi - delta * ln.sigma2[1] - nu1 * z[1] - nu2 * (z[1] ^ 2 - 1) # значение случайной ошибки уравнения RV для первого момента времени
  
  for(t in 2:n)
  {
    # рекурсивные расчеты
    ln.sigma2[t] <- omega + alpha * ln.x[t - 1] + beta * ln.sigma2[t - 1]
    z[t] <-  eps[t] / sqrt(exp(ln.sigma2[t]))
    u[t] <- ln.x[t] - xi - delta * ln.sigma2[t] - nu1 * z[t] - nu2 * (z[t] ^ 2 - 1)
    
  }
  
  sigma <- sqrt(exp(ln.sigma2))
  L_vector <- dnorm(eps, mean = 0, sd = sigma) * dnorm(u, mean = 0, sd = lambda)
  
  
  lnL_value <- log(L_vector)               # считаем значение логарифма 
                                           # функции правдоподобия
  
  
  if(!is_aggregate)                        # возвращаем значение логарифма 
  {                                        # функции правдоподобия при        
    return(lnL_value)                      # заданных значениях параметров
  }
  
  return(sum(lnL_value))                   # возвращаем значения для каждого
                                           # наблюдения по отдельности
}
x0 <- c(mu, omega, alpha, beta, xi, delta, nu1, nu2, lambda)
lnL(x = x0, data = data)                   # значение функции правдоподобия в точках истинных параметров
lnL(x = package_est_vector, data = data)   # значение функции правдоподобия в точках оценок пакета для сравнения
likelihood(fit)                            # значение правдоподобия из пакета rugarch

                ###################################################
                # Теперь необходимо запрограммировать оптимизатор #
                ###################################################

GARCH <- function(data, x0)                              # в функцию подаются сгенерированные данные
  
{
  data <- as.matrix(data)                                # данные
  
  result <- optim(par = x0,                              # начальная точка (лучше тоже аргументом подавать, но для простоты пока прямо в теле функции)
                  method = "BFGS",                       # численный метод оптимизации
                  fn = lnL,                              # максимизируемая функция правдоподобия
                  control = list(maxit = 10000,          # чтобы минимизационную задачу превратить
                                 fnscale = -1,           # в максимизационную умножаем функцию на -1
                                 reltol = 1e-10),        # установим достаточно высокую точность          
                  hessian = FALSE,                       # вернем Гессиан функции
                  data = data)                           # аргументы оптимизируемой функции 
  
  # Для получения более точного результата
  # дополнительно используем генетический
  # алгоритм глобальной оптимизации
  ga_result <- ga(
    type = "real-valued",                   # воспринимаем оптимизируемые параметры
    # как вещественные числа
    fitness = lnL,                          # максимизируемая функция
    lower = -abs(result$par) * 10,          # векторы верхних и нижних границ,
    upper = abs(result$par) * 10,           # в которых ищутся оптимальные значения параметров
    suggestions = result$par,               # начальная точка
    seed = "8",                             # в целях воспроизводимости
    optim = TRUE,
    optimArgs = list(method = "Nelder-Mead", 
                     poptim = 0.05,
                     pressel = 0.5,
                     control = list(maxit = 1000,
                                    fnscale = -1,
                                    reltol = 1e-10)),
    data = data,
    maxiter = 100)                             
  ga_result_summary <- summary(ga_result)
  ga_result_summary$solution
  
  gamma_est <- ga_result_summary$solution                # оценки коэффициентов
  
  
  return_list <- list("gamma" = gamma_est,               # возвращаем оценки коэффициентов и
                      # асимптотической ковариационной матрицы
                      "data" = data,                     # возвращаем использовавшийся датафрейм
                      "lnL" = result$value)              # возвращаем логарифм функции правдоподобия
  
  
  class(return_list) <- "GARCH"                          # для удобства назначим класс
  # возвращаемой из функции переменной
  return(return_list)                                    # возвращаем результат                               
}
get_res <- GARCH(data = data, x0 = x0)
get_res$gamma

# Сравним полученные оценки с истинными значениями параметров и оценками rugarch
cbind(est = as.vector(get_res$gamma), true = c(mu, omega, alpha, beta, xi, delta, nu1, nu2, lambda), rugarch = package_est_vector)
