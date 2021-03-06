# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied econometrics 5, Copulas for portfolio modeling Орлов
# Евгений, 21.11.2014

# ДАННЫЕ ####
gazp.df <- read.csv(file.path("data", "ГАЗПРОМ ао [Price].txt"), header = TRUE)
gmkn.df <- read.csv(file.path("data", "ГМКНорНик [Price].txt"), header = TRUE)
tail(gazp.df)
# Цены закрытия начиная с 01.01.2010 по 07.11.2014
# включительно
gazp.close <- gazp.df[gazp.df$X.DATE. >= 20100101, "X.CLOSE."]
gmkn.close <- gmkn.df[gmkn.df$X.DATE. >= 20100101, "X.CLOSE."]
# Даты
gazp.dates <- gazp.df[gazp.df$X.DATE. >= 20100101, "X.DATE."]
gazp.dates <- as.Date(as.character(gazp.dates), format = "%Y%m%d")
gmkn.dates <- gmkn.df[gmkn.df$X.DATE. >= 20100101, "X.DATE."]
gmkn.dates <- as.Date(as.character(gmkn.dates), format = "%Y%m%d")
date.comp <- gazp.dates == gmkn.dates
if (sum(date.comp)/length(date.comp) == TRUE) {
  print("Даты совпадают!")
} else {
  print("Даты не совпадают!")
}
# Размер выборок
samp.len <- length(gazp.close)
print(paste0("Размер выборки: ", samp.len))
# Векторы доходностей
gazp.ret <- gazp.close[2:samp.len]/gazp.close[1:(samp.len - 1)] - 1
gmkn.ret <- gmkn.close[2:samp.len]/gmkn.close[1:(samp.len - 1)] - 1

# Графики цен закрытия
plot(gazp.dates, gazp.close, type = "l", xlab = "index", ylab = "close", main = "Gazprom:Close")
plot(gmkn.dates, gmkn.close, type = "l", xlab = "index", ylab = "close", main = "NorilskNickel:Close")
# Графики доходностей
plot(gazp.dates[-1], gazp.ret, type = "l", xlab = "index", ylab = "return", main = "Gazprom:Return")
plot(gmkn.dates[-1], gmkn.ret, type = "l", xlab = "index", ylab = "return", main = "NorilskNickel:Return")

# МОДЕЛИРОВАНИЕ КОПУЛЫ ####

library("ghyp", quietly = FALSE)
# моделирование частных функций
# распределения
gazp.fit <- stepAIC.ghyp(gazp.ret, dist = c("gauss", "t", "ghyp"), symmetric = NULL, 
  silent = TRUE)$best.model
gmkn.fit <- stepAIC.ghyp(gmkn.ret, dist = c("gauss", "t", "ghyp"), symmetric = NULL, 
  silent = TRUE)$best.model
# расчёт значений F1(u) и F2(u)
gazp.cdf <- pghyp(gazp.ret, object = gazp.fit)
gmkn.cdf <- pghyp(gmkn.ret, object = gmkn.fit)
cdf <- array(c(gazp.cdf, gmkn.cdf), dim = c(length(gazp.ret), 2))

# МОДЕЛИРОВАНИЕ КОПУЛЫ ####
library("copula", quietly = TRUE)
# объявление копулы
norm.cop <- normalCopula(dim = 2, param = 0.5, dispstr = "un")
stud5.cop <- tCopula(dim = 2, param = 0.5, df = 5, df.fixed = TRUE, dispstr = "un")
stud10.cop <- tCopula(dim = 2, param = 0.5, df = 10, df.fixed = TRUE, dispstr = "un")
stud20.cop <- tCopula(dim = 2, param = 0.5, df = 20, df.fixed = TRUE, dispstr = "un")
gumb.cop <- gumbelCopula(dim = 2, param = 2)
clay.cop <- claytonCopula(dim = 2, param = 2)
# подгонка копулы
norm.fit <- fitCopula(cdf, copula = norm.cop)
stud5.fit <- fitCopula(cdf, copula = stud5.cop)
stud10.fit <- fitCopula(cdf, copula = stud10.cop)
stud20.fit <- fitCopula(cdf, copula = stud20.cop)
gumb.fit <- fitCopula(cdf, copula = gumb.cop)
clay.fit <- fitCopula(cdf, copula = clay.cop)

# РАСЧЕТ VaR, ES НА ВСЕЙ ВЫБОРКЕ С ПОМОЩЬЮ
# КОПУЛЫ ####

# выбор оптимальной копулы
norm.fit@loglik
stud5.fit@loglik
stud10.fit@loglik
stud20.fit@loglik
gumb.fit@loglik
clay.fit@loglik
AIC(norm.fit)
AIC(stud5.fit)
AIC(stud10.fit)
AIC(stud20.fit)
AIC(gumb.fit)
AIC(clay.fit)

# значения частных функций распределения
N <- 10^3
stud5.sim <- rCopula(n = N, copula = stud5.fit@copula)

# доходности активов
gazp.sim <- qghyp(stud5.sim[, 1], object = gazp.fit)
gmkn.sim <- qghyp(stud5.sim[, 2], object = gmkn.fit)
w <- c(0.5, 0.5)
port.sim <- w[1] * gazp.sim + w[2] * gmkn.sim

# измерители риска
alpha <- 0.1
port.sim <- sort(port.sim)
VaR <- port.sim[alpha * N]
ES <- mean(port.sim[1:(alpha * N - 1)])

# КРИВАЯ VaR ####

# Параметры для расчета VaR
T1 <- 500
T2 <- length(gazp.ret) - T1
alpha <- 0.05
N <- 10^4
w <- c(0.5, 0.5)  # Веса активов в портфеле
VaR.copula <- numeric()

h <- T1  # Длина обучающей выборки
for (i in (T1 + 1):(T1 + T2)) {
  print((i - T1)/T2)
  # Обучающая выборка
  gazp.train <- gazp.ret[(i - h):(i - 1)]
  gmkn.train <- gmkn.ret[(i - h):(i - 1)]
  # Подгон частных функций распределения
  gazp.fit <- stepAIC.ghyp(gazp.train, dist = c("gauss", "t"), symmetric = TRUE, 
    silent = TRUE)$best.model
  gmkn.fit <- stepAIC.ghyp(gmkn.train, dist = c("gauss", "t"), symmetric = TRUE, 
    silent = TRUE)$best.model
  gazp.cdf <- pghyp(gazp.train, object = gazp.fit)
  gmkn.cdf <- pghyp(gmkn.train, object = gmkn.fit)
  cdf <- array(c(gazp.cdf, gmkn.cdf), dim = c(length(gazp.ret), 2))
  # Подгон копулы
  stud5.fit <- fitCopula(cdf, copula = stud5.cop)
  # Симуляция доходности портфеля
  stud5.sim <- rCopula(n = N, copula = stud5.fit@copula)
  gazp.sim <- qghyp(stud5.sim[, 1], object = gazp.fit)
  gmkn.sim <- qghyp(stud5.sim[, 2], object = gmkn.fit)
  port.sim <- w[1] * gazp.sim + w[2] * gmkn.sim
  port.sim <- sort(port.sim)
  # Расчет VaR
  VaR.copula[i - T1] <- port.sim[alpha * N]
}
# График кривой VaR
port.test <- w[1] * gazp.ret[(T1 + 1):(T1 + T2)] + w[2] * gmkn.ret[(T1 + 1):(T1 + 
  T2)]
plot(port.test, type = "l", main = "VaR curve")
lines(VaR.copula, col = "blue", lwd = 2)
legend("bottomleft", c("Returns", "VaR stud5 copula"), col = c("black", "blue"), 
  lty = c(1, 1), lwd = c(1, 2))

# Частота пробоев
kupiec.test <- function(ret, VaR, alpha) {
  # Тест Купика: H0: модельная и эмпирическая
  # частоты пробоя VaR совпадают
  K <- sum(ret < VaR)
  T2 <- length(ret)
  alpha0 <- K/T2
  S <- -2 * log((1 - alpha)^(T2 - K) * alpha^K) + 2 * log((1 - alpha0)^(T2 - K) * 
    alpha0^K)
  p.value <- 1 - pchisq(S, df = 1)
  return(c(alpha0, p.value))
}

# Тест VaR-кривой на частоту пробоев
print(paste0("Kupiec test, alpha = ", alpha))
kup <- kupiec.test(port.test, VaR.copula, alpha)
print(paste0("GPD: alpha0 = ", round(kup[1], 6), ", p-value  = ", round(kup[2], 6))) 
