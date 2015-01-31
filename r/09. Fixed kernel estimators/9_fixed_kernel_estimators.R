# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied econometrics 9, nonparametric modeling Part 1. Fixed kernel estimators
# Орлов Евгений, 22.01.2015

# НУЖНО ДОПОЛНИТЕЛЬНОЕ ТЕСТИРОВАНИЕ

# ДАННЫЕ ####
start.date <- 20100101
rts.df <- read.csv(file.path("data", "RTSI [Price].txt"), header = TRUE)
tail(rts.df)
# Цены закрытия начиная с 01.01.2010 по 07.11.2014
# включительно
rts.close <- rts.df[rts.df$X.DATE. >= start.date, "X.CLOSE."]
# Вектор доходностей
rcl <- length(rts.close)
print(paste0("Размер выборки: ", rcl))
rts.ret <- rts.close[2:rcl]/rts.close[1:(rcl - 1)] - 1

rts.dates <- rts.df[rts.df$X.DATE. >= start.date, "X.DATE."]
rts.dates <- as.Date(as.character(rts.dates), format = "%Y%m%d")
# График цен закрытия
plot(rts.dates, rts.close, type = "l", xlab = "index", ylab = "close", main = "RTS Index:Close")
# График доходностей
plot(rts.dates[-1], rts.ret, type = "l", xlab = "index", ylab = "return", main = "RST Index:Return")


# Построение гистограммы ####
hist(rts.ret, probability = TRUE)

# Простая непараметрическая оценка
# плотности ####

L <- 10^3
N <- length(rts.ret)
h <- 0.01  # ширина интервала

# В точках x будет оцениваться плотность
x <- seq(-0.15, 0.1, length = L)
f.naive <- numeric()  # нулевой (пока) вектор оценок

# Считаем количество элементов в интервалах
# xi ± h/2
for (i in 1:L) {
  f.naive[i] <- sum(1 * ((rts.ret > x[i] - h/2) & (rts.ret < x[i] + h/2)))
}
f.naive <- f.naive/(N * h)  # нормируем оценку

# График простой оценки
plot(x, f.naive, type = "l", main = "Naive estimate", xlab = "y", ylab = "Density")
rug(rts.ret)


# Ядерные оценки ####
library("np", quietly = TRUE)

f.fix <- npudens(tdat = rts.ret, edat = x, ckertype = "gaussian", bwtype = "fixed")
plot(x, f.fix$dens, type = "l", main = "Gaussian kernel, fixed bandwidth", xlab = "y", 
  ylab = "Density")


# Поиск квантиля оценки распределения ####

# Метод деления пополам
q.bisec <- function(CDF, alpha = 0.05) {
  a <- 1
  b <- CDF$nobs
  ab <- trunc((a + b)/2)
  while ((b - a) > 2) {
    if (CDF$dist[ab] <= alpha) {
      a <- ab
    }
    if (CDF$dist[ab] >= alpha) {
      b <- ab
    }
    ab <- trunc((a + b)/2)
  }
  q.fix <- x[ab]
  # print(q.fix)
}
F.fix <- npudist(tdat = rts.ret, edat = x, ckertype = "gaussian", bwtype = "fixed")
q.bisec(F.fix, 0.05)

# Метод Монте-Карло
q.mc <- function(f.fix, alpha = 0.05, M = 10^6) {
  y.fix.sim <- sample(x, prob = f.fix$dens, size = M, replace = TRUE)
  q.fix <- sort(y.fix.sim)[alpha * M]
  # print(q.fix)
}
f.fix <- npudens(tdat = rts.ret, edat = x, ckertype = "gaussian", bwtype = "fixed")
q.mc(f.fix, 0.05, M)

# Кривая VaR, поиск квантиля методом деления
# пополам ####
alpha <- 0.05
T1 <- length(rts.ret)/2
T2 <- length(rts.ret) - T1
VaR.np1 <- numeric()
h <- T1

for (i in (T1 + 1):(T1 + T2)) {
  print((i - T1)/T2)
  # Training set
  rts.train <- rts.ret[(i - h):(i - 1)]
  # Ядерная оценка
  F.fix <- npudist(tdat = rts.train, edat = x, ckertype = "gaussian", bwtype = "fixed")
  # Поиск квантиля методом деления пополам
  VaR.np1[i - T1] <- q.bisec(F.fix, alpha)
}
# График кривой VaR
rts.test <- rts.ret[(T1 + 1):(T1 + T2)]
plot(rts.test, type = "l", main = "VaR curve")
lines(VaR.np1, col = "blue", lwd = 2)
legend("bottomleft", c("Returns", "VaR, Fixed Kernel, Bisection"), col = c("black", 
  "blue"), lty = c(1, 1), lwd = c(1, 2))


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
kup <- kupiec.test(rts.test, VaR.np1, alpha)
print(paste0("NP, Fixed kernel, ver1: alpha0 = ", round(kup[1], 6), ", p-value  = ", 
  round(kup[2], 6)))


# Кривая VaR, поиск квантиля методом
# Монте-Карло ####
alpha <- 0.05
M <- 10^6
T1 <- length(rts.ret)/2
T2 <- length(rts.ret) - T1
VaR.np2 <- numeric()
h <- T1

for (i in (T1 + 1):(T1 + T2)) {
  print((i - T1)/T2)
  # Обучающая выборка
  rts.train <- rts.ret[(i - h):(i - 1)]
  # Ядерная оценка
  f.fix <- npudens(tdat = rts.train, edat = x, ckertype = "gaussian", bwtype = "fixed")
  # Поиск квантиля методом Монте-Карло
  VaR.np2[i - T1] <- q.mc(f.fix, alpha, M)
}
# График кривой VaR
plot(rts.test, type = "l", main = "VaR curve")
lines(VaR.np2, col = "blue", lwd = 2)
legend("bottomleft", c("Returns", "VaR, Fixed Kernel, Monte Carlo"), col = c("black", 
  "blue"), lty = c(1, 1), lwd = c(1, 2))

# Тест VaR-кривой на частоту пробоев
print(paste0("Kupiec test, alpha = ", alpha))
kup <- kupiec.test(rts.test, VaR.np2, alpha)
print(paste0("NP, Fixed kernel, ver2: alpha0 = ", round(kup[1], 6), ", p-value  = ", 
  round(kup[2], 6))) 
