# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied econometrics 4, VaR curve using GARCH Орлов Евгений,
# 20.11.2014

# ДАННЫЕ ####
sber.df <- read.csv(file.path("data", "Сбербанк [Price].txt"), header = TRUE)
tail(sber.df)
# Цены закрытия начиная с 01.01.2010 по 07.11.2014
# включительно
sber.close <- sber.df[sber.df$X.DATE. >= 20100101, "X.CLOSE."]
# Вектор доходностей
scl <- length(sber.close)
print(paste0("Размер выборки: ", scl))
# Вектор доходностей
sber.ret <- sber.close[2:scl]/sber.close[1:(scl - 1)] - 1

sber.dates <- sber.df[sber.df$X.DATE. >= 20100101, "X.DATE."]
sber.dates <- as.Date(as.character(sber.dates), format = "%Y%m%d")
# График цен закрытия
plot(sber.dates, sber.close, type = "l", xlab = "index", ylab = "close", main = "Sberbank:Close")
# График доходностей
plot(sber.dates[-1], sber.ret, type = "l", xlab = "index", ylab = "return", main = "Sberbank:Return")

# ТЕСТЫ НА СТАЦИОНАРНОСТЬ ####
library("tseries", quietly = TRUE)
# ADF-тест
adf.test(sber.ret)
# PP-тест
pp.test(sber.ret)
# KPSS-тест
kpss.test(sber.ret, null = "Level")

# ОПТИМАЛЬНЫЙ ПОРЯДОК МОДЕЛИ ARMA(m, n) НА ОСНОВЕ
# КРИТЕРИЯ АКАИКЕ ####
arma.coef.aic <- function(ret, a, b) {
  # Выбор оптимальных коэффициентов ARMA(a, b) для
  # ret
  coefs <- expand.grid(a = a, b = b)
  aic <- apply(coefs, 1, function(coef) AIC(arima(ret, order = c(coef[1], 0, coef[2]))))
  ind <- which(aic == min(aic))
  return(c(coefs$a[ind], coefs$b[ind]))
}
## Наилучшие параметры m (1:5) и n (1:5) для модели
## ARMA(m,n)
arma.coef.aic(sber.ret, 1:5, 1:5)

# ТЕСТ НА ARCH-ЭФФЕКТЫ #### Тест множителей
# Лагранжа (LM-тест)
library("FinTS", quietly = TRUE)
ArchTest(sber.ret, lags = 12)

# VaR НА ОСНОВЕ МОДЕЛИ ARMA-GARCH ####
library("fGarch", quietly = TRUE)
# Оптимизация параметров p, q GARCH(p, q)
arma_garch.coef.aic <- function(ret, a, b, c, d) {
  # Выбор оптимальных коэффициентов ARMA(a,
  # b)-GARCH(c, d) для ret
  coefs <- expand.grid(a = a, b = b, c = c, d = d)
  aic <- apply(coefs, 1, function(coef) garchFit(substitute(~arma(a, b) + garch(c, 
    d), list(a = coef[1], b = coef[2], c = coef[3], d = coef[4])), data = ret, 
    trace = FALSE)@fit$ics[1])
  ind <- which(aic == min(aic))
  return(c(coefs$a[ind], coefs$b[ind], coefs$c[ind], coefs$d[ind]))
}
arma_garch.coef.aic(sber.ret, 2, 2, 1:5, 1:5)

# Наилучшая модель
sber.gfit <- garchFit(~arma(2, 2) + garch(2, 1), data = sber.ret, trace = FALSE)
plot(sber.gfit, which = 13)
plot(sber.gfit, which = 10)

# Расчет VaR по всей выборке наблюдений
sber.pred <- predict(sber.gfit, n.ahead = 1)
alpha <- 0.05
VaR <- as.numeric(sber.pred[1] + sber.pred[3] * qged(alpha, mean = 0, sd = 1))
VaR

# КРИВАЯ VaR ####
T1 <- 500
T2 <- length(sber.ret) - T1
VaR.arma_garch <- numeric()

h <- T1  # Длина обучающей выборки
for (i in (T1 + 1):(T1 + T2)) {
  sber.train <- sber.ret[(i - h):(i - 1)]
  sber.gfit <- garchFit(formula = ~arma(1, 1) + garch(1, 1), data = sber.train, 
    delta = 2, include.delta = FALSE, include.shape = FALSE, include.skew = FALSE, 
    trace = FALSE)
  sber.pred <- predict(sber.gfit, n.ahead = 1)
  VaR.arma_garch[i - T1] <- as.numeric(sber.pred[1] + sber.pred[3] * qsged(alpha, 
    mean = 0, sd = 1))
}
# График кривой VaR
sber.test <- sber.ret[(T1 + 1):(T1 + T2)]
plot(sber.test, type = "l", main = "VaR curve")
lines(VaR.arma_garch, col = "blue", lwd = 1)
legend("bottomleft", c("Returns", "VaR ARMA-GARCH"), col = c("black", "blue"), lty = c(1, 
  1))

# Верификация VaR Частота пробоев
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

# Глубина пробоев
lopez.lf <- function(ret, VaR) {
  # Функция потерь Лопеса
  K <- sum(ret < VaR)
  value <- sum((ret - VaR)^2 * (ret < VaR))/K
}

blanco.lf <- function(ret, VaR) {
  # Функция потерь Бланко-Ила
  K <- sum(ret < VaR)
  value <- sum((ret - VaR)/VaR * (ret < VaR))/K
}

# Тест VaR-кривой на частоту пробоев
print(paste0("Kupiec test, alpha = ", alpha))
kup <- kupiec.test(sber.test, VaR.arma_garch, alpha)
print(paste0("ARMA(1, 1)-GARCH(1, 1): alpha0 = ", round(kup[1], 6), ", p-value  = ", 
  round(kup[2], 6)))
# Оценка глубины пробоев
print(paste0("Lopez loss function: ", round(lopez.lf(sber.test, VaR.arma_garch), 
  6)))
print(paste0("Blanco loss function: ", round(blanco.lf(sber.test, VaR.arma_garch), 
  6))) 
