# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied econometrics 3, Generalized hyperbolic distribution
# Орлов Евгений, 20.11.2014

# ДАННЫЕ ####
rts.df <- read.csv("RTSI [Price].txt", header=TRUE)
tail(rts.df)
# Цены закрытия начиная с 01.01.2010 по 07.11.2014 включительно
rts.close <- rts.df[rts.df$X.DATE. >= 20100101, "X.CLOSE."]
# Вектор доходностей
rcl <- length(rts.close)
print(paste0("Размер выборки: ", rcl))
rts.ret <- rts.close[2:rcl] / rts.close[1:(rcl-1)] - 1

rts.dates <- rts.df[rts.df$X.DATE. >= 20100101, "X.DATE."]
rts.dates <- as.Date(as.character(rts.dates), format="%Y%m%d")
# График цен закрытия
plot(rts.dates, rts.close, 
     type='l', xlab='index', ylab='close', main='RTS Index:Close')
# График доходностей
plot(rts.dates[-1], rts.ret, 
     type='l', xlab='index', ylab='return', main='RST Index:Return')

# ОБОБЩЕННОЕ ГИПЕРБОЛИЧЕСКОЕ РАСПРЕДЕЛЕНИЕ ####
library(ghyp)

# Подгон параметров распределения
# Эмпирическая плотность
rts.emp_dens<-density(rts.ret, bw='nrd')
# Нормальное распределение
rts.gauss <- fit.gaussuv(rts.ret)
# t-распределение Стьюдента
rts.t <- fit.tuv(rts.ret, silent=TRUE)
# Нормально-обратное гауссовское
rts.nig <- fit.NIGuv(rts.ret, silent=TRUE)
# Variance-Gamma
rts.vg <- fit.VGuv(rts.ret, silent=TRUE)
# Гиперболическое
rts.hyp <- fit.hypuv(rts.ret, silent=TRUE)
# Обобщенное гиперболическое
rts.ghyp <- fit.ghypuv(rts.ret, silent=TRUE)

# Визуализация подобранных распределений
colors <- palette(rainbow(6))
plot(rts.emp_dens, type="l", lwd=3, 
     main="Fitted density functions", xlab='Return')
lines(rts.gauss, type="l", col=colors[1])
lines(rts.t, type="l", col=colors[2])
lines(rts.nig, type="l", col=colors[3])
lines(rts.vg, type="l", col=colors[4])
lines(rts.hyp, type="l", col=colors[5])
lines(rts.ghyp, type="l", col=colors[6])
legend('topleft', c("emp_dens", "gauss", "t", "nig", "vg", "hyp", "ghyp"), 
       col=c("black", colors), lty=rep(1, 7))

# Выбор наилучшей модели на основе критерия Акаике
aic.uv <- stepAIC.ghyp(rts.ret, dist=c("gauss","t","hyp","ghyp","nig", "vg"),
                       silent=TRUE)
rts.best <- aic.uv$best.model
summary(rts.best)  # Информация о наилучшей модели
aic.uv$fit.table  # Сводная информация об оцененных моделях

# Сравнение наилучшей модели с нормальным распределением

# Гистограмма и график квантиль-квантиль 
hist(rts.best)
qqghyp(rts.best)

# Тест Колмогорова-Смирнова для выбранной модели и нормального распределения
ks.test(rts.ret, rghyp(n=100*length(rts.ret),object=rts.best))
ks.test(rts.ret, rghyp(n=100*length(rts.ret),object=rts.gauss))

# ОЦЕНКИ РИСКА МЕТОДОМ МОНТЕ-КАРЛО НА ОСНОВЕ НАИЛУЧШЕЙ МОДЕЛИ ####

# Параметры расчета VaR и ES
conf.level <- 0.90  # доверительный уровень VaR
alpha <- 1 - conf.level  # уровень квантиля
N <- 10^6  # кол-во сгенерированных доходностей

# Симуляция и вывод результатов
rts.best.sim <- sort(rghyp(n=N, object=rts.best))
rts.gauss.sim <- sort(rghyp(n=N, object=rts.gauss))
VaR1.best <- rts.best.sim[alpha*N]
VaR1.gauss <- rts.gauss.sim[alpha*N]
VaR2.best <- qghyp(alpha, object=rts.best)
VaR2.gauss <- qghyp(alpha, object=rts.gauss)
ES.best <- mean(rts.best.sim[1:(alpha*N-1)])
ES.gauss <- mean(rts.gauss.sim[1:(alpha*N-1)])
print(paste0("VaR1 best model/gauss: ", round(VaR1.best, 4), 
             "/", round(VaR1.gauss, 4)))
print(paste0("VaR2 best model/gauss: ", round(VaR2.best, 4), 
             "/", round(VaR2.gauss, 4)))
print(paste0("ES best model/gauss: ", round(ES.best, 4), 
             "/", round(ES.gauss, 4)))

# VaR-КРИВЫЕ ####
T1 <- 500
T2 <- length(rts.ret)-T1
VaR.hyp <- numeric()  # VaR-кривая на основе гиперболического распределения
VaR.gauss <- numeric()  # VaR-кривая на основе нормального распределения

# Построение VaR-кривых
h <- T1  # длина обучающей выборки
for (i in (T1+1):(T1+T2)) {
  rts.train <- rts.ret[(i-h):(i-1)]
  rts.hyp.fit <- fit.hypuv(rts.train, symmetric=TRUE, silent=TRUE)
  rts.gauss.fit <- fit.gaussuv(rts.train)
  
  VaR.hyp[i-T1] <- qghyp(alpha, object=rts.hyp.fit)
  VaR.gauss[i-T1] <- qghyp(alpha, object=rts.gauss.fit)
}
# График
rts.test <- rts.ret[(T1+1):(T1+T2)]
plot(rts.test, type="l", main="VaR curves")
lines(VaR.hyp, col="blue", lwd=2)
lines(VaR.gauss, col="red", lwd=2)
legend('bottomleft', c("Returns", "VaR hyp", "VaR gauss"), 
       col=c("black", "blue", "red"), lty=rep(1, 3))

# СРАВНЕНИЕ МОДЕЛЕЙ И ВЫВОДЫ ####

# Частота пробоев
kupiec.test <- function(ret, VaR, alpha) {
  # Тест Купика:
  # H0: модельная и эмпирическая частоты пробоя VaR совпадают
  K <- sum(ret < VaR)
  T2 <- length(ret)
  alpha0 <- K / T2
  S <- -2*log((1-alpha)^(T2-K) * alpha^K) + 2*log((1-alpha0)^(T2-K) * alpha0^K)
  p.value <- 1-pchisq(S, df=1)
  return(c(alpha0, p.value))
}

# Глубина пробоев
lopez.lf <- function(ret, VaR) {
  # Функция потерь Лопеса
  K <- sum(ret < VaR)
  value <- sum((ret-VaR)^2*(ret < VaR)) / K
}

blanco.lf <- function(ret, VaR) {
  # Функция потерь Бланко-Ила
  K <- sum(ret < VaR)
  value <- sum((ret-VaR)/VaR * (ret < VaR)) / K
}

# Тест VaR-кривых на частоту пробоев
print(paste0("Kupiec test, alpha = ", alpha))
kup <- kupiec.test(rts.test, VaR.hyp, alpha)
print(paste0("hyp: alpha0 = ", round(kup[1], 6), 
             ", p-value  = ", round(kup[2], 6)))
kup <- kupiec.test(rts.test, VaR.gauss, alpha)
print(paste0("gauss: alpha0 = ", round(kup[1], 6), 
             ", p-value  = ", round(kup[2], 6)))
# Сравнение глубины пробоев
print(paste0("Lopez loss function hyp/gauss: ", 
             round(lopez.lf(rts.test, VaR.hyp), 6), 
             "/", round(lopez.lf(rts.test, VaR.gauss), 6)))
print(paste0("Blanco loss function hyp/gauss: ", 
             round(blanco.lf(rts.test, VaR.hyp), 6), 
             "/", round(blanco.lf(rts.test, VaR.gauss), 6)))
