# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied econometrics 8, Multivariate extreme value theory
# Орлов Евгений, 19.01.2015

# ПОСМОТРЕТЬ ЗА СОСТОЯНИЕМ ПЕРЕМеннОЙ alpha
# .h -> .train
# Необходим профайлинг медленных команд
# (посмотреть презентации R in Finance 2014)

# Рассчитать оценки риска для портфеля из двух биржевых индексов
# с помощью многомерных версий метода блочных максим
# и метода превышения порогового значения
# построить кривую VaR для портфеля и проверить качество оценок

# Поскольку в теории экстремальных значений моделируются не доходности а
# убытки, то кривая VaR ограничивает их величину сверху.
# Будьте внимательны при проведении теста Купика и следите
# за знаками неравенств при выявлении пробоев.

# ДАННЫЕ ####
start.date <- 20090804
gazp.df <- read.csv(file.path("data", "ГАЗПРОМ ао [Price].txt"), header=TRUE)
gmkn.df <- read.csv(file.path("data", "ГМКНорНик [Price].txt"), header=TRUE)
tail(gazp.df)
# Цены закрытия начиная с 04.08.2009 по 07.11.2014 включительно
gazp.close <- gazp.df[gazp.df$X.DATE. >= start.date, "X.CLOSE."]
gmkn.close <- gmkn.df[gmkn.df$X.DATE. >= start.date, "X.CLOSE."]
# Даты
gazp.dates <- gazp.df[gazp.df$X.DATE. >= start.date, "X.DATE."]
gazp.dates <- as.Date(as.character(gazp.dates), format="%Y%m%d")
gmkn.dates <- gmkn.df[gmkn.df$X.DATE. >= start.date, "X.DATE."]
gmkn.dates <- as.Date(as.character(gmkn.dates), format="%Y%m%d")
date.comp <- gazp.dates == gmkn.dates
if (sum(date.comp) / length(date.comp) == TRUE) {
  print("Даты совпадают!")
} else {
  print("Даты не совпадают!")
}
# Размер выборок
samp.len <- length(gazp.close)
print(paste0("Размер выборки: ", samp.len))
# Векторы доходностей
gazp.ret <- gazp.close[2:samp.len] / gazp.close[1:(samp.len-1)] - 1
gmkn.ret <- gmkn.close[2:samp.len] / gmkn.close[1:(samp.len-1)] - 1

# Графики цен закрытия
plot(gazp.dates, gazp.close,
     type='l', xlab='index', ylab='close', main='Gazprom:Close')
plot(gmkn.dates, gmkn.close,
     type='l', xlab='index', ylab='close', main='NorilskNickel:Close')
# Графики доходностей
plot(gazp.dates[-1], gazp.ret,
     type='l', xlab='index', ylab='return', main='Gazprom:Return')
plot(gmkn.dates[-1], gmkn.ret,
     type='l', xlab='index', ylab='return', main='NorilskNickel:Return')

# Векторы убытков
gazp.loss <- -gazp.ret
gmkn.loss <- -gmkn.ret
losses <- cbind(gazp.loss, gmkn.loss)

# МЕТОД БЛОЧНЫХ МАКСИМ ####

# Расчет максим
n <- 66  # кол-во дней в максиме (~ кол-во рабочих дней в квартале)
m <- 20  # кол-во максим
T <- m * n
Mn <- rep(0, times=m*2)
dim(Mn)  <- c(m, 2)
for	(i	in	1:2) {
  for	(j	in	1:m) {
    Mn[j, i] <- max(losses[((j-1)*n+1):(j*n), i])
  }
}

# Частные распределения на основе GED
library("evd", quietly=TRUE)
fit1 <- fgev(Mn[, 1])
fit2 <- fgev(Mn[, 2])

# Экстремальные копулы
library("copula", quietly=TRUE)
gumb.cop <- gumbelCopula(2)
gal.cop  <- galambosCopula(2)

# Значения частных функкций распределения
cdf1 <- pgev(Mn[, 1], loc=fit1$estimate[1], scale=fit1$estimate[2],
             shape=fit1$estimate[3])
cdf2 <- pgev(Mn[, 2], loc=fit2$estimate[1], scale=fit2$estimate[2],
             shape=fit2$estimate[3])
cdf  <- cbind(cdf1, cdf2)

# Подгонка копулы
gumb.fit <- fitCopula(cdf, copula=gumb.cop)
gal.fit  <- fitCopula(cdf, copula=gal.cop)
gumb.fit@loglik  # копула Гумбеля показала лучший результат
gal.fit@loglik

# Модельные значения максим
N <- 10^5
cdf.sim  <- rCopula(n=N, copula=gumb.fit@copula)
sim1 <- qgev(cdf.sim[, 1], loc=fit1$estimate[1], scale=fit1$estimate[2],
             shape=fit1$estimate[3])
sim2 <- qgev(cdf.sim[, 2], loc=fit2$estimate[1], scale=fit2$estimate[2],
             shape=fit2$estimate[3])

# Модельные убытки портфеля
w <- c(0.5, 0.5)  # Веса активов в портфеле
port.loss <- sort(w[1]*sim1 + w[2]*sim2)

# Расчет мер риска
k <- 4
alpha <- 1 - 1/k  # в среднем происходит раз в год
VaR  <- port.loss[alpha*N]
ES	<- mean(port.loss[(alpha*N+1):N])
VaR
ES

# Построение кривой VaR
T1 <- length(gazp.loss) / 2
T2 <- length(gazp.loss) - T1
VaR.h <- numeric()
m.h  <- m / 2
h <- T1
for  (i in (T1+1):(T1+T2))	{
  # Векторы убытков
  h.gazp <- gazp.loss[(i-h):(i-1)]
  h.gmkn <- gmkn.loss[(i-h):(i-1)]
  losses.h <- cbind(h.gazp, h.gmkn)
  # Расчет максим
  Mn.h <- rep(0, times=m.h*2)
  dim(Mn.h)	<- c(m.h, 2)
  for (k in	1:2) {
    for	(j in 1:m.h) {
      Mn.h[j, k] <- max(losses.h[((j-1)*n+1):(j*n), k])
    }
  }
  # Частные распределения на основе GED
  fit1.h <- fgev(Mn.h[, 1], std.err=FALSE)
  fit2.h <- fgev(Mn.h[, 2], std.err=FALSE)
  # Экстремальные копулы
  gumb.cop.h <- gumbelCopula(2)
  #gal.cop.h	<- galambosCopula(2)
  # Значения частных функкций распределения
  cdf1.h <- pgev(Mn.h[, 1], loc=fit1.h$estimate[1], scale=fit1.h$estimate[2],
                 shape=fit1.h$estimate[3])
  cdf2.h <- pgev(Mn.h[, 2], loc=fit2.h$estimate[1], scale=fit2.h$estimate[2],
                 shape=fit2.h$estimate[3])
  cdf.h	<- cbind(cdf1.h, cdf2.h)
  # Подгонка копулы
  gumb.fit.h <- try(fitCopula(cdf.h, copula=gumb.cop.h))
  if (class(gumb.fit.h) == "try-error") {
    gumb.fit.h <- fitCopula(cdf.h, copula=gumb.cop.h, method="itau")
  }
  #gal.fit.h	<- try(fitCopula(cdf.h, copula=gal.cop.h))
  #if (class(gal.fit.h) == "try-error") {
  #  gal.fit.h <- fitCopula(cdf.h, copula=gal.cop.h, method="itau")
  #}
  cdf.sim.h	<- rCopula(n=N, copula=gumb.fit.h@copula)
  # Модельные значения максим
  sim1.h <- qgev(cdf.sim.h[, 1], loc=fit1.h$estimate[1],
                 scale=fit1.h$estimate[2], shape=fit1.h$estimate[3])
  sim2.h <- qgev(cdf.sim.h[, 2], loc=fit2.h$estimate[1],
                 scale=fit2.h$estimate[2], shape=fit2.h$estimate[3])
  # Модельные убытки портфеля
  loss.h <- sort(w[1]*sim1.h + w[2]*sim2.h)
  # VaR для данного шага
  VaR.h[i-T1] <- loss.h[alpha*N]
}
# График кривой VaR
port.test <- w[1]*gazp.ret[(T1+1):(T1+T2)] + w[2]*gmkn.ret[(T1+1):(T1+T2)]
plot(port.test, type="l")
lines(-VaR.h, col="red")

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
# Тест VaR-кривой на частоту пробоев
print(paste0("Kupiec test, alpha = ", 1-alpha))
kup <- kupiec.test(port.test, -VaR.h, 1-alpha)
print(paste0("GPD: alpha0 = ", round(kup[1], 6),
             ", p-value  = ", round(kup[2], 6)))


# МЕТОД ПРЕВЫШЕНИЯ ПОРОГА ####

# Выборка значений, превышающих многомерный порог
thresh <- 0.9
u <- c(sort(gazp.loss)[thresh*T], sort(gmkn.loss)[thresh*T])
t.losses <- losses[(losses[, 1] > u[1]) & (losses[, 2] > u[2]), ]
# Частные распределения на основе GED
fit1 <- fpot(t.losses[,1], threshold=u[1], model="gpd", method="SANN")
fit2 <- fpot(t.losses[,2], threshold=u[2], model="gpd", method="SANN")
# Значения частных функций распределения
cdf1 <- pgpd(t.losses[, 1], loc=u[1], scale=fit1$par[1], shape=fit1$par[2])
cdf2 <- pgpd(t.losses[, 2], loc=u[2], scale=fit2$par[1], shape=fit2$par[2])
cdf	<- cbind(cdf1,cdf2)
# Подгонка копулы
gumb.fit <- fitCopula(cdf, copula=gumb.cop)
gal.fit	<- fitCopula(cdf, copula=gal.cop)
gumb.fit@loglik
gal.fit@loglik
# Модельные значения убытков
cdf.sim	<- rCopula(n=N, copula=gal.fit@copula)
sim1 <- qgpd(cdf.sim[, 1], loc=u[1], scale=fit1$par[1], shape=fit1$par[2])
sim2 <- qgpd(cdf.sim[, 2], loc=u[2], scale=fit2$par[1], shape=fit2$par[2])
# Убытки по портфелю
loss <- sort(w[1]*sim1 + w[2]*sim2)
# Расчет мер риска
Fu <- nrow(t.losses) / T
alpha <- 1-1/(260*Fu)
VaR	<- loss[alpha*N]
ES	<- mean(loss[(alpha*N + 1):N])
VaR
ES

# Построение кривой VaR
VaR.gpd <- numeric()
N <- 10^4
h <- T1
for	(i in (T1+1):(T1+T2)) {
  print((i - T1) / T2)
  # Выборка значений, превышающих многомерный порог
  h.gazp <- gazp.loss[(i-h):(i-1)]
  h.gmkn	<- gmkn.loss[(i-h):(i-1)]
  losses.h <- cbind(h.gazp, h.gmkn)
  h.u <- c(sort(h.gazp)[thresh*h], sort(h.gmkn)[thresh*h])
  h.t.losses <- losses.h[(losses.h[, 1]> h.u[1]) & (losses.h[, 2] > h.u[2]), ]
  # Частные распределения на основе GED
  fit1.h <- fpot(h.t.losses[, 1], threshold=h.u[1],
                 model="gpd", method="SANN", std.err=FALSE)
  fit2.h <- fpot(h.t.losses[, 2], threshold=h.u[2],
                 model="gpd", method="SANN", std.err=FALSE)
  # Значения частных функций распределения
  cdf1.h <- pgpd(h.t.losses[, 1], loc=h.u[1], scale=fit1.h$par[1],
                 shape=fit1.h$par[2])
  cdf2.h <- pgpd(h.t.losses[, 2], loc=h.u[2], scale=fit2.h$par[1],
                 shape=fit2.h$par[2])
  cdf.h	<- cbind(cdf1.h, cdf2.h)
  # Подгонка копулы
  gumb.fit.h <- fitCopula(cdf.h, copula=gumb.cop)
  gal.fit.h	<- fitCopula(cdf.h, copula=gal.cop)
  # Модельные значения убытков
  cdf.sim.h	<- rCopula(n=N, copula=gal.fit.h@copula)
  sim1.h <- qgpd(cdf.sim.h[, 1], loc=h.u[1], scale=fit1.h$par[1],
                 shape=fit1.h$par[2])
  sim2.h <- qgpd(cdf.sim.h[, 2], loc=h.u[2], scale=fit2.h$par[1],
                 shape=fit2.h$par[2])
  # Убытки по портфелю
  loss.h <- sort(w[1]*sim1.h + w[2]*sim2.h)
  # Расчет VaR
  Fu <- nrow(h.t.losses) / T
  alpha	<- 1 - 1/(260*Fu)
  VaR.gpd[i-T1] <- loss.h[alpha*N]
}
# График кривой VaR
port.test <- w[1]*gazp.ret[(T1+1):(T1+T2)] + w[2]*gmkn.ret[(T1+1):(T1+T2)]
plot(port.test, type="l")
lines(-VaR.gpd, col="red")
# Тест VaR-кривой на частоту пробоев
print(paste0("Kupiec test, alpha = ", 1-alpha))
kup <- kupiec.test(port.test, -VaR.gpd, 1-alpha)
print(paste0("GPD: alpha0 = ", round(kup[1], 6),
             ", p-value  = ", round(kup[2], 6)))
