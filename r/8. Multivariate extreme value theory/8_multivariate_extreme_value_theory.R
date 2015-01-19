# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied econometrics 8, Multivariate extreme value theory
# Орлов Евгений, 19.01.2015

# ДАННЫЕ ####
gazp.df <- read.csv(file.path("data", "ГАЗПРОМ ао [Price].txt"), header=TRUE)
gmkn.df <- read.csv(file.path("data", "ГМКНорНик [Price].txt"), header=TRUE)
tail(gazp.df)
# Цены закрытия начиная с 04.08.2009 по 07.11.2014 включительно
gazp.close <- gazp.df[gazp.df$X.DATE. >= 20090804, "X.CLOSE."]
gmkn.close <- gmkn.df[gmkn.df$X.DATE. >= 20090804, "X.CLOSE."]
# Даты
gazp.dates <- gazp.df[gazp.df$X.DATE. >= 20090804, "X.DATE."]
gazp.dates <- as.Date(as.character(gazp.dates), format="%Y%m%d")
gmkn.dates <- gmkn.df[gmkn.df$X.DATE. >= 20090804, "X.DATE."]
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
sim1 <- qgev(cdf.sim[,1], loc=fit1$estimate[1], scale=fit1$estimate[2], 
             shape=fit1$estimate[3])
sim2 <- qgev(cdf.sim[,2], loc=fit2$estimate[1], scale=fit2$estimate[2], 
             shape=fit2$estimate[3])

# Модельные убытки портфеля
w <- c(0.5, 0.5)  # Веса активов в портфеле
port.loss <- sort(w[1]*sim1 + w[2]*sim2)

# Расчет мер риска
k <- 4
alpha <- 1 - 1/k
VaR  <- port.loss[alpha*N]
ES	<- mean(port.loss[(alpha*N+1):N])
VaR
ES

# Построение кривой VaR
T1 <- length(gazp.loss) / 2
T2 <- length(gazp.loss) - T1
VaR.h <- numeric()
m.h  <- 10
h <- T1
for  (i in (T1+1):(T1+T2))	{
  # Векторы убытков
  h.gazp <- gazp.loss[(i-h):(i-1)]
  h.gmkn	<- gmkn.loss[(i-h):(i-1)]
  losses.h	<- cbind(h.gazp, h.gmkn)
  # Расчет максим
  Mn.h	<- rep(0, times=m.h*2)
  dim(Mn.h)	<- c(m.h, 2)
  for (k in	1:2) {
    for	(j in 1:m.h) {
      Mn.h[j, k] <- max(losses.h[((j-1)*n+1):(j*n), k])
    }
  }
  # Частные распределения на основе GED
  fit1.h  <- fgev(Mn.h[, 1],	std.err=FALSE)
  fit2.h  <- fgev(Mn.h[, 2],	std.err=FALSE)
  # Экстремальные копулы
  gumb.cop.h <- gumbelCopula(2)
  gal.cop.h	<- galambosCopula(2)
  # Значения частных функкций распределения
  cdf1.h <- pgev(Mn.h[, 1], loc=fit1$estimate[1], scale=fit1$estimate[2], 
                 shape=fit1$estimate[3])
  cdf2.h <- pgev(Mn.h[, 2], loc=fit2$estimate[1], scale=fit2$estimate[2], 
                 shape=fit2$estimate[3])
  cdf.h	<- cbind(cdf1.h, cdf2.h)
  # Подгонка копулы
  gumb.fit.h <- fitCopula(cdf, copula=gumb.cop)
  gal.fit.h	<- fitCopula(cdf, copula=gal.cop)
  cdf.sim.h	<- rCopula(n=N, copula=gumb.fit.h@copula)
  # Модельные значения максим
  sim1.h <- qgev(cdf.sim.h[,1], loc=fit1$estimate[1], scale=fit1$estimate[2], shape=fit1$estimate[3])
  sim2.h <- qgev(cdf.sim.h[,2], loc=fit2$estimate[1], scale=fit2$estimate[2], shape=fit2$estimate[3])
  # Модельные убытки портфеля
  loss.h <- sort(w[1]*sim1.h + w[2]*sim2.h)
  # VaR для данного шага
  VaR.h[i-T1] <- loss.h[alpha*N]
}
# График кривой VaR
port <- w[1]*gazp.ret[(T1+1):(T1+T2)] + w[2]*gmkn.ret[(T1+1):(T1+T2)]
plot(port, type="l")
lines(-VaR.h, col="red")

# Тест Купика
K <- sum(fact>VaR.h)
alpha0 <- K/T2
S <- -2*log((1-alpha)^(T2-K)*alpha^K)+2*log((1-alpha0)^(T2-K)*alpha0^K)
p.value	<- 1-pchisq(S, df=1)
p.value
u <- c(sort(ixic)[0.9*T],sort(nya)[0.9*T])	
t.ESM <- ESM[(ESM[,1]>u[1])&(ESM[,2]>u[2]), ]
fit1 <- fpot(t.ESM[,1], threshold=u[1], model="gpd", method="SANN")
fit2 <- fpot(t.ESM[,2], threshold=u[2], model="gpd", method="SANN")
cdf1 <- pgpd(t.ESM[, 1], loc=u[1], scale=fit1$par[1], shape=fit1$par[2])
cdf2 <- pgpd(t.ESM[, 2], loc=u[2], scale=fit2$par[1], shape=fit2$par[2])
cdf	<- cbind(cdf1,cdf2)
gumb.fit <- fitCopula(cdf, copula=gumb.cop)
gal.fit	<- fitCopula(cdf, copula=gal.cop)
cdf.sim	<- rCopula(n=N, copula=gal.fit@copula)
sim1 <- qgpd(cdf.sim[, 1], loc=u[1], scale=fit1$par[1], shape=fit1$par[2])
sim2 <- qgpd(cdf.sim[, 2], loc=u[2], scale=fit2$par[1], shape=fit2$par[2])
loss <- sort(w[1]*sim1 + w[2]*sim2)
Fu <- nrow(t.ESM) / T
alpha <- 1-1/(260*Fu)
VaR	<- loss[alpha*N]
ES	<- mean(loss[(alpha*N+1):N])
# 4
# VaR Curve
VaR.h <- numeric()
m.h	<- 10
h <- T1
for	(i in (T1+1):(T1+T2)) {
  h.ixic <- ixic[(i-h):(i-1)]
  h.nya	<- nya[(i-h):(i-1)]
  h.u <- c(sort(h.ixic)[0.9*T], sort(h.nya)[0.9*T])	
  h.t.ESM <- ESM[(ESM[,1]>h.u[1])&(ESM[,2]>h.u[2]),]
  fit1.h <- fpot(h.t.ESM[, 1], threshold=h.u[1], model="gpd", method="SANN")
  fit2.h <- fpot(h.t.ESM[, 2], threshold=h.u[2], model="gpd", method="SANN")
  cdf1.h <- pgpd(h.t.ESM[, 1], loc=h.u[1], scale=fit1$par[1], shape=fit1$par[2])
  cdf2.h <- pgpd(h.t.ESM[, 2], loc=h.u[2], scale=fit2$par[1], shape=fit2$par[2])
  cdf.h	<- cbind(cdf1.h,cdf2.h)
  gumb.fit.h <- fitCopula(cdf.h, copula=gumb.cop)
  gal.fit.h	<- fitCopula(cdf.h, copula=gal.cop)
  cdf.sim.h	<- rCopula(n=N, copula=gumb.fit.h@copula)
  sim1.h <- qgpd(cdf.sim.h[, 1], loc=h.u[1], scale=fit1$par[1], shape=fit1$par[2])
  sim2.h <- qgpd(cdf.sim.h[, 2], loc=h.u[2], scale=fit2$par[1], shape=fit2$par[2])
  loss.h <- sort(w[1]*sim1.h + w[2]*sim2.h)
  Fu <- nrow(h.t.ESM) / T
  alpha	<- 1 - 1/(260*Fu)
  VaR.h[i-T] <- loss.h[alpha*N]
}
fact <- ixic[(T1+1):(T1+T2)]
plot(fact, type="l")
lines(VaR.h, col="red")
#Kupiec	Test
K <- sum(fact>VaR.h)
alpha0 <- K/T2
S <- -2*log((1-alpha)^(T2-K)*alpha^K)+2*log((1-alpha0)^(T2-K)*alpha0^K)
p.value	<- 1 - pchisq(S, df=1)
p.value

# 1
#ixic.data <- read.csv("C:/Users/Ruslan/Downloads/IXIC07.csv", header=TRUE, sep=",")
#nya.data <- read.csv("C:/Users/Ruslan/Downloads/NYA07.csv", header=TRUE, sep=",")
#p.ixic <- ixic.data[, 7]
#p.nya <- nya.data[, 7]
#ixic <- -diff(p.ixic) / tail(p.ixic, -1)
#nya  <- -diff(p.nya) / tail(p.nya, -1)
#ixic <- -ixic
#nya <- -nya
#w <- c(0.5, 0.5)
#r <- w[1]*ixic + w[2]*nya
#ESM <- cbind(ixic,nya)
#n <- 90
#m <- 20
#T <- m * n
#Mn <- rep(0, times=m*2)
#dim(Mn)	<- c(m,2)
#for	(i	in	1:2) {
#  for	(j	in	1:m) {
#    Mn[j, i] <- max(ESM[((j-1)*n+1):(j*n), i])
#  }
#}
#library("evd")
#fit1 <- fgev(Mn[, 1])
#fit2 <- fgev(Mn[, 2])
#library("copula")
#gumb.cop <- gumbelCopula(2)
#gal.cop	<- galambosCopula(2)
#cdf1 <- pgev(Mn[, 1], loc=fit1$estimate[1], scale=fit1$estimate[2], shape=fit1$estimate[3])
#cdf2 <- pgev(Mn[, 2], loc=fit2$estimate[1], scale=fit2$estimate[2], shape=fit2$estimate[3])
#cdf	<- cbind(cdf1, cdf2)
#gumb.fit <- fitCopula(cdf, copula=gumb.cop)
#gal.fit	<- fitCopula(cdf, copula=gal.cop)
#gumb.fit@loglik
#gal.fit@loglik
#N <- 10^5
# 2
#cdf.sim	<- rCopula(n=N,copula=gumb.fit@copula)
#sim1 <- qgev(cdf.sim[,1], loc=fit1$estimate[1], scale=fit1$estimate[2], shape=fit1$estimate[3])
#sim2 <- qgev(cdf.sim[,2], loc=fit2$estimate[1], scale=fit2$estimate[2], shape=fit2$estimate[3])
#w <- c(0.5, 0.5)
#loss <- sort(w[1]*sim1 + w[2]*sim2)
#k <- 4
#alpha <- 1-1/k
#VaR	<- loss[alpha*N]
#ES	<- mean(loss[(alpha*N+1):N])
#VaR
#ES

# VaR Curve
#T1 <- (length(ixic)-1) / 2
#T2 <- length(ixic) - T1
#VaR.h <- numeric()
#m.h  <- 10
#h <- T1
#for	(i in (T1+1):(T1+T2))	{
#  h.ixic <- ixic[(i-h):(i-1)]
#  h.nya	<- nya[(i-h):(i-1)]
#  ESM.h	<- cbind(h.ixic, h.nya)
#  Mn.h	<- rep(0, times=m.h*2)
#  dim(Mn.h)	<- c(m.h, 2)
#  for (i in	1:2) {
#	for	(j in 1:m.h) {
#	  Mn.h[j, i] <- max(ESM.h[((j-1)*n+1):(j*n), i])
#    }
#  }
#  fit1.h  <- fgev(Mn.h[,1],	std.err=FALSE)
#  fit2.h  <- fgev(Mn.h[,2],	std.err=FALSE)
#  gumb.cop.h <- gumbelCopula(2)
#  gal.cop.h	<- galambosCopula(2)
#  cdf1.h <- pgev(Mn.h[,1], loc=fit1$estimate[1], scale=fit1$estimate[2], shape=fit1$estimate[3])
#  cdf2.h <- pgev(Mn.h[,2], loc=fit2$estimate[1], scale=fit2$estimate[2], shape=fit2$estimate[3])
#  cdf.h	<- cbind(cdf1.h, cdf2.h)
#  gumb.fit.h <- fitCopula(cdf, copula=gumb.cop)
#  gal.fit.h	<- fitCopula(cdf, copula=gal.cop)
#  cdf.sim.h	<- rCopula(n=N, copula=gumb.fit.h@copula)
# 3
#  sim1.h <- qgev(cdf.sim.h[,1], loc=fit1$estimate[1], scale=fit1$estimate[2], shape=fit1$estimate[3])
#  sim2.h <- qgev(cdf.sim.h[,2], loc=fit2$estimate[1], scale=fit2$estimate[2], shape=fit2$estimate[3])
#  loss.h <- sort(w[1]*sim1.h + w[2]*sim2.h)
#  VaR.h[i-T1] <- loss.h[alpha*N]
#}
# Plotting
#fact <- ixic[(T1+1):(T1+T2)]
#plot(fact, type="l")
#lines(VaR.h, col="red")

