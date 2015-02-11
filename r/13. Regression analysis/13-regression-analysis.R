# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied econometrics 13, Regression Analysis
# Орлов Евгений, 29.12.2014
#

# ПОСТРОИТЬ ПРОГНОЗ В ПОСЛЕДНЕМ СЛУЧАЕ
# (ГОМОСКЕДАСТИЧНОСТЬ, ГЕТЕРСКЕДАСТИЧНОСТЬ)

# Линейная регрессия ====

# исходные данные
library(datasets)
ozone <- airquality$Ozone
rad <- airquality$Solar.R

rem <- is.na(ozone) | is.na(rad)
ozone <- ozone[!rem]
rad <- rad[!rem]

# разделим выборку на обучающую и экзаменующую
N <- length(ozone)
E <- 20
T <- N - E
train.obs <- (1:T)
eval.obs <- ((T + 1):N)

t.rad <- rad[train.obs]
t.ozone <- ozone[train.obs]
e.rad <- rad[eval.obs]
e.ozone <- ozone[eval.obs]

# регрессионная модель
fit.par <- lm(ozone ~ radiation,
              data = data.frame(radiation = t.rad, ozone = t.ozone),
              weights = NULL)

# другой вариант
fit.par2 <- lm(t.ozone ~ t.rad)

# анализ качества модели
summary(fit.par)
fit.par$coefficients
fit.par$residuals
fit.par$fitted.values

plot(t.rad, t.ozone, pch = 16, xlab = "radiation", ylab = "ozone")
z <- order(t.rad)
lines(t.rad[z], fit.par$fitted.values[z], col = "blue", lwd = 3)

# анализ остатков модели
res <- fit.par$residuals
hist(res)
plot(res, type = "l")
# тесты на нормальность
library("fBasics")
shapiro.test(res)
jarqueberaTest(res)

# переформулировка модели
fit.par <- lm(log(ozone) ~ rad + rad2,
              data = data.frame(rad = t.rad, rad2 = t.rad^2, ozone = t.ozone))
summary(fit.par)

plot(t.rad, t.ozone, pch = 16, xlab= "radiation", ylab = "ozone")
z <- order(t.rad)
lines(t.rad[z], exp(fit.par$fitted.values[z]), col = "blue", lwd = 3)

# анализ остатков модели
res <- fit.par$residuals
hist(res)
plot(res, type = "l")
# тесты на нормальность
shapiro.test(res)
jarqueberaTest(res)

# Тесты на гетероскедастичность
library("lmtest")
# тест Бреуша-Паган
bptest(fit.par, varformula = NULL, data = NULL, studentize = FALSE)
# тест Голфелда-Куандта
gqtest(fit.par, fraction = 25, alternative = "two.sided")

# Учёт гетероскедастичности
# Регрессия с гетероскедастичностью
oz <- lm(t.ozone ~ t.rad)
bptest(oz, varformula = NULL, data = NULL, studentize = FALSE)
# Оценки стандартных отклонений
e.sq <- oz$residuals^2
sigma.hat <- lm(e.sq ~ t.rad)$fitted.values^0.5
# Используем взвешенный МНК
oz.wgt <- lm(t.ozone ~ t.rad, weights = 1 / sigma.hat)
summary(oz)
summary(oz.wgt)

# Тест на причинность
# тест Гранжера
grangertest(ozone ~ rad, order = 2)

# Тест на автокорреляцию
# тест Дарбина-Ватсона
dwtest(fit.par, alternative = "two.sided")

# Переформулировка модели
ar1 <- log(t.ozone[1:(T - 1)])
fit.par <- lm(log(ozone) ~ rad + rad2 + ar1,
              data = data.frame(ozone = t.ozone[2:T], rad = t.rad[2:T],
                                rad2 = t.rad[2:T]^2))
summary(fit.par)
tr <- t.rad[2:T]
tz <- t.ozone[2:T]
z <- order(tr)
plot(tr[z], tz[z], type = "l", lty = "dashed",
     xlab= "radiation", ylab = "ozone")
lines(tr[z], exp(fit.par$fitted.values[z]), col = "blue", lwd = 3)
# Анализ остатков модели
res <- fit.par$residuals
# нормальность
shapiro.test(res)
jarqueberaTest(res)
# гетероскедастичность
bptest(fit.par, varformula = NULL, data = NULL, studentize = FALSE)
gqtest(fit.par, fraction = 25, alternative = "two.sided")
# автокорреляция
dwtest(fit.par, alternative = "two.sided")

# Построение прогноза
frc.par <- predict(fit.par,
                   newdata = data.frame(rad = e.rad[2:E], rad2 = e.rad[2:E]^2,
                                        ar1 = log(e.ozone[1:(E - 1)])),
                   se.fit = TRUE,
                   interval = "prediction",
                   level = 0.90)
frc.par$fit  # прогнозные значения и интервалы
frc.par$se.fit  # стандартные ошибки

er <- e.rad[2:E]
ez <- e.ozone[2:E]
z <- order(er)
plot(er[z], ez[z], pch = 16, ylim = c(0, 150),
     xlab= "radiation", ylab = "ozone")
lines(er[z], exp(frc.par$fit[z, 1]), col = "blue", lwd = 2)
lines(er[z], exp(frc.par$fit[z, 2]), col = "red", lwd = 2, lty = "dashed")
lines(er[z], exp(frc.par$fit[z, 3]), col = "red", lwd = 2, lty = "dashed")


# Непараметрическая регрессия - одномерный случай ====

# Ядерная оценка Надарая-Уотсона

# расчет величины h
library("np")
bw <- npregbw(ozone ~ rad, ckertype = "gaussian", bwtype = "fixed",
              data = data.frame(ozone = t.ozone, rad = t.rad))
h <- bw$bw

# ядро и функция Надарая-Уотсона
kern <- function(x) {
  exp(-(x^2 / 2)) / sqrt(2 * pi)
}

NW <- function(x, x.dat, y.dat, h) {
  K1 <- 0
  K2 <- 0
  N <- length(y.dat)
  for (i in 1:N) {
    K1 <- K1 + kern((x - x.dat[i]) / h) * y.dat[i]
    K2 <- K2 + kern((x - x.dat[i]) / h)
  }
  K1 / K2
}

# График оценки
plot(t.rad, t.ozone, pch = 16)
z <- order(t.rad)
lines(t.rad[z], NW(t.rad, t.rad, t.ozone, h)[z], col = "blue", lwd = 3)

# График оценки (другой вариант)
plot(t.rad, t.ozone, pch = 16)
fit.npar <- npreg(bw)
ozone.hat <- predict(fit.npar)
lines(t.rad[z], ozone.hat[z], col = "blue", lwd = 3)

# Построение прогноза, sigma^2 = const
e <- t.ozone - ozone.hat
s2 <- var(e)
g <- (4 / (3 * T))^(1/5) * sqrt(s2)
# симулированные значения ошибок (метод бутстрапа)
b <- 10^4
e.star <- e[sample(1:T, size = b, replace = TRUE)] + g * rnorm(b)
e.star <- sort(e.star)
# прогноз и доверительные границы
alpha <- 0.1
y <- predict(fit.npar, newdata = data.frame(rad = e.rad))
bottom <- y + e.star[alpha / 2 * b]
top <- y + e.star[(1 - alpha / 2) * b]
# график прогноза и доверительных границ
z <- order(e.rad)
plot(e.rad, e.ozone, ylim = range(c(top, bottom)))
lines(e.rad[z], y[z], col = "blue", lwd = 2)
lines(e.rad[z], top[z], col = "red", lty = "dashed")
lines(e.rad[z], bottom[z], col = "red", lty = "dashed")

# Построение прогноза, sigma^2 = sigma^2(z)
# моделирование условной дисперсии
s2 <- (e - mean(e))^2  # оценка условной дисперсии остатков
h.res <- npregbw(s2 ~ rad, ckertype = "gaussian", bwtype = "fixed",
                 data = data.frame(rad = t.rad))
res.npar <- npreg(h.res)  # модель оценки условной дисперсии ошибок

# метод бутстрапа
b <- 10^4
alpha <- 0.1
top <- numeric(E)
bottom <- numeric(E)
y <- predict(fit.npar, newdata = data.frame(rad = e.rad))

for (i in 1:E) {
  s2.hat <- predict(res.npar, newdata = data.frame(rad = e.rad[i]))
  g <- (4 / (3 * T))^(1/5) * sqrt(s2.hat)

  e.star <- e[sample(1:T, size = b, replace = TRUE)] + g * rnorm(b)
  e.star <- sort(e.star)

  bottom[i] <- y[i] + e.star[alpha / 2 * b]
  top[i] <- y[i] + e.star[(1 - alpha / 2) * b]
}
# график прогноза и доверительных границ
z <- order(e.rad)
plot(e.rad, e.ozone, ylim = range(c(top, bottom)))
lines(e.rad[z], y[z], col = "blue", lwd = 2)
lines(e.rad[z], top[z], col = "red", lty = "dashed")
lines(e.rad[z], bottom[z], col = "red", lty = "dashed")


# Непараметрическая регрессия - двумерный случай ====

# ядерная оценка Надарая-Уотсона
kern2 <- function(x, dim = 2) {
  exp(-sum(x^2) / 2) * (2 * pi)^(-dim / 2)
}

NW2 <- function(x, x.dat, y.dat, h, dim = 2) {
  N <- length(y.dat)
  M <- nrow(x)
  K1 <- numeric(M)
  K2 <- numeric(M)
  for (j in 1:M) {
    for (i in 1:M) {
      K1[j] <- K1[j] + kern2((x[j, ] - x.dat[i, ]) / h, dim) * y.dat[i]
      K2[j] <- K2[j] + kern2((x[j, ] - x.dat[i, ]) / h, dim)
    }
  }
  K1 / K2
}

# добавочная объясняющая переменная
temp <- airquality$Temp
temp <- temp[!rem]
t.temp <- temp[train.obs]
e.temp <- temp[eval.obs]

# регрессионная модель
bw2 <- npregbw(ozone ~ rad + temp, ckertype = "gaussian", bwtype = "fixed",
               data = data.frame(ozone = t.ozone, rad = t.rad, temp = t.temp))
rad.temp <- cbind(t.rad, t.temp)
ozone.hat <- NW2(rad.temp, rad.temp, t.ozone, bw2$bw, 2)

# альтернативный вариант
ozone.reg <- npreg(bw2)
ozone.hat <- predict(ozone.reg)

# график
zt <- order(t.rad)
plot(t.ozone[zt], pch = 16)
lines(ozone.hat[zt], col = "blue", lwd = 2)

# Построение прогноза



