# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied econometrics 7, Univariate extreme value theory Орлов
# Евгений, 21.11.2014

# ДАННЫЕ ####
sber.df <- read.csv(file.path("data", "Сбербанк [Price].txt"), header = TRUE)
tail(sber.df)
# Цены закрытия начиная с 01.08.2009 по 07.11.2014
# включительно
start.date <- 20090801
sber.close <- sber.df[sber.df$X.DATE. >= start.date, "X.CLOSE."]

sber.dates <- sber.df[sber.df$X.DATE. >= start.date, "X.DATE."]
sber.dates <- as.Date(as.character(sber.dates), format = "%Y%m%d")
# График цен закрытия
plot(sber.dates, sber.close, type = "l", xlab = "index", ylab = "close", main = "Sberbank:Close")

# МЕТОД БЛОЧНЫХ МАКСИМ ####
n <- 65  # кол-во дней в максиме (~ кол-во рабочих дней в квартале)
m <- 20  # кол-во максим
size <- n * m
sber.loss <- sber.close[2:(size + 1)]/sber.close[1:size] - 1
sber.loss <- -sber.loss  # Вектор убытков
# Расчет максим
Mn <- rep(0, times = m)
for (i in 1:m) {
  Mn[i] <- max(sber.loss[((i - 1) * n + 1):(i * n)])
}
# Распределение максим на основе GEV
library("evd", quietly = TRUE)
Mn.fit <- fgev(Mn)
plot(Mn.fit, which = 2)  # график 'квантиль-квантиль'
plot(Mn.fit, which = 3)  # график эмпирической плотности

mu <- Mn.fit$estimate[1]
sigma <- Mn.fit$estimate[2]
xi <- Mn.fit$estimate[3]
k <- 4
u <- 0.09
r.nk <- mu + sigma/xi * ((-log(1 - 1/k))^(-xi) - 1)
k.nr <- 1/(1 - pgev(u, loc = mu, scale = sigma, shape = xi))
print(paste0("Уровень потерь, который будет пройден в среднем 1 раз в год: ", 
  round(r.nk, 6)))
print(paste0("Средний период наступления убытка > 9% (в кварталах): ", 
  round(k.nr, 2)))

# VaR, ES ПО ВСЕЙ ВЫБОРКЕ С ИСПОЛЬЗОВАНИЕМ
# РАСПРЕДЕЛЕНИЯ ПАРЕТО ####

alpha1 <- 0.95
# Пороговое значение
u <- sort(sber.loss)[(alpha1) * length(sber.loss)]
# Подбор параметров
gpd.fit <- fpot(sber.loss, threshold = u, model = "gpd", method = "SANN")
plot(gpd.fit, which = 2)
plot(gpd.fit, which = 3)
# Оценки параметров
beta <- gpd.fit$estimate[1]
xi <- gpd.fit$estimate[2]
Fu <- gpd.fit$pat
alpha2 <- 1 - 1/260  # соответствует 1 превышению в год
VaR1 <- u + beta/xi * (((1 - alpha1)/Fu)^(-xi) - 1)  # 95%
VaR2 <- u + beta/xi * (((1 - alpha2)/Fu)^(-xi) - 1)  # раз в год
ES1 <- (VaR1 + beta - xi * u)/(1 - xi)  # 95%
ES2 <- (VaR2 + beta - xi * u)/(1 - xi)  # раз в год

# КРИВАЯ VaR С ИСПОЛЬЗОВАНИЕМ РАСПРЕДЕЛЕНИЯ
# ПАРЕТО ####
T1 <- 500
T2 <- length(sber.loss) - T1
alpha0 <- 0.95  # Персентиль порогового значения
alpha <- 0.95  # VaR
VaR.gpd <- numeric()

h <- T1  # Длина обучающей выборки
for (i in (T1 + 1):(T1 + T2)) {
  sber.train <- sber.loss[(i - h):(i - 1)]
  u <- sort(sber.train)[(alpha0) * T1]
  gpd.fit <- fpot(sber.train, threshold = u, model = "gpd", method = "SANN", std.err = FALSE)
  beta <- gpd.fit$estimate[1]
  xi <- gpd.fit$estimate[2]
  Fu <- gpd.fit$pat
  if (xi != 0) {
    VaR.gpd[i - T1] <- -(u + beta/xi * (((1 - alpha)/Fu)^(-xi) - 1))
  } else {
    VaR.gpd[i - T1] <- -(u - beta * log((1 - alpha)/Fu))
  }
}
# График кривой VaR
sber.test <- -sber.loss[(T1 + 1):(T1 + T2)]
plot(sber.test, type = "l", main = "VaR curve")
lines(VaR.gpd, col = "blue", lwd = 2)
legend("bottomleft", c("Returns", "VaR GPD"), col = c("black", "blue"), lty = c(1, 
  1), lwd = c(1, 2))

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
print(paste0("Kupiec test, alpha = ", 1 - alpha))
kup <- kupiec.test(sber.test, VaR.gpd, 1 - alpha)
print(paste0("GPD: alpha0 = ", round(kup[1], 6), ", p-value  = ", round(kup[2], 6))) 
