# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied econometrics 10, nonparametric modeling Part 2. Adaptive and
# multivariate kernel estimators Орлов Евгений, 26.01.2015

# рассчитать оценки риска для портфеля из
# двух биржевых индексов с помощью
# многомерного адаптивного метода (с
# величинами λ_i) построить кривую VaR для
# портфеля и проверить качество оценок

# ДАННЫЕ ####
start.date <- 20090804
gazp.df <- read.csv(file.path("data", "ГАЗПРОМ ао [Price].txt"), header = TRUE)
gmkn.df <- read.csv(file.path("data", "ГМКНорНик [Price].txt"), header = TRUE)
tail(gazp.df)
# Цены закрытия начиная с 04.08.2009 по 07.11.2014
# включительно
gazp.close <- gazp.df[gazp.df$X.DATE. >= start.date, "X.CLOSE."]
gmkn.close <- gmkn.df[gmkn.df$X.DATE. >= start.date, "X.CLOSE."]
# Даты
gazp.dates <- gazp.df[gazp.df$X.DATE. >= start.date, "X.DATE."]
gazp.dates <- as.Date(as.character(gazp.dates), format = "%Y%m%d")
gmkn.dates <- gmkn.df[gmkn.df$X.DATE. >= start.date, "X.DATE."]
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


# Одномерный адаптивный метод ####

hist(gazp.ret, probability = TRUE)
min(gazp.ret)
max(gazp.ret)

L <- 10^3
N <- length(gazp.ret)
x <- seq(-0.15, 0.1, length = L)

library("np", quietly = TRUE)
f.fix <- npudens(tdat = gazp.ret, edat = x, ckertype = "gaussian", bwtype = "fixed")

# Адаптивный метод с λi
pilot <- npudens(tdat = gazp.ret, ckertype = "gaussian", bwtype = "fixed")
h <- pilot$bws$bw  # оценка глобальной составляющей интервала

# среднегеометрическое пилотных оценок
g <- 1
for (i in 1:N) {
  g <- g * pilot$dens[i]^(1/N)
}

# расчёт локальной концентрации наблюдений
alpha <- 0.5
lambda <- (g/pilot$dens)^alpha

kern <- function(u) {
  exp(-u^2/2)/sqrt(2 * pi)  # ядро Гаусса
}

# расчёт оценок плотности
f.ada <- numeric(L)
for (i in 1:L) {
  f.ada[i] <- sum(kern((x[i] - gazp.ret)/(h * lambda))/(h * lambda))
}
f.ada <- f.ada/N

# Сравнение адаптивной и фиксированной
# оценок
plot(x, f.fix$dens, type = "l", lty = "dashed", main = "Fixed and adaptive estimates", 
  xlab = "y", ylab = "Density")
lines(x, f.ada)
rug(gazp.ret)

# Нахождение квантилей оценки распределения

# оценка функции распределения
F.fix <- npudist(tdat = gazp.ret, edat = x, ckertype = "gaussian", bwtype = "fixed")

# для адаптивного варианта
F.ada <- rep(0, times = L)
dx <- x[2] - x[1]
for (i in 1:L) {
  F.ada[i] <- sum(f.ada[1:i]) * dx
}

# поиск квантиля методом деления пополам
alpha <- 0.01
a <- 1
b <- L
ab <- trunc((a + b)/2)
while ((b - a) > 2) {
  if (F.fix$dist[ab] <= alpha) {
    a <- ab
  }
  if (F.fix$dist[ab] >= alpha) {
    b <- ab
  }
  ab <- trunc((a + b)/2)
}
q.fix1 <- x[ab]

# для адаптивного варианта
a <- 1
b <- L
ab <- trunc((a + b)/2)
while ((b - a) > 2) {
  if (F.ada[ab] <= alpha) {
    a <- ab
  }
  if (F.ada[ab] >= alpha) {
    b <- ab
  }
  ab <- trunc((a + b)/2)
}
q.ada1 <- x[ab]

# Генератор случайных чисел (метод
# Монте-Карло)

# фиксированный интервал
M <- 10^6
y.fix.sim <- sample(x, prob = f.fix$dens, size = M, replace = TRUE)
q.fix2 <- sort(y.fix.sim)[alpha * M]

# для адаптивного варианта
y.ada.sim <- sample(x, prob = f.ada, size = M, replace = TRUE)
q.ada2 <- sort(y.ada.sim)[alpha * M]


# Многомерный адаптивный метод ####

port.ret <- cbind(gazp.ret, gmkn.ret)
N <- nrow(port.ret)
min(port.ret[, 1])
max(port.ret[, 1])
min(port.ret[, 2])
max(port.ret[, 2])

# сетка для расчёта оценок плотности
L <- 50
u <- seq(-0.15, 0.15, length = L)
v <- seq(-0.15, 0.15, length = L)
uv <- expand.grid(u, v)

# оценка плотности
f2.fix <- npudens(tdat = port.ret, edat = uv, ckertype = "gaussian", bwtype = "fixed")

# графики оценки
w <- f2.fix$dens
dim(w) <- c(L, L)
persp(u, v, w, theta = 30, main = "Bivariate kernel estimate, 3D plot", xlab = "Gazprom returns", 
  ylab = "NorNickel returns", zlab = "Density")
contour(u, v, w, nlevel = 7, main = "Bivariate kernel estimate, contour plot", xlab = "Gazprom returns", 
  ylab = "NorNickel returns")

# Адаптивный метод с λi, аналогично
# одномерному случаю
pilot2 <- npudens(tdat = port.ret, ckertype = "gaussian", bwtype = "fixed")
h2 <- pilot$bws$bw

g2 <- 1
for (i in 1:N) {
  g2 <- g2 * pilot2$dens[i]^(1/N)
}

alpha2 <- 0.5
lambda2 <- (g2/pilot2$dens)^alpha2

kern2 <- function(x) {
  exp(-(x[1]^2 + x[2]^2)/2)/(2 * pi)  # 2-мерное ядро Гаусса
}

f2.ada <- rep(0, times = L^2)
for (i in 1:(L^2)) {
  for (j in 1:N) {
    f2.ada[i] <- f2.ada[i]
    +kern2((uv[i, ] - port.ret[j, ])/(h2 * lambda2[j]))/lambda2[j]^2
  }
  if (length(h2) > 1) {
    # 2
    f2.ada[i] <- f2.ada[i]/(N * h2[1] * h2[2])
  } else {
    # 1
    f2.ada[i] <- f2.ada[i]/(N * h2[1]^2)
  }
  print(paste(i, f2.ada[i], sep = " "))
}

# Расчёт функций распределения

# фиксированный метод
F2.fix <- npudist(tdat = port.ret, edat = uv, ckertype = "gaussian", bwtype = "fixed")

# адаптивный метод
du <- u[2] - u[1]
dv <- v[2] - v[1]
w <- f2.ada
dim(w) <- c(L, L)
F2.ada <- rep(0, times = L^2)
for (i in 1:L) {
  for (j in 1:L) {
    F2.ada[j + (i - 1) * L] <- sum(w[1:j, 1:i]) * du * dv
  }
}
persp(u, v, F2.ada, theta = 30, main = "CDF estimate, 3D plot", xlab = "Gazprom returns", 
  ylab = "NorNickel returns", zlab = "CDF")
contour(u, v, F2.ada, nlevel = 7, main = "CDF estimate, contour plot", xlab = "Gazprom returns", 
  ylab = "NorNickel returns")

# Генератор случайных чисел (метод
# Монте-Карло)

# для адаптивного метода
alpha <- 0.01
M <- 5000
smpl.ind <- sample(1:(L^2), prob = f2.ada, size = M, replace = TRUE)
port.ada.sim <- uv[smpl.ind, ]
plot(port.ada.sim, xlab = "Gazprom returns", ylab = "NorNickel returns")

# Рисование графиков с перекрывающими друг
# друга точками
plot(port.ada.sim, col = rgb(0, 0, 1, alpha = 0.2))
smoothScatter(port.ada.sim) 
