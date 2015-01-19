library("np", quietly=TRUE)
f.fix <- npudens(tdat=y, edat=x, ckertype="gaussian", bwtype="fixed")

# Одномерный случай
# Пример 1. Острова

# Адаптивный метод с λi
pilot <- npudens(tdat=y, ckertype="gaussian", bwtype="fixed")
h <- pilot$bws$bw  # оценка глобальной составляющей интервала

# среднегеометрическое пилотных оценок
g <- 1
for (i in 1:N) {
  g <- g*pilot$dens[i]^(1/N) 
}

# расчёт локальной концентрации наблюдений
alpha <- 0.5
lambda <- (g/pilot$dens)^alpha
kern <- function(u) {
  exp(-u^2/2)/sqrt(2*pi)  # ядро Гаусса
}

# расчёт оценок плотности
f <- numeric(L)
for (i in 1:L) {
  f[i] <- sum(kern((x[i]-y)/(h*lambda))/(h*lambda)
}
f <- f / N

# Сравнение адаптивной и фиксированной оценок
plot(x, f.fix$dens,type="l", lty="dashed", ylim=c(0, 0.4), 
     main="Fixed and adaptive estimates", xlab="y",ylab="Density")
lines(x, f)

# Нахождение квантилей оценки распределения, F −1 α

# оценка функции распределения
F.fix <- npudist(tdat=y,edat=x,ckertype="gaussian",bwtype="fixed")

# для адаптивного варианта
F <- rep(0, times=L)
for (i in 1:L) {
  F[i] <- sum(f[1:i])*dx
}

# поиск квантиля методом деления пополам
alpha <- 0.99
a <- 1
b <- L
ab <- trunc((a+b)/2)
while ((b-a) > 2) {
  if (F.fix$dist[ab] <= alpha) {
    a <- ab
  }
  if (F.fix$dist[ab] >= alpha) {
    b <- ab
  }
  ab <- trunc((a+b)/2)
}
q.fix <- x[ab]

# Генератор случайных чисел

# фиксированный интервал
M <- 10^6
y.fix.sim <- sample(x, prob=f.fix$dens, size=M, replace=TRUE)
q.fix <- sort(y.fix.sim)[alpha*M]

# для адаптивного варианта
y.ada.sim <- sample(x, prob=f, size=M, replace=TRUE)
q.ada <- sort(y.ada.sim)[alpha*M]

# Многомерный случай
# Пример 2. Старый служака
y <- faithful
N <- nrow(y)

# сетка для расчёта оценок плотности
L <- 50
u <- seq(0, 7, length=L)
v <- seq(30, 110, length=L)
uv <- expand.grid(u, v)

# оценка плотности
f.fix <- npudens(tdat=y, edat=uv, ckertype="gaussian", bwtype="fixed")

# графики оценки
w <- f.fix$dens
dim(w) <- c(L, L)
persp(u, v, w, theta=30, main="Bivariate kernel estimate, 3D plot", 
      xlab="Eruption time", ylab="Waiting time", zlab="Density")

contour(u, v, w, nlevel=7, 
        main="Bivariate kernel estimate, contour plot",
        xlab="Eruption time", ylab="Waiting time")

# Адаптивный метод с λi, аналогично одномерному случаю
pilot <- npudens(tdat=y, ckertype="gaussian", bwtype="fixed")
h <- pilot$bws$bw

g <- 1
for (i in 1:N) {
  g <- g*pilot$dens[i]^(1/N)
}

alpha <- 0.5
lmbd <- (g/pilot$dens)^alpha

kern <- function(x) {
  exp(-(x[1]^2+x[2]^2)/2)/(2*pi)
}

f <- rep(0, times=L^2)
for (i in 1:(L^2)) {
  for (j in 1:N) {
    f[i] <- f[i]+kern((uv[i,]-y[j,])/(h*lmbd[j]))/lmbd[j]^2  
  }
  f[i] <- f[i]/(N*h[1]*h[2])
}

# Расчёт функций распределения

# фиксированный метод
F.fix <- npudist(tdat=y, edat=uv, ckertype="gaussian", bwtype="fixed")

# адаптивный метод
du <- u[2] - u[1]
dv <- v[2] - v[1]
w <- f
dim(w) <- c(L, L)
F <- rep(0, times=L^2)
for (i in 1:L) {
  for (j in 1:L) {
    F[j+(i-1)*L] <- sum(w[1:j, 1:i])*du*dv
  }
}

# Генератор случайных чисел

# для адаптивного метода
alpha <- 0.99
M <- 5000
smpl.ind <- sample(1:(L^2), prob=f, size=M, replace=TRUE)
y.ada.sim <- uv[smpl.ind, ]
plot(y.ada.sim, xlab="Eruption", ylab="Waiting time")

# Рисование графиков с перекрывающими друг друга точками
plot(y.ada.sim, col=rgb(0, 0, 1, alpha=0.2))
smoothScatter(y.ada.sim)
