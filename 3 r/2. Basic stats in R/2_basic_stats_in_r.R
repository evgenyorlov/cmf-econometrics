# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied econometrics 2, Basic stats in R
# Орлов Евгений, 20.11.2014

# ДАННЫЕ ####
sber.df <- read.csv("Сбербанк [Price].txt", header=TRUE)
tail(sber.df)
# Цены закрытия начиная с 01.01.2010 по 07.11.2014 включительно
sber.close <- sber.df[sber.df$X.DATE. >= 20100101, "X.CLOSE."]
# Вектор доходностей
scl <- length(sber.close)
print(paste0("Размер выборки: ", scl))
# Вектор доходностей
sber.ret <- sber.close[2:scl] / sber.close[1:(scl-1)] - 1

sber.dates <- sber.df[sber.df$X.DATE. >= 20100101, "X.DATE."]
sber.dates <- as.Date(as.character(sber.dates), format="%Y%m%d")
# График цен закрытия
plot(sber.dates, sber.close, 
     type='l', xlab='index', ylab='close', main='Sberbank:Close')
# График доходностей
plot(sber.dates[-1], sber.ret, 
     type='l', xlab='index', ylab='return', main='Sberbank:Return')

# Тесты на нормальность
shapiro.test(sber.ret)
ks.test(sber.ret, 'pnorm', mean=mean(sber.ret), sd=sd(sber.ret))

# Тест на симметричность распределения
library(lawstat)
symmetry.test(sber.ret)

# График нормальный квантиль-квантиль
qqnorm(sber.ret)
qqline(sber.ret)

# Эмпирическая плотность распределения
plot(density(sber.ret))
rug(sber.ret)
library(fBasics)
histPlot(timeSeries(sber.ret))
