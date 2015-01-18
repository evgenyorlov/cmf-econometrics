# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Домашнее задание №1 по курсу "Advanced applied econometrics"
# Орлов Евгений, 12.11.2014

# ДАННЫЕ ####
setwd("./Downloads/programming//r/econometrics/")
sber.df <- read.csv("./raw_data/Сбербанк [Price].txt", header=TRUE)
tail(sber.df)
# Цены закрытия начиная с 01.01.2010 по 07.11.2014 включительно
sber.close <- sber.df[which(sber.df$X.DATE. >= 20100101), "X.CLOSE."]
# Вектор доходностей
scl <- length(sber.close)
scl
sber.ret <- sber.close[2:scl] / sber.close[1:(scl-1)] - 1
# График цен закрытия
plot(sber.close, type='l', xlab='index', ylab='close', main='Сбербанк, ао')
# График доходностей
plot(sber.ret, type='l', xlab='index', ylab='return', main='Сбербанк, ао')

# Тесты на нормальность
plot(density(sber.ret))
rug(sber.ret)
qqnorm(sber.ret)
qqline(sber.ret)
shapiro.test(sber.ret)
ks.test(sber.ret, 'pnorm', mean=mean(sber.ret), sd=sd(sber.ret))

# Сводная статистическая информация
library(fBasics)
basicStats(sber.ret)
histPlot(timeSeries(sber.ret))
acf(sber.ret)

# Тест на симметричность распределения
library(lawstat)
symmetry.test(sber.ret)

# VaR-КРИВАЯ (ОБОБЩЕННОЕ ГИПЕРБОЛИЧЕСКОЕ РАСПРЕДЕЛЕНИЕ) ####
var.ghyp <- numeric()


# VaR-КРИВАЯ (ОБОБЩЕННОЕ РАСПРЕДЕЛЕНИЕ ПАРЕТО) ####


# VaR-КРИВАЯ (НЕПАРАМЕТРИЧЕСКОЕ МОДЕЛИРОВАНИЕ) ####


# VaR-КРИВАЯ (GARCH-МОДЕЛЬ) ####


# ПРОВЕРКА КАЧЕСТВА ОЦЕНОК РИСКА ####

kupic.test <- function() {
  # Тест Купика
  #
}
