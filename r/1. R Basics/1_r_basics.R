# Encoding UTF-8

# ЦМФ МГУ, математические финансы, осень 2014
# Applied Econometrics 1, R Basics
# Орлов Евгений, 20.11.2014

# ДАННЫЕ ####
sber.df <- read.csv(file.path("data", "Сбербанк [Price].txt"), header=TRUE)
gazp.df <- read.csv(file.path("data", "ГАЗПРОМ ао [Price].txt"), header=TRUE)
rosn.df <- read.csv(file.path("data", "Роснефть [Price].txt"), header=TRUE)
tail(sber.df)
# Цены закрытия начиная с 01.01.2010 по 07.11.2014 включительно
sber.close <- sber.df[sber.df$X.DATE. >= 20100101, "X.CLOSE."]
gazp.close <- gazp.df[gazp.df$X.DATE. >= 20100101, "X.CLOSE."]
rosn.close <- rosn.df[rosn.df$X.DATE. >= 20100101, "X.CLOSE."]
# Вектора доходностей
scl <- length(sber.close)
gcl <- length(gazp.close)
rcl <- length(rosn.close)
print(paste0("Размер выборки (Роснефть): ", rcl))
sber.ret <- sber.close[2:scl] / sber.close[1:(scl-1)] - 1
gazp.ret <- gazp.close[2:gcl] / gazp.close[1:(gcl-1)] - 1
rosn.ret <- rosn.close[2:rcl] / rosn.close[1:(rcl-1)] - 1

# ГРАФИКИ ####
sber.dates <- sber.df[sber.df$X.DATE. >= 20100101, "X.DATE."]
sber.dates <- as.Date(as.character(sber.dates), format="%Y%m%d")
gazp.dates <- gazp.df[gazp.df$X.DATE. >= 20100101, "X.DATE."]
gazp.dates <- as.Date(as.character(gazp.dates), format="%Y%m%d")
rosn.dates <- rosn.df[rosn.df$X.DATE. >= 20100101, "X.DATE."]
rosn.dates <- as.Date(as.character(rosn.dates), format="%Y%m%d")
# Графики цен закрытия
plot(sber.dates, sber.close, 
     type='l', xlab='index', ylab='close', main='Sberbank:Close')
plot(gazp.dates, gazp.close, 
     type='l', xlab='index', ylab='close', main='Gazprom:Close')
plot(rosn.dates, rosn.close, 
     type='l', xlab='index', ylab='close', main='Rosneft:Close')
# Графики доходностей
plot(sber.dates[-1], sber.ret, 
     type='l', xlab='index', ylab='return', main='Sberbank:Return')
plot(gazp.dates[-1], gazp.ret, 
     type='l', xlab='index', ylab='return', main='Gazprom:Return')
plot(rosn.dates[-1], rosn.ret, 
     type='l', xlab='index', ylab='return', main='Rosneft:Return')

# ПОЛЬЗОВАТЕЛЬСКАЯ ФУНКЦИЯ #####
calculate.stats <- function(returns) {
  list(mean=mean(returns), sd=sd(returns))
}
print("Сбербанк:")
calculate.stats(sber.ret)
print("Газпром:")
calculate.stats(gazp.ret)
print("Роснефть:")
calculate.stats(rosn.ret)

