# ---------------------------------------------------------
#Pablo Benavides Herrera
# Utilizar el algoritmo de PSO para minimizar el MAPE de una
# recta con ruido
# ---------------------------------------------------------

# Paqueterias necesarias y limpiar el entorno####
rm(list=ls())
setwd(paste("C:/Users/behep/OneDrive - ITESO/PhD/",
            "Tesis/Semestre 1/Simulaciones R/",sep = ""))
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if (!require("psoptim")) install.packages("psoptim")
library(psoptim)
source("fn_pso_pb.R")

var=.03
# Función de MAPE para una recta ####
mape_recta <- function(X){
  n <- length(y)
  x <- seq(from = 1, to = n,by = 1)
  M <- X[,1]
  B <- X[,2]
  y_hat <- M%*%t(x)+B%*%t(rep(1,n)) 
  ym <- rep(1,length(X[,1]))%*%t(y)
  
  mape_error <- 1/n*rowSums(abs((ym-y_hat)/
                                  ym))
  return(-mape_error) #Se le agrega un -, porque la fn
  #psoptim maximiza
}

# Funcion montecarlo####
montecarloFun <- function(s=300,m.l=1000,w=0.85,c1=0.1,
                c2=0.1, xmin=c(-20,-20),xmax=c(20,20),
                vmax=c(4,4), seed=sample(c(1:900000),1)){
  s <- 50
  m.l <- 100
  w <- 0.85
  c1 <- 0.1
  c2 <- 0.1
  xmin <- c(0,-20)
  xmax <- c(5,20)
  vmax <- c(4, 4)

  # el algoritmo de PSO
  pso_psoptim<- pso_pb(FUN=mape_recta, n=s, max.loop=m.l,
                      w=w, c1=c1, c2=c2,xmin=xmin, 
                      xmax=xmax, vmax=vmax, seed = 
                        sample(c(1:90000),1), anim=FALSE)
  # los
  
  return(c(pso_psoptim))  
}

# Funcion creación recta ####
rectaFun <- function(q=100){
  x<-seq(from = 1, to = q,by = 1)
  
  n <- length(x)
  
  b <- 8
  
  m <- 3
  
  z<- rnorm(n,sd=var)
  y <- m*x + b + z
  return(cbind(y,x))
  
}

# Funcion MCO ####
mco<- function(x,y){
  #Calcular la regresión lineal
  lin <- lm(y~x)
  #asignar valores estimados de m ^ b
  m_mco <- lin$coefficients[2]
  b_mco <- lin$coefficients[1]
  #recta estimada de MCO-OLS
  y_mco <- m_mco*x + b_mco
  #matriz para incluir todos los valores de m^b
  X_mco <- matrix(c(m_mco, b_mco),ncol = 2)
  #calcular el RMSE de MCO-OLS
  mco_rss <- c(crossprod(lin$residuals))
  mco_mse <- mco_rss / (length(lin$residuals)-2)
  mco_rmse <- sqrt(mco_mse)
  return(c(mape_recta(X_mco),mco_rmse))
}

# Simulacion PSO ####
montecarlo_sim_mape <- c()
mco_sim_mape <- c()
y_reales <- c()
n <- 2000
start_time <- Sys.time()
for (i in 1:n){
  l <- rectaFun()
  y<-l[,1]
  x<-l[,2]
  montecarlo_sim_mape[[i]]<-data.frame(montecarloFun())
  mco_sim_mape[[i]]<-mco(x,y)
  y_reales[[i]] <- l
}

y_mapes <- c()
mape_resids <- c()
mape_rss <- c()
mape_mse <- c()

mapes_pso <-c()
rmse_pso <- c()

for (i in 1:n){
  mapes_pso[i]<-abs(montecarlo_sim_mape[[i]]$val)
  y_mapes[[i]] <- montecarlo_sim_mape[[i]]$sol.x1*x +
    montecarlo_sim_mape[[i]]$sol.x2
  mape_resids[[i]] <- y_reales[[i]][,1]-y_mapes[[i]]
  mape_rss[[i]] <- c(crossprod(mape_resids[[i]]))
  mape_mse[[i]] <- mape_rss[[i]]/(length(mape_resids[[i]])-2)
  rmse_pso[[i]] <- sqrt(mape_mse[[i]])
}

mapes_mco <- c()
rmse_mco <- c()
for (i in 1:n){
  mapes_mco[i] <- abs(mco_sim_mape[[i]][1])
  rmse_mco[i] <- mco_sim_mape[[i]][2]
}

end_time <- Sys.time()
end_time - start_time
# Graficar scatterplot de PSO y MCO para RMSE vs. MAPE ####

modelo <-  data.frame(mapes_pso,rmse_pso,mapes_mco,rmse_mco)

ggplot(data = modelo)+
    geom_point(aes(x= rmse_pso, y = mapes_pso,colour="PSO"))+
    geom_point(aes(x= rmse_mco,y = mapes_mco, colour="MCO"))+

  guides(col = guide_legend(title = ""))+
  xlab("RMSE") + ylab("MAPE")+
  scale_colour_manual(values = c(PSO="blue", MCO ="red"))
  
# # Grafica con xlim y ylim ####
# ggplot(data = modelo)+
#   geom_point(aes(x= rmse_pso, y = mapes_pso,colour="PSO"))+
#   geom_point(aes(x= rmse_mco,y = mapes_mco, colour="MCO"))+
#   
#   guides(col = guide_legend(title = ""))+
#   xlab("RMSE") + ylab("MAPE")+
#   scale_colour_manual(values = c(PSO="blue", MCO ="red"))+
#   xlim(10,30) + ylim(0,1)

# Validar diferencias entre mape y rmse ####
  
modelo2<- data.frame(mapes = c(mapes_pso,mapes_mco),
                     rmse = c(rmse_pso,rmse_mco), 
                      metodo = c(rep("PSO",n),
                                 rep("MCO",n)))
#boxplot para mapes
ggplot(data = modelo2, aes(x = metodo, y = mapes))+
  geom_boxplot(fill = "grey80", col = "blue")+
  scale_x_discrete()+ xlab("Método") + ylab("MAPES")

#boxplot para rmse
ggplot(data = modelo2, aes(x = metodo, y = rmse))+
  geom_boxplot(fill = "grey80", col = "blue")+
  scale_x_discrete()+ xlab("Método") + ylab("RMSE")

mod_mapes <- lm(modelo2$mapes ~ modelo2$metodo, 
                data = modelo2)

summary(mod_mapes)
confint(mod_mapes)
anova(mod_mapes)
confint(mod_mapes)

# par(mfrow=c(2,2))
# plot(mod_mapes)

mod_rmse <- lm(modelo2$rmse ~ modelo2$metodo, 
               data = modelo2) 
summary(mod_rmse)
anova(mod_rmse)
confint(mod_rmse)

# Código anterior a funciones ####


# Funcion mape para nlm ####

mape_recta_nlm <- function(X){
  n <- length(y)
  x <- seq(from = 1, to = n,by = 1)
  M <- X[1]
  B <- X[2]
  y_hat <- M%*%t(x)+B%*%t(rep(1,n)) 
  ym <- rep(1,length(M))%*%t(y)
  
  mape_error <- 1/n*rowSums(abs((ym-y_hat)/
                                  ym))
  return(mape_error)
}

# Funcion montecarlo nlm ####
montecarloFun <- function(b=8,m=5){
  ans=nlm(mape_recta_nlm,c(b,m))
  return(c(ans$estimate,ans$minimum
))  
}

# Simulación nlm ####
# Inicialización de variables - simulación
montecarlo_nlm <- c()
mco_sim_mape <- c()
y_reales <- c()
valores <-c()
n <- 10000 #simulaciones
# Simulación
for (i in 1:n){
  l <- rectaFun()
  y<-l[,1]
  x<-l[,2]
  ans<-montecarloFun()
  montecarlo_nlm[i]<-ans[3]
  valores[[i]]<-data.frame(ans[1:2])
  mco_sim_mape[[i]]<-mco(x,y)
  y_reales[[i]] <- l
}

# Inicialización de variables para mape y rmse - nlm
y_mapes <- c()
mape_resids <- c()
mape_rss <- c()
mape_mse <- c()
mapes_nlm <-c()
rmse_nlm <- c()
# Cálculo de MAPE y RMSE para método nlm
for (i in 1:n){
  mapes_nlm[i]<-abs(montecarlo_nlm[i])
  y_mapes[[i]] <- valores[[i]][1,]*x +
    valores[[i]][2,]
  mape_resids[[i]] <- y_reales[[i]][,1]-y_mapes[[i]]
  mape_rss[[i]] <- c(crossprod(mape_resids[[i]]))
  mape_mse[[i]] <- mape_rss[[i]]/(length(mape_resids[[i]])-2)
  rmse_nlm[[i]] <- sqrt(mape_mse[[i]])
}

# Inicialización de variables para mape y rmse - MCO
mapes_mco <- c()
rmse_mco <- c()
# Cálculo de mape y rmse para MCO
for (i in 1:n){
  mapes_mco[i] <- abs(mco_sim_mape[[i]][1])
  rmse_mco[i] <- mco_sim_mape[[i]][2]
}

# Graficar scatterplot de NLM y MCO para RMSE vs. MAPE ####

modelo <-  data.frame(mapes_nlm,rmse_nlm,mapes_mco,rmse_mco)

ggplot(data = modelo)+
  geom_point(aes(x= rmse_nlm, y = mapes_nlm,colour="NLM"))+
  geom_point(aes(x= rmse_mco,y = mapes_mco, colour="MCO"))+
  
  guides(col = guide_legend(title = ""))+
  xlab("RMSE") + ylab("MAPE")+
  scale_colour_manual(values = c(NLM="blue", MCO ="red"))


#   # Creación de la recta con ruido ####
# x<-seq(from = 1, to = 100,by = 1)
# 
# n <- length(x)
# 
# b <- 8
# 
# m <- 3
# 
# z<- rnorm(n,sd=20)
# z <- read.csv(paste("C:/Users/behep/OneDrive - ITESO/",
#         "PhD/Tesis/z.csv",sep = ""), header = FALSE)
# # z <- z$V1
# y <- m*x + b + z
# 
# r <- m*x + b

#   # PSO con paquetería "psoptim" ####

# s <- 500
# m.l <- 500
# w <- 0.85
# c1 <- 0.1
# c2 <- 0.1
# xmin <- c(-20, -20)
# xmax <- c(20, 20)
# vmax <- c(4, 4)
# 
# 
# pso_psoptim<- psoptim(FUN=mape_recta, n=s, max.loop=m.l,
#                     w=w, c1=c1, c2=c2,xmin=xmin, xmax=xmax,
#                     vmax=vmax, seed = sample(c(1:1000),1),
#                       anim=FALSE)
# 
# print(pso_psoptim)
# 


#   # Modelo de regresión lineal MCO-OLS ####
# 
# lin <- lm(y~x)
# 
# summary(lin)
# 
# m_mco <- lin$coefficients[2]
# 
# b_mco <- lin$coefficients[1]
# 
# y_mco <- m_mco*x + b_mco
# 
# X_mco <- matrix(c(m_mco, b_mco),ncol = 2)
# # Aplicar el mape a la regresión de MCO-OLS
# print(pso_psoptim)
# 
# mape_recta(matrix(c(3,8),1,2))
# 
# print(mape_recta(X_mco))
# 

#   # Graficar resultados ####
# y_real <- m*x + b
# y_pso <-  pso_psoptim$sol[,1]*x+pso_psoptim$sol[,2]
# 
# modelo <- data.frame(x,y,y_pso,y_mco,y_real)
# 
# ggplot(data = modelo, aes(x=x))+
#   geom_point(aes(y = y,colour="y"))+
#   geom_line(aes(y = y_pso, colour="y_pso"))+
# geom_line(aes(y=y_real, colour= "y_real"))+
# geom_line(aes(y=y_mco, colour="y_mco"))+
# 
# guides(col = guide_legend(title = ""))+
# ylab("")+
# scale_colour_manual(values = c(y="light blue", 
#                                       y_pso ="red"))
# y_real="green", y_mco="orange"


