# ------------------------------------------------------
# Minimización del MAPE de una recta con ruido
# Comparación de MAPE vs. RMSE
# Determinación de si reducir RMSE (en MCO) implica reducir
#                                   MAPE o no lo implica
# ------------------------------------------------------
# install.packages("ggplot2", dep = TRUE) para evitar un
# warning
# Paqueterias necesarias y limpiar el entorno####
rm(list=ls())
setwd("C:/Users/behep/OneDrive - ITESO/PhD")
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if (!require("ggpubr")) install.packages("ggpubr")
library(ggpubr)

# Datos para el análisis ####
start_time <- Sys.time() #start the timer
var <- 30
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
  return(-mape_error) #Se le agrega un -, porque la
  #                               fn psoptim maximiza
}

# Función creación recta ####
rectaFun <- function(q=100,var=0.5){
  x<-seq(from = 1, to = q,by = 1)
  
  n <- length(x)
  
  b <- 8
  
  m <- 3
  
  z<- rnorm(n,sd=var)
  y <- m*x + b + z
  return(cbind(y,x))
  
}

# Función mape para nlm ####
mape_recta_nlm <- function(X){
  n <- length(y)
  x <- seq(from = 1, to = n,by = 1)
  M <- X[1]
  B <- X[2]
  y_hat <- M%*%t(x)+B%*%t(rep(1,n)) 
  ym <- rep(1,length(M))%*%t(y)
  temp<-rowSums(abs((ym-y_hat)/
                ym))
  
  mape_error <- 1/n*temp
  return(mape_error)
}

# Función montecarlo nlm ####
montecarloFun <- function(b=8,m=5){
  ans=nlm(mape_recta_nlm,c(b,m))
  return(c(ans$estimate,ans$minimum
  ))  
}

# Función MCO ####
mco <- function(x,y){
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
  return(c(mape_recta(X_mco),mco_rmse,mco_mse))
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
  l <- rectaFun(var = var)
  y<-l[,1]
  x<-l[,2]
  tryCatch({
    ans<-montecarloFun()
  }, error=function(e){cat("ERROR :",conditionMessage(e), 
                                                "\n")})
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
mse_mco <- c()
# Cálculo de mape y rmse para MCO
for (i in 1:n){
  mapes_mco[i] <- abs(mco_sim_mape[[i]][1])
  rmse_mco[i] <- mco_sim_mape[[i]][2]
  mse_mco[i] <- mco_sim_mape[[i]][3]
}

# Graficar scatterplot de NLM y MCO para RMSE vs. MAPE ####

modelo <-  data.frame(mapes_nlm,rmse_nlm,mape_mse,
                      mapes_mco,rmse_mco,mse_mco)

scatter_mape_v_rmse <- ggplot(data = modelo)+
  geom_point(aes(x= rmse_nlm, y = mapes_nlm,colour="NLM"))+
  geom_point(aes(x= rmse_mco,y = mapes_mco, colour="MCO"))+
  ggtitle(paste("Análisis con Desv. Std= ",var))+
  guides(col = guide_legend(title = ""))+
  xlab("RMSE") + ylab("MAPE")+
  scale_colour_manual(values = c(NLM="blue", MCO ="red")) #+
  # xlim(20,40)+ylim(0,25)
scatter_mape_v_rmse

scatter_mape_v_mse <- ggplot(data = modelo)+
  geom_point(aes(x= mape_mse, y = mapes_nlm,colour="NLM"))+
  geom_point(aes(x= mse_mco,y = mapes_mco, colour="MCO"))+
  ggtitle(paste("Análisis con Desv. Std= ",var))+
  guides(col = guide_legend(title = ""))+
  xlab("MSE") + ylab("MAPE")+
  scale_colour_manual(values = c(NLM="blue", MCO ="red")) #+
# xlim(20,40)+ylim(0,25)
scatter_mape_v_mse

# Validar diferencias significativas MAPE vs. RMSE ####
modelo2<- data.frame(mapes = c(mapes_nlm,mapes_mco),
                     rmse = c(rmse_nlm,rmse_mco), 
                     mse = c(mape_mse,mse_mco),
                     metodo = c(rep("NLM",n),
                                rep("MCO",n)))
#boxplot para mapes
boxplot_mapes <- ggplot(data = modelo2, aes(x = metodo,
                                              y = mapes))+
  geom_boxplot(fill = "grey80", col = "blue")+
  scale_x_discrete()+ xlab("Método") + ylab("MAPES")
boxplot_mapes
#boxplot para rmse
boxplot_rmse <- ggplot(data = modelo2, aes(x = metodo,
                                              y = rmse))+
  geom_boxplot(fill = "grey80", col = "blue")+
  scale_x_discrete()+ xlab("Método") + ylab("RMSE")
boxplot_rmse

#boxplot para mse
boxplot_mse <- ggplot(data = modelo2, aes(x = metodo,
                                           y = mse))+
  geom_boxplot(fill = "grey80", col = "blue")+
  scale_x_discrete()+ xlab("Método") + ylab("MSE")
boxplot_mse


mod_mapes <- lm(modelo2$mapes ~ modelo2$metodo, 
                data = modelo2)

summary(mod_mapes)
confint(mod_mapes)

t.test(mapes_nlm,mapes_mco)
t.test(rmse_nlm,rmse_mco)
t.test(mape_mse,mse_mco)

end_time <- Sys.time()
end_time - start_time
# Guardar los resultados ####
save.image(paste("C:/Users/behep/OneDrive - ITESO/PhD/",
"Tesis/Semestre 1/Simulaciones R/nlm10k it sd ",var,
".RData",sep = ""))