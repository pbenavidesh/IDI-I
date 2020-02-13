# ------------------------------------------------------
# Minimizacion del MAPE de una recta con ruido
# Comparacion de MAPE vs. RMSE
# Determinacion de si reducir RMSE (en MCO) implica reducir
#                                   MAPE o no lo implica
# ------------------------------------------------------

# Paqueterias necesarias y limpiar el entorno####
rm(list=ls())
setwd("C:/Users/behep/OneDrive - ITESO/PhD")
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)


# Datos para el analisis ####
start_time <- Sys.time() #start the timer
# var <- 30
var <- c(0.0001,0.001,0.01,0.05,0.1,0.5,seq(1,10),
         seq(15,30,by=5),50,100)
# Funcion de MAPE para una recta ####
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
  #                                 fn psoptim maximiza
}

# Funcion creacion recta ####
rectaFun <- function(q=100,var=0.5){
  x<-seq(from = 1, to = q,by = 1)
  
  n <- length(x)
  
  b <- 8
  
  m <- 3
  
  z<- rnorm(n,sd=var)
  y <- m*x + b + z
  return(cbind(y,x))
  
}

# Funcion mape para nlm ####
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

# Funcion montecarlo nlm ####
montecarloFun <- function(b=8,m=5){
  ans=nlm(mape_recta_nlm,c(b,m))
  return(c(ans$estimate,ans$minimum
  ))  
}

# Funcion MCO ####
mco <- function(x,y){
  #Calcular la regresi?n lineal
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

# Simulacion nlm ####
# Inicializacion de variables - simulacion
montecarlo_nlm <- c()
mco_sim_mape <- c()
y_reales <- c()
valores <-c()
n <- 10000 #simulaciones

# Simulaci?n
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

# Inicializaci?n de variables para mape y rmse - nlm
y_mapes <- c()
mape_resids <- c()
mape_rss <- c()
mape_mse <- c()
mapes_nlm <-c()
rmse_nlm <- c()
# C?lculo de MAPE y RMSE para m?todo nlm
for (i in 1:n){
  mapes_nlm[i]<-abs(montecarlo_nlm[i])
  y_mapes[[i]] <- valores[[i]][1,]*x +
    valores[[i]][2,]
  mape_resids[[i]] <- y_reales[[i]][,1]-y_mapes[[i]]
  mape_rss[[i]] <- c(crossprod(mape_resids[[i]]))
  mape_mse[[i]] <- mape_rss[[i]]/(length(mape_resids[[i]])-2)
  rmse_nlm[[i]] <- sqrt(mape_mse[[i]])
}

# Inicializaci?n de variables para mape y rmse - MCO
mapes_mco <- c()
rmse_mco <- c()
# C?lculo de mape y rmse para MCO
for (i in 1:n){
  mapes_mco[i] <- abs(mco_sim_mape[[i]][1])
  rmse_mco[i] <- mco_sim_mape[[i]][2]
}

# Graficar scatterplot de NLM y MCO para RMSE vs. MAPE ####

modelo <-  data.frame(mapes_nlm,rmse_nlm,mapes_mco,rmse_mco)

scatter_mape_v_rmse <- ggplot(data = modelo)+
  geom_point(aes(x= rmse_nlm, y = mapes_nlm,colour="NLM"))+
  geom_point(aes(x= rmse_mco,y = mapes_mco, colour="OLS"))+
  ggtitle(paste("Simulation with SD= ",var))+
  guides(col = guide_legend(title = ""))+
  xlab("RMSE") + ylab("MAPE")+
  scale_colour_manual(values = c(NLM="blue", OLS ="red")) #+
  # xlim(20,40)+ylim(0,25)
scatter_mape_v_rmse
# Validar diferencias significativas MAPE vs. RMSE ####
modelo2<- data.frame(mapes = c(mapes_nlm,mapes_mco),
                     rmse = c(rmse_nlm,rmse_mco), 
                     metodo = c(rep("NLM",n),
                                rep("OLS",n)))
#boxplot para mapes
boxplot_mapes <- ggplot(data = modelo2, aes(x = metodo,
                                                y = mapes))+
  geom_boxplot(fill = "grey80", col = "blue")+
  scale_x_discrete()+ xlab("Method") + ylab("MAPES")
boxplot_mapes
#boxplot para rmse
boxplot_rmse <- ggplot(data = modelo2, aes(x = metodo,
                                              y = rmse))+
  geom_boxplot(fill = "grey80", col = "blue")+
  scale_x_discrete()+ xlab("Method") + ylab("RMSE")
boxplot_rmse

mod_mapes <- lm(modelo2$mapes ~ modelo2$metodo, 
                data = modelo2)

# The linear regression
# summary(mod_mapes)
xtable(summary(mod_mapes),display = c("s","e",
                                "e","f","g"),
       caption = paste("Linear regression with SD=",var),
       label = paste("lin_reg_sd_",var,sep = ""))

# Confidence intervals
confint(mod_mapes)



# Welch t-tests
t.test(mapes_nlm,mapes_mco)
t.test(rmse_nlm,rmse_mco)

end_time <- Sys.time()
end_time - start_time
# Guardar los resultados ####
save.image(paste("C:/Users/behep/OneDrive - ITESO/PhD/",
"Tesis/Semestre 1/Simulaciones R/nlm10k it sd ",var,
".RData",sep = ""))
