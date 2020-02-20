# ------------------------------------------------------
# Minimizacion del MAPE de una recta con ruido
# Comparacion de MAPE vs. RMSE
# Determinacion de si reducir RMSE (en MCO) implica reducir
#                                   MAPE o no lo implica
# ------------------------------------------------------

# Paqueterias necesarias y limpiar el entorno####
rm(list=ls())
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)


# Datos para el analisis ####
start_time <- Sys.time() #start the timer
var <- c(seq(0.5,20, by = 0.5))
# var <- c(0.0001,0.001,0.01,0.05,0.1,0.5,seq(1,10),
#          seq(15,50,by=5))
len_var <- length(var)
iter_per_var <- 1000

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
rectaFun <- function(q = 100,var = 0.5, m = 3, b = 8){
  x <- seq(from = 1, to = q,by = 1)
  
  n <- length(x)
  
  z <- rnorm(n, sd = var)
  y <- m * x + b + z
  return(cbind(y, x))
  
}

# Funcion mape para nlm ####
mape_recta_nlm <- function(X){
  n <- length(y)
  x <- seq(from = 1, to = n,by = 1)
  M <- X[1]
  B <- X[2]
  y_hat <- M %*% t(x) + B %*% t(rep(1,n)) 
  ym <- rep(1, length(M)) %*% t(y)
  temp <- rowSums(abs((ym - y_hat)/
                ym))
  
  mape_error <- 1/n * temp
  return(mape_error)
}

# Funcion montecarlo nlm ####
montecarloFun <- function(m = 3, b = 8){
  ans=nlm(mape_recta_nlm,c(b,m))
  return(c(ans$estimate, ans$minimum))  
}

# Funcion MCO ####
mco <- function(x,y){
  #Calcular la regresi?n lineal
  lin <- lm(y ~ x)
  #asignar valores estimados de m ^ b
  m_mco <- lin$coefficients[2]
  b_mco <- lin$coefficients[1]
  #recta estimada de MCO-OLS
  y_mco <- m_mco * x + b_mco
  #matriz para incluir todos los valores de m^b
  X_mco <- matrix(c(m_mco, b_mco), ncol = 2)
  #calcular el RMSE de MCO-OLS
  mco_rss <- c(crossprod(lin$residuals))
  mco_mse <- mco_rss / (length(lin$residuals)-2)
  mco_rmse <- sqrt(mco_mse)
  return(c(mape_recta(X_mco), mco_rmse))
}

# Simulacion nlm ####
# Inicializacion de variables - simulacion
montecarlo_nlm <- c()
mco_sim_mape <- c()
y_reales <- c()
valores <-c()
n <- len_var * iter_per_var #simulaciones

# Simulacion
for (i in 1:n){
  l <- rectaFun(var = var[(i %/% (iter_per_var+0.01)+1)])
  y<-l[,1]
  x<-l[,2]
  tryCatch({
    ans<-montecarloFun()
  }, error=function(e){cat("ERROR :",conditionMessage(e), 
                                      "\n")})
  montecarlo_nlm[i] <- ans[3]
  valores[[i]] <- data.frame(ans[1:2])
  mco_sim_mape[[i]] <- mco(x,y)
  y_reales[[i]] <- l
}

# Inicializacion de variables para mape y rmse - nlm
y_mapes <- c()
mape_resids <- c()
mape_rss <- c()
mape_mse <- c()
mapes_nlm <-c()
rmse_nlm <- c()
# Calculo de MAPE y RMSE para metodo nlm
for (i in 1:n){
  mapes_nlm[i] <- abs(montecarlo_nlm[i])
  y_mapes[[i]] <- valores[[i]][1,] * x +
    valores[[i]][2,]
  mape_resids[[i]] <- y_reales[[i]][,1]-y_mapes[[i]]
  mape_rss[[i]] <- c(crossprod(mape_resids[[i]]))
  mape_mse[[i]] <- mape_rss[[i]]/(length(mape_resids[[i]])-2)
  rmse_nlm[[i]] <- sqrt(mape_mse[[i]])
}

# Inicializacion de variables para mape y rmse - MCO
mapes_mco <- c()
rmse_mco <- c()
# Calculo de mape y rmse para MCO
for (i in 1:n){
  mapes_mco[i] <- abs(mco_sim_mape[[i]][1])
  rmse_mco[i] <- mco_sim_mape[[i]][2]
}

# Graficar scatterplot de NLM y MCO para RMSE vs. MAPE ####
sd_level <- c()
for (i in 1:8){
  sd_level <- c(sd_level,rep(i,iter_per_var))
}

modelo <-  data.frame(mapes_nlm,rmse_nlm,mapes_mco,rmse_mco,
                      sd = rep(var,iter_per_var),
                      sd_level)

scatter_mape_v_rmse <- ggplot(data = modelo)+
  geom_point(aes(x= rmse_nlm, y = mapes_nlm,colour="NLM"),
             alpha = 0.1)+
  geom_point(aes(x= rmse_mco,y = mapes_mco, colour="OLS"),
             alpha = 0.1)+
  # ggtitle(paste("Simulation with SD= ",var))+
  guides(col = guide_legend(title = ""))+
  xlab("RMSE") + ylab("MAPE")+
  scale_colour_manual(values = c(NLM="blue", OLS ="red")) +
  xlim(0,25)+ylim(0,5)
scatter_mape_v_rmse

# prueba colores por SD
ggplot(data = modelo)+
  geom_point(aes(x= rmse_nlm, y = mapes_nlm,
                 colour=factor(sd_level)),
             alpha = 0.2)+
  # geom_point(aes(x= rmse_mco,y = mapes_mco, colour=sd),
  #            alpha = 0.1)+
  # ggtitle(paste("Simulation with SD= ",var))+
  guides(col = guide_legend(title = ""))+
  xlab("RMSE") + ylab("MAPE")+
  # scale_colour_manual(values = c(NLM="blue", OLS ="red")) +
  xlim(0,25)+ylim(0,0.5) +
  scale_color_hue()



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
save.image(paste("nlm10k it sd ",var,
".RData",sep = ""))


# Analisis de MAPEs y RMSEs por desv. std ####
analisis <- modelo %>%
  group_by(sd) %>% 
  summarise( 
            MAPE = mean(mapes_nlm),
            RMSE = mean(rmse_nlm))  
analisis

ggplot(data = analisis) + 
  geom_line(aes(x = sd, y = MAPE),color = "blue") + 
  ggtitle("MAPE by sd")

ggplot(data = analisis) + 
  geom_line(aes(x = sd, y = RMSE),color = "red") + 
  ggtitle("RMSE by sd")

ggplot(data = modelo) + 
  geom_boxplot(aes(x = factor(sd_level), y = mapes_nlm),
               alpha = 0.5) + 
  ggtitle("MAPE by sd")  