# Otras funciones de PSO

# #Con paquetería "metaheuristicOpt" ####
# 
# if (!require("metaheuristicOpt")) install.packages("metaheuristicOpt")
# library(metaheuristicOpt)
# 
# numvar <- 5
# ci <- 0.1
# cg <- 0.1
# numpop <- 1000
# rangeVar <- matrix(c(-10,10), nrow=2)
# vmax <- 2
# w <- 0.5
# maxiter <- 1000
# 
# pso_meta <- PSO(mape_recta, optimType = "MIN", numVar = numvar, 
#                 numPopulation = numpop, maxIter = maxiter,rangeVar, 
#                 Vmax = vmax, ci = ci, cg = cg, w = w)
# # Error in X[, 1] : incorrect number of dimensions


#Con paquetería "pso" ####
# detach("package:psoptim", unload=TRUE)
# 
# if (!require("pso")) install.packages("pso")
# library(pso)
# 
# 
# pso_pso <- psoptim(rep(NA,2), fn =  mape_recta,lower = 0, 
#                  upper =  20,control = list(trace=1,
#                                  REPORT=1, trace.stats=TRUE, s=2000) )
# 
# 
# 



# PSO de Ponchito ####

np<-20; #N?mero de particulas
#inicializaci?n
x1p<-list()
for(j in 1:length(seq(np))){
  x1p[[j]]<-c(0,0,0)
  
} 

for(j in 1:length(seq(np))){
  x1p[[j]][1]<-runif(1, min=-.3, max=0)
  x1p[[j]][2]<-runif(1, min=0, max=.5)
  x1p[[j]][3]<-runif(1, min=0, max=1)
} 

x1p[[1]]<-c(-.03,.2,.25)
x1pg<-c(0,0,0)
vx1<-list()
for(j in 1:length(seq(np))){
  vx1[[j]]<-c(0,0,0)
}
x1pL<-x1p

fxpg<-1000 #desempe?o valor inicial del mejor global
fxpL<-list()
for(j in 1:length(seq(np))){
  fxpL[[j]]<-c(fxpg) #desempe?o delos mejores locales
}
c1<-0.3 #Velocidad de convergencia al  mejor global
c2<-0.3 #velocidad de convergencia al mejor local
#iteraciones
for(k in 1:length(seq(20))){
  fx<-list()
  a<- -1000
  for(i in 1:length(seq(np))){
    t<-trading_strategy(Historico,x1p[[i]][1],x1p[[i]][2],x1p[[i]][2])
    fx[[i]]<- -(t[[2]]+a*max(x1p[[i]][1],0)+a*max(-x1p[[i]][2],0)+a*max(x1p[[i]][2]-1,0)+a*max(-x1p[[i]][3],0)+a*max(x1p[[i]][3]-1,0))
  }
  ind<-which.min(fx)
  val<-fx[[ind]]
  if(val<fxpg){
    x1pg<-x1p[[ind]]
    fxpg<-val;
  }
  for(p in 1:seq((length(np)))){
    if(fx[[p]]<fxpL[[p]]){
      x1pL[[p]]<-x1p[[p]]
    }
  }
  for(p in 1:seq(length(np))){
    vx1[[p]]=vx1[[p]]+c1*runif(3, min=0, max=1)*(x1pg-x1p[[p]])+c2*runif(3, min=0, max=1)*(x1pL[[p]]-x1p[[p]])
  } 
}
optime_result<-trading_strategy(Historico,x1pg[1],x1pg[2],x1pg[3])
toc()


#mejor resultado segun corrida en iteso  #1
result<-trading_strategy(Historico,-.0455,.4033,.2691)