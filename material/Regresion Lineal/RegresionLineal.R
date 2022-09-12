###########################################################################
### Regresion Lineal Simple
###########################################################################

###########################################################################
### Datos simulados
dev.new(width=4,height=4)

n = 100 ; B0 = 1 ; B1 = 5 ; sigma2 = 2.5 ;   # valores de los parametros
set.seed(12345)                            # semilla
x = runif(n, 0,3)                          # variable explicativa
Eps = rnorm(n, mean=0,sd=sqrt(sigma2))     # errores
y = B0 + B1*x + Eps                        # variable respuesta
plot(x,y)                                  # grafica los puntos (x,y)
modelo <- lm(y ~ x)                        # regresion lineal simple
lines(x,modelo$fit, lwd=2)                 # grafica la linea ajustada
summary(modelo)

###########################################################################
### Residuales

idx = c(1,35,65,99)
x1 = x[order(x)]                      # ordena x, guarda en x1 las 
x1 = x1[idx]                          # observaciones     
y1 = y[order(x)]                      # ordena y con base en x, guarda en 
y1 = y1[idx]                          # y1 las observaciones
fit1 = modelo$fitted.values[order(x)] # ordena los datos ajustados con base
fit1 = fit1[idx]                      # en x, guarda en fit1 las 
                                      # observaciones

plot(x1,y1, xlab="x",ylab="y",        # grafica (x1,y1)
  xlim=c(min(x)-0.5,max(x)+0.5),ylim=c(min(y),max(y)+1))
lines(x,modelo$fitted.values)         # grafica la linea ajustada
abline(reg = modelo)
for(i in 1:4){                        # grafica las lineas entre los puntos
  lines(c(x1[i],x1[i]),c(y1[i],fit1[i]))   # (x1,y1) y (x1,fit1)     
}
points(x1,fit1, pch=19)               # grafica los puntos (x1,fit1)
text(x1,fit1, labels=round(fit1,1),pos=4)  # grafica la coordenada fit1
text(x1,y1, labels=round(y1,1),pos=4)      # grafica la coordenada y1

plot(x1,y1, xlab="x",ylab="y",        # grafica (x1,y1)
  xlim=c(min(x)-0.5,max(x)+0.5),ylim=c(min(y),max(y)+1))
lines(x,modelo$fitted.values)         # grafica la linea ajustada
abline(reg = modelo)
for(i in 1:4){                        # grafica las lineas entre los puntos
  lines(c(x1[i],x1[i]),c(y1[i],fit1[i]))   # (x1,y1) y (x1,fit1)  
  text(x1[i],(fit1[i]+y1[i])/2,       # grafica el valor y1-fit1
    labels=bquote(epsilon[.(idx[i])]==.(round(y1[i]-fit1[i],1))),pos=3)
}

###########################################################################

x2 = data.frame(x=x[order(x)])
conf <- predict(modelo,newdata=x2, interval="confidence")
pred <- predict(modelo,newdata=x2, interval="prediction")

dev.new(width=4,height=4)
plot(x,y, xlab="x",ylab="y", main="", cex=0.5)
lines(x2$x,conf[,1], lwd=2)
lines(x2$x,conf[,2], lwd=2,lty=2)
lines(x2$x,conf[,3], lwd=2,lty=2)
lines(x2$x,pred[,2], lwd=2,lty=3)
lines(x2$x,pred[,3], lwd=2,lty=3)

legend("topleft", legend=c("Ajuste","Intervalo de E[Y|X]","Intervalo de Y*"), cex=0.7,lty=c(1,2,3),lwd=2)

###########################################################################
### Ejemplo

set.seed(12345)
n = 100
x = runif(n, 0, 1)
x2 = rnorm(n,0,1)
B0 = 2
B1 = 3
sigma2 = 1
E = rnorm(n, 0, sqrt(sigma2))
y = B0 + B1*x + E

plot(x,y)

modelo <- lm(y ~ x)
abline(modelo, col="blue", lwd=3)

summary(modelo)

confint(modelo, parm=c("(Intercept)","x"), level=0.95)

anova(modelo)

predict(modelo)
points(x,predict(modelo), col="red", pch=19)

#x2 = data.frame(x=c(0.2,0.4,0.6,0.8))
x2 = data.frame(x=x[order(x)])
conf <- predict(modelo, newdata=x2, interval="confidence", level=0.95)
pred <- predict(modelo, newdata=x2, interval="prediction", level=0.95)

plot(x,y, xlab="x",ylab="y", main="", cex=0.5)
lines(x2$x,conf[,1], lwd=2, col="red")
lines(x2$x,conf[,2], lwd=2,lty=2, col="blue")
lines(x2$x,conf[,3], lwd=2,lty=2, col="blue")
lines(x2$x,pred[,2], lwd=2,lty=3, col="magenta")
lines(x2$x,pred[,3], lwd=2,lty=3, col="magenta")

###########################################################################
###########################################################################
###########################################################################
###########################################################################


