##################################################
### Transformaciones 
##################################################
# Para ejemplificar esto estudiaremos la dependencia 
# entre la altura y el diametro de los cedro (arboles)
# de cierta region de Estados Unidos. 

library(alr3)
library(car)
data(ufcgf)
attach(ufcgf)

# Dbh: diametro del arbol en milimetros (mm)
# Height: altura del arbol den decimetros (dm). 

plot(Dbh,Height,pch=20,cex=1.8)
lines(lowess(Height~Dbh,f=2/3,iter=1), lty=2,lwd=2,col="red")

# Se buscara con transformaciones, 
# lograr un modelo ajustado lo mas parecido posible a la linea roja, 
# con las variables transformadas `adecuadamente'.

lhat <- inv.tran.estimate(Dbh,Height,family="box.cox")
DbhPower <- powtran(Dbh,lhat$lambda, family="bcPower")
DbhPowerScaled <- powerTransform(Dbh,lhat$lambda,family="box.cox")
par(mfrow=c(1,2))
plot(DbhPower,Height,pch=20,cex=1.8)
lines( lowess(Height~DbhPower,f=2/3,iter=1), lty=2,lwd=2,col="red")
box(lwd=2,col='black')
plot(DbhPowerScaled,Height,pch=20,cex=1.8)
lines( lowess(Height~DbhPowerScaled,f=2/3,iter=1), lty=2,lwd=2,col="red")

lam <-c(1,0,-1) 
inv.tran.plot(Dbh,Height,lambda = lam,pch=20,cex=1.8, col=c(1,2,3,4))

##################################################

library(MASS)
m1<-lm(Height~Dbh)
boxcox(m1,xlab=expression(lambda[y])) # transformar Y
boxTidwell(Height~Dbh) # transformar X
lbox<-1.35 
Hbox<-powtran(Height,lbox,family="box.cox",modified=T)  
plot(Dbh,Hbox,pch=20,cex=1.8,font.lab=7) 
lines(lowess(Hbox~Dbh,f=2/3,iter=1), lty=2,lwd=2,col="red") 

##################################################
# Analisis de Residuos

qqnorm(m1$res,pch=19,cex=1)
qqline(m1$res)

plot(m1$fit,rstudent(m1),pch=19)


ad.test(m1$res)
sf.test(m1$res)
bptest(m1)
ncvTest(m1)
bgtest(m1)
dwtest(m1)
durbinWatsonTest(m1)

##################################################

