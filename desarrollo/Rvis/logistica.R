rm(list = ls(all.names = TRUE))
gc()
#--------------------------------------*----------------------------------------
#------------------------------- Proyecto 1  -----------------------------------
#--------------------------------------*----------------------------------------
datos <- read.csv("~/GitHub/Proyecto/desarrollo/python/datos_patronesejercicio.csv", header = TRUE)
str(datos) # 1067 observaciones
summary(datos)

# Tratamos a la mayoria de las variables como categoricas/factor
datos$id_sexo <- factor(datos$id_sexo,
                        levels = c('M', 'F'),
                        labels = c("hombre", "mujer"))
datos$id_gestud <- factor(datos$id_gestud)
datos$ejer_0 <- factor(datos$ejer_0)
datos$ejer_1 <- factor(datos$ejer_1)
datos$ejer_5 <- factor(datos$ejer_5)
datos$ejer_10 <- factor(datos$ejer_10)
datos$ejer0_10 <- factor(datos$ejer0_10)
datos$obesity <- factor(datos$obesity)
str(datos)
summary(datos)

ggplot(data = datos, aes(x = ejer0_10, fill=obesity)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Frecuencia de patrones") + 
  ylab("Ocurrencias") +
  xlab("Patrón") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = datos, aes(x = obesity, fill=ejer_0)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Participantes por sexo") + 
  ylab("Cantidad") +
  xlab("Clase Obesidad o No Obesidad") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = datos, aes(x = obesity, fill=ejer_1)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Participantes por sexo") + 
  ylab("Cantidad") +
  xlab("Clase Obesidad o No Obesidad") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

library(lessR)
###############################################################################
# Por grupo de edad coloreados por IMC
obesidad <- datos$obesity
ejer0_10 <- datos$ejer0_10
ejer_0 <- datos$ejer_0
ejer_1 <- datos$ejer_1
BarChart(x = obesidad, by = ejer0_10, stack100 = TRUE)
BarChart(x = obesidad, by = ejer_0, stack100 = TRUE)
BarChart(x = ejer_0, by = obesidad, stack100 = TRUE)
BarChart(x = obesidad, by = ejer_1, stack100 = TRUE)



#--------------------------------------*----------------------------------------
#---------------------- Conteos Naive   ----------------------------------------
#--------------------------------------*----------------------------------------
nrow(datos[datos$obesity == 1, ])
nrow(datos[datos$ejer_0 == "B", ])
nrow(datos[datos$ejer_1 == "M", ])
nrow(datos[datos$ejer_5 == "M", ])
nrow(datos[datos$ejer_10 == "M", ])
nrow(datos[datos$obesity == 1 & datos$ejer_0 == "B", ])
nrow(datos[datos$obesity == 1 & datos$ejer_1 == "M", ])
nrow(datos[datos$obesity == 1 & datos$ejer_5 == "M", ])
nrow(datos[datos$obesity == 1 & datos$ejer_10 == "M", ])
nrow(datos[datos$obesity == 0, ])


#--------------------------------------*----------------------------------------
#---------------------- Regresion Logística  -----------------------------------
#--------------------------------------*----------------------------------------
# 2/3 para entrenar y 1/3 para pruebas
train <- datos[1:711,]
test <- datos[712:1067,]

#--------------------------------------*----------------------------------------
#---------------------------- Primer Intento -----------------------------------
#----------------------- Modelo con interacciones ------------------------------
#--------------------------------------*----------------------------------------
# La variable referencia es A
fit <- glm(obesity ~ ejer_0 * ejer_1 * ejer_5 * ejer_10, family=binomial(link='logit'), data=train)
summary(fit)

anova(fit, test="Chisq")

fitted.results <- predict(fit, newdata = subset(test, select=c("ejer_0", "ejer_1", "ejer_5", "ejer_10")), type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test$obesity)
print(paste('Accuracy', 1 - misClasificError)) # Accuracy 0.806179775280899



library(ROCR)
#p <- predict(fit, newdata=subset(test,select=c("ejer_0", "ejer_1", "ejer_5", "ejer_10")), type="response")
p <- predict(fit, test, type="response")
pr <- prediction(p, test$obesity)
pr <- prediction(p, test$obesity)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf, main="Curva ROC mod con interacciones")
abline(coef = c(0, 1),
       col = "red",
       lwd = 1)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc ### Area bajo la curva ROC 0.5791799
legend("bottomright", inset=.02, title="AUC", legend = c(auc), col = c(auc))

#--------------------------------------*----------------------------------------
#--------------------------- Segundo  Intento ----------------------------------
#----------------------- Modelo sin interacciones ------------------------------
#--------------------------------------*----------------------------------------


modelo <- glm(obesity ~ ejer_0 + ejer_1 + ejer_5 + ejer_10, family=binomial(link='logit'), data=train)
summary(modelo)

anova(modelo, test="Chisq")

#fitted.results <- predict(modelo, newdata = subset(test, select=c("ejer_0", "ejer_1", "ejer_5", "ejer_10")), type='response')
fitted.results <- predict(modelo, test, type="response")
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test$obesity)
print(paste('Accuracy', 1 - misClasificError)) # "Accuracy 0.806179775280899"



library(ROCR)
#p <- predict(fit, newdata=subset(test,select=c("ejer_0", "ejer_1", "ejer_5", "ejer_10")), type="response")
p <- predict(modelo, test, type="response")
pr <- prediction(p, test$obesity)
pr <- prediction(p, test$obesity)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf, main="Curva ROC mod sin interacciones")
abline(coef = c(0, 1),
       col = "red",
       lwd = 1)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc ### Area bajo la curva ROC 0.5823613 vs 0.5791799 del anterior
legend("bottomright", inset=.02, title="AUC", legend = c(auc), col = c(auc))


#--------------------------------------*----------------------------------------
#--------------------------- Tercer  Intento -----------------------------------
#----------------------- Modelo solo con habito actual -------------------------
#--------------------------------------*----------------------------------------


fit2 <- glm(obesity ~ ejer_0 + ejer_1 + ejer_5 + ejer_10 +id_sexo + id_gestud, family=binomial(link='logit'), data=train)
summary(fit2)

anova(fit2, test="Chisq")

#fitted.results <- predict(fit2, newdata = subset(test, select=c("ejer_0", "ejer_1", "ejer_5", "ejer_10")), type='response')
fitted.results <- predict(fit2, test, type="response")
fitted.results <- ifelse(fitted.results > 0.5,1,0)
misClasificError <- mean(fitted.results != test$obesity)
print(paste('Accuracy', 1 - misClasificError)) # "Accuracy 0.808988764044944" vs "Accuracy 0.806179775280899" 



library(ROCR)
#p <- predict(fit, newdata=subset(test,select=c("ejer_0", "ejer_1", "ejer_5", "ejer_10")), type="response")
p <- predict(fit2, test, type="response")
pr <- prediction(p, test$obesity)
pr <- prediction(p, test$obesity)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf, main="Curva ROC mod sin interacciones más sexo y grado de estudios")
abline(coef = c(0, 1),
       col = "red",
       lwd = 1)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc ### Area bajo la curva ROC 0.6177347 vs 0.5545877 vs 0.5823613 vs 0.5791799 del anterior

legend("bottomright", inset=.02, title="AUC", legend = c(auc), col = c(auc))


#library(bestglm)
# Prepare data
#X.matrix <- as.data.frame(model.matrix(obesity ~ ejer_0 + ejer_1 + ejer_5 + ejer_10 + id_sexo + id_gestud, datos))

## Perform
#res.best.logistic <-
#  bestglm(X.matrix,
#         family = binomial,          # binomial family for logistic
#          IC = "AIC",                 # Information criteria for
#          method = "exhaustive")

