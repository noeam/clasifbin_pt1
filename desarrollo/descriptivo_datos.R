rm(list = ls(all.names = TRUE))
gc()
#--------------------------------------*----------------------------------------
#------------------------------- Proyecto 1  -----------------------------------
#--------------------------------------*----------------------------------------
#datos <- read.csv("data/data_histories.csv", header = TRUE)
#variables <- c('dp_folio', 'id_sexo', 'Aedad', 'AAedad', 'AIMC', 'id_gestud',
#           'ejer_act', 'ejer1', 'ejer5', 'ejer10', 'ejer20', 'ejer30',
#           'salud_act', 'salud1', 'salud5', 'salud10', 'salud20', 'salud30',
#           'peso_act', 'peso1', 'peso5', 'peso10', 'peso20', 'peso30',
#           'peso_acc', 'peso_ehoy', 'peso_edes')
#df <- datos[,variables]
#write.csv(df, "datos_proyecto.csv")
datos <- read.csv("datos_proyecto.csv", header = TRUE)
str(datos)
summary(datos)

datos$id_sexo <- factor(datos$id_sexo,
                    levels = c('M', 'F'),
                    labels = c("hombre", "mujer"))


datos$AAedad <- factor(datos$AAedad,
                        levels = c(1,2,3,4,5,6,7,8),
                        labels = c("19-27","28-32","33-37","38-42",
                                   "43-47","48-52","53-58","59-81"))

names <- c('AIMC', 'id_gestud',
           'ejer30', 'ejer20', 'ejer10', 'ejer5', 'ejer1', 'ejer_act',
           'salud_act', 'salud1', 'salud5', 'salud10', 'salud20', 'salud30',
           'peso_act', 'peso1', 'peso5', 'peso10', 'peso20', 'peso30',
           'peso_acc', 'peso_ehoy', 'peso_edes')
datos[,names] <- lapply(datos[,names] , factor)
str(datos)
summary(datos)

#-------------------------- GRAFICAS -------------------------------------------
library(tidyverse)
library(ggpubr)

# Por sexo
(bsex <- ggplot(data = datos, aes(x = id_sexo, fill=id_sexo)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Participantes por sexo") + 
  ylab("Cantidad") +
  xlab("Sexo") + 
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_stack(0.5),
             colour="white"))

# Por grupo de edad
(bedad <- ggplot(data = datos, aes(x = AAedad, fill = AAedad)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Participantes por categoria de edad") + 
  ylab("Cantidad") +
  xlab("Categoria de edad (en años)") + 
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_stack(0.5),
             colour="white"))

# Por grado de estudios
(bestud <- ggplot(data = datos, aes(x = id_gestud, fill = id_gestud)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Participantes por grado de estudios") + 
  ylab("Cantidad") +
  xlab("Categoria de estudios") + 
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_stack(0.5),
             colour="white") +
  scale_x_discrete(limits = c("Prim", "Sec", "Bach", "CarTec", "Lic",
                              "Mast", "Doc", "PDoc", "Otro")))

# Juntamos sexo, edad y grado de estudios
ggpubr::ggarrange(bsex, bedad, bestud, 
                  labels = c("1", "2","3"),
                  ncol = 3, nrow = 1)

# Por IMC
ggplot(data = datos, aes(x = AIMC, fill = AIMC)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Participantes por grupo de Indice de Masa Corporal (IMC)") + 
  ylab("Cantidad") +
  xlab("Categoria de IMC") + 
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_stack(0.5),
             colour="white")


library(lessR)
###############################################################################
# Por grupo de edad coloreados por IMC
AAedad <- datos$AAedad
AIMC <- datos$AIMC
b1 <- BarChart(x = AAedad, by = AIMC, stack100 = TRUE)

ggplot(data = datos, aes(x = AAedad, fill = AIMC)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Participantes por grupo de edad e IMC") + 
  ylab("Cantidad") +
  xlab("Categoria de edad") + 
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_stack(0.5),
             colour="white") 

AAedad <- datos$AAedad
peso_ehoy <- datos$peso_ehoy
b <- BarChart(x = AAedad, by = peso_ehoy, stack100 = TRUE)

ggplot(data = datos, aes(x = AAedad, fill = peso_ehoy)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Participantes por grupo de edad e IMC") + 
  ylab("Cantidad") +
  xlab("Categoria de edad") + 
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_stack(0.5),
             colour="white")

AAedad <- datos$AAedad
peso_edes <- datos$peso_edes
b <- BarChart(x = AAedad, by = peso_edes, stack100 = TRUE)

ggplot(data = datos, aes(x = AAedad, fill = peso_edes)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Participantes por grupo de edad e IMC") + 
  ylab("Cantidad") +
  xlab("Categoria de edad") + 
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_stack(0.5),
             colour="white")

# Quieren bajar de peso
AAedad <- datos$AAedad
peso_acc <- datos$peso_acc
b2 <- BarChart(x = AAedad, by = peso_acc, stack100 = TRUE)

ggplot(data = datos, aes(x = AAedad, fill = peso_acc)) +
  geom_bar(stat = "count", colour = "white") +
  ggtitle("Acciones a realizar por parte de los participantes por grupo de edad") + 
  ylab("Cantidad") +
  xlab("Categoria de edad") + 
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_stack(0.5),
             colour="white")

# horas de ejercicio a la semana actualmente
ggplot(data = datos, aes(x = ejer_act, fill = AAedad)) +
  geom_bar() +
  ggtitle("Cantidad de personas que realizan ejercicio categorizado en horas por semana") + 
  ylab("Cantidad") +
  xlab("Horas de ejercicio por semana")
  #stat_count(geom = "text", 
  #           aes(label = stat(count)),
  #           position=position_stack(0.5),
  #           colour="white")
bact <- ggplot(data = datos, aes(x = ejer_act, fill = AAedad)) +
  geom_bar() +
  ylab("Cantidad") +
  xlab("Horas de ejercicio por semana")
#stat_count(geom = "text", 
#           aes(label = stat(count)),
#           position=position_stack(0.5),
#           colour="white")

b1 <- ggplot(data = datos, aes(x = ejer1, fill = AAedad)) +
  geom_bar() +
  ylab("Cantidad") +
  xlab("Horas de ejercicio por semana")

b5 <- ggplot(data = datos, aes(x = ejer5, fill = AAedad)) +
  geom_bar() +
  ylab("Cantidad") +
  xlab("Horas de ejercicio por semana")

b10 <- ggplot(data = datos, aes(x = ejer10, fill = AAedad)) +
  geom_bar() +
  ylab("Cantidad") +
  xlab("Horas de ejercicio por semana")

b20 <- ggplot(data = datos, aes(x = ejer20, fill = AAedad)) +
  geom_bar() +
  ylab("Cantidad") +
  xlab("Horas de ejercicio por semana")

b30 <- ggplot(data = datos, aes(x = ejer30, fill = AAedad)) +
  geom_bar() +
  ylab("Cantidad") +
  xlab("Horas de ejercicio por semana")

ggpubr::ggarrange(bact, b1, b5, b10, b20, b30, 
                  labels = c("Actual","1", "5","10", "20", "30"),
                  ncol = 3, nrow = 2)

ggpubr::ggarrange(b1, b2,
                  ncol = 2, nrow = 1)


