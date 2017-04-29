# RF_y_Parsimonia-
..
## Alejandra Ruiz.
## 25-Abril-2017
## Editado, corregido y mejorado por el profesor Daniel Miranda-Esquivel
## Con RF se comparan arboles de parsimonia lineal y piwi (2,5,9,18,36,72) a partir de matrices (con muestreo incompleto de genes).    

### Inicio. 
rm(list=ls())

listadoTopologias <- read.csv("/media/aleu/AleDocs/CONGRESO/parsimonia/tabla_score.csv", sep = ",")

head(listadoTopologias)

listadoNumeros <- read.csv("/media/aleu/AleDocs/CONGRESO/parsimonia/terminales")

head(listadoNumeros)

todo<- data.frame(listadoTopologias,listadoNumeros)

### R puede manejar cosas comprimidas con gz. asi ahoramos espacio
### por ejemplo pasamos de 1.1 m a 29.2 k

distanciaRFasim <- read.csv("/media/aleu/AleDocs/CONGRESO/parsimonia/intentoasimetrico.gz", sep = ",",header = T)

head(distanciaRFasim)

### Localiza la diagonal 

rfasim <-as.matrix(distanciaRFasim)
rfasim[1,1]

### Remplaza la diagonal por NA

diag(rfasim) <-  NA
rfasim[1,1]

### Histograma de yo(pesos Iguales) vs todos (piwi), con la diagonal de 0 

yovsTodo <- rfasim[which(lower.tri(rfasim, diag = TRUE))]

head(yovsTodo)
tail(yovsTodo)

length(yovsTodo)

newTodos <- yovsTodo[!is.na(yovsTodo)] 

length(newTodos)
if (!is.null(dev.list())){
  dev.off()}

### Parametros de la grafica 
par(mfrow = c(1,2))
hist(yovsTodo)

hist(newTodos)

### FuturoNA: arriba de la diagonal pondremos NA, pero la diagonal es de 0.  

FuturoNA <- which(x = upper.tri(rfasim, diag = F))

rfasim[FuturoNA] <- NA

tail(rfasim)

### Extraer pesos iguales 

listadoIguales <- which(todo$peso=="iguales")

### Comparar pesos Iguales vs Iguales

IgualvsIgual <- rfasim[listadoIguales,listadoIguales]

tail(IgualvsIgual)
hist(IgualvsIgual)
mean(IgualvsIgual,
     na.rm = TRUE)

### Comparar pesos iguales con peso 2 

nivelesPeso <- levels(todo$peso)

contNiveles  <- length(nivelesPeso)

newCont <- contNiveles%/%2+1

##dev.off()

par(mfrow = c(newCont,2))

for (i in nivelesPeso){

pesoTmp <- which(todo$peso==i)                    #Lo extraigo de la tabla
pesAdo <- rfasim[pesoTmp,listadoIguales]          # Extraigo de RF los que son pesos (2, 5, 9...) y pesos iguales

hist(pesAdo, main = i, xlab = "RF Asym")

mean(pesAdo,
     na.rm = TRUE)

}
## Graficos
## 27 de Abril 2017
## El resumen de todos los histogramas (IgualesvsCadaGen) con media. (Boxplot.default para ignorar los NAs.)

par(mfrow = c(2,1))

boxplot.default(rfasim[pesoTmp,listadoIguales],
                rfasim[which(listadoTopologias$peso==2),listadoIguales],
                rfasim[which(listadoTopologias$peso==5),listadoIguales],
                rfasim[which(listadoTopologias$peso==9),listadoIguales],
                rfasim[which(listadoTopologias$peso==18),listadoIguales],
                rfasim[which(listadoTopologias$peso==36),listadoIguales],
                rfasim[which(listadoTopologias$peso==72),listadoIguales],
                pch = 1:10, frame.plot = T,
                main = "Comparación de IgualvsTodos los pesos",
                xlab = "Pesos (desde Igual a 72)", ylab = "(Frecuencia) Media de RF")

## Distribucion general de datos

boxplot(rfasim[which(listadoTopologias$peso==i)], 
        rfasim[which(listadoTopologias$peso==2)],
        rfasim[which(listadoTopologias$peso==5)],
        rfasim[which(listadoTopologias$peso==9)],
        rfasim[which(listadoTopologias$peso==18)],
        rfasim[which(listadoTopologias$peso==36)],
        rfasim[which(listadoTopologias$peso==72)],
        pch = 1:10, frame.plot = T,
        names = c("iguales","2","5","9","18","36","72"),
        main = "Distribucion general de los datos",
        xlab = "Peso", ylab = "(Frecuencia) Media de RF")

## Violin Plots, para la misma grafica, comportamiento general de los datos

install.packages("vioplot")
library(vioplot)
  x1 <- rfasim[which(listadoTopologias$peso==i)]
  x2 <- rfasim[which(listadoTopologias$peso==2)]
  x3 <- rfasim[which(listadoTopologias$peso==5)]
  x4 <- rfasim[which(listadoTopologias$peso==9)]
  x5 <- rfasim[which(listadoTopologias$peso==18)]
  x6 <- rfasim[which(listadoTopologias$peso==36)]
  x7 <- rfasim[which(listadoTopologias$peso==72)] 

# Cambio na por 0
x1[is.na(x1)] <- 0
x2[is.na(x2)] <- 0
x3[is.na(x3)] <- 0
x4[is.na(x4)] <- 0
x5[is.na(x5)] <- 0
x6[is.na(x6)] <- 0
x7[is.na(x7)] <- 0

vioplot(x1,x2,x3,x4,x5,x6,x7,
        names = c("iguales","2","5","9","18","36","72"),
        col= c("grey")
)
title("Distribucion general de los datos")


