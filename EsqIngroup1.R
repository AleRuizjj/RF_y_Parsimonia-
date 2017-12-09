# Alejandra Ruiz.
# Resultados de PAR (RF ) 
# Par_Esqumas de ET en ingroup
# Junio 28 DE 2017
# Modif Julio 22

rm(list=ls()) 
#W <- setwd("/CONGRESO/Datos_Rivera&Daza/src/R/Par/")
#Lu <- setwd("/media/aleu/AleDocs/CONGRESO/Datos_Rivera&Daza/src/R/Par/")

# Cargo las matrices 
#setwd("/media/aleu/AleDocs/CONGRESO/Datos_Rivera&Daza/src/R/Par")
RF_asim <- as.matrix(read.csv("EsqIngroup_RFasim.gz"))

ListTopo <- read.csv("ListaMyIngroup_topo.csv", sep = ",")

head(ListTopo)
tail(ListTopo)

#"Eliminar" (Poner NA) por encima de la diagonal 
FuturoNA <- which(x = upper.tri(RF_asim, diag = F))
RF_asim[FuturoNA] <- NA

tail(RF_asim)

# Genes(y esquemas) vs ET
# Genes 12s, 16s, coI, 12coI, 16coI, 1216, ET(12+16+coI)

#####
# Gen 12s vs ET
# Creo una tabla con los siguientes titulos: RFasim, nter, Gen, Peso 

# ArbolIguales son aquellos arboles con peso-igual y tipo-arbol (Consenso NO)
 
ArbolIguales <- ListTopo$peso=="iguales"&ListTopo$tipo=="arbol"

# Así, ArbolP$ son aquellos arboles con peso $ y tipo- arbol (Consenso NO)

ArbolP2 <- ListTopo$peso=="2"&ListTopo$tipo=="arbol"
ArbolP5 <- ListTopo$peso=="5"&ListTopo$tipo=="arbol"
ArbolP9 <- ListTopo$peso=="9"&ListTopo$tipo=="arbol"
ArbolP18 <- ListTopo$peso=="18"&ListTopo$tipo=="arbol"
ArbolP36 <- ListTopo$peso=="36"&ListTopo$tipo=="arbol"
ArbolP72 <- ListTopo$peso=="72"&ListTopo$tipo=="arbol"

# Peso Iguales 12s ####

s12Ig <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolIguales&ListTopo$Gen=="ETde3"), 
                                                   which(ArbolIguales&ListTopo$Gen=="12s")]), rep(139,50)))

s12Ig <- cbind(s12Ig, "12s", "iguales")
colnames(s12Ig) <- c("RFasim","nter", "Gen", "Peso")

head(s12Ig)
str(s12Ig)

# Peso 2 de 12s

s12P2 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP2&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP2&ListTopo$Gen=="12s")]), rep(139,1)))

s12P2 <- cbind(s12P2, "12s", "2")
colnames(s12P2) <- c("RFasim","nter", "Gen", "Peso")

head(s12P2)
str(s12P2)

# Peso 5 de 12s

s12P5 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP5&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP5&ListTopo$Gen=="12s")]), rep(139,1)))

s12P5 <- cbind(s12P5, "12s", "5")
colnames(s12P5) <- c("RFasim","nter", "Gen", "Peso")

head(s12P5)
str(s12P5)

# Peso 9 de 12s

s12P9 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP9&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP9&ListTopo$Gen=="12s")]), rep(139,1)))

s12P9 <- cbind(s12P9, "12s", "9")
colnames(s12P9) <- c("RFasim","nter", "Gen", "Peso")

head(s12P9)
str(s12P9)

# Peso 18 de 12s

s12P18 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP18&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP18&ListTopo$Gen=="12s")]), rep(139,1)))

s12P18 <- cbind(s12P18, "12s", "18")
colnames(s12P18) <- c("RFasim","nter", "Gen", "Peso")

head(s12P18)
str(s12P18)

# Peso 36 de 12s

s12P36 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP36&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP36&ListTopo$Gen=="12s")]), rep(139,1)))

s12P36 <- cbind(s12P36, "12s", "36")
colnames(s12P36) <- c("RFasim","nter", "Gen", "Peso")

head(s12P36)
str(s12P36)

# Peso 72 de 12s

s12P72 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP72&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP72&ListTopo$Gen=="12s")]), rep(139,1)))

s12P72 <- cbind(s12P72, "12s", "72")
colnames(s12P72) <- c("RFasim","nter", "Gen", "Peso")

head(s12P72)
str(s12P72)

# Peso Iguales del gen 16s

s16Ig <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolIguales&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolIguales&ListTopo$Gen=="16s")]), rep(175,80)))

s16Ig <- cbind(s16Ig, "16s", "iguales")
colnames(s16Ig) <- c("RFasim","nter", "Gen", "Peso")

head(s16Ig)
str(s16Ig)

# Peso 2 de 16s

s16P2 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP2&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP2&ListTopo$Gen=="16s")]), rep(175,3)))

s16P2 <- cbind(s16P2, "16s", "2")
colnames(s16P2) <- c("RFasim","nter", "Gen", "Peso")

head(s16P2)
str(s16P2)

# Peso 5 de 16s 

s16P5 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP5&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP5&ListTopo$Gen=="16s")]), rep(175,3)))

s16P5 <- cbind(s16P5, "16s", "5")
colnames(s16P5) <- c("RFasim","nter", "Gen", "Peso")

head(s16P5)
str(s16P5)

# Peso 9 de 16s 

s16P9 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP9&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP9&ListTopo$Gen=="16s")]), rep(175,1)))

s16P9 <- cbind(s16P9, "16s", "9")
colnames(s16P9) <- c("RFasim","nter", "Gen", "Peso")

head(s16P9)
str(s16P9)

# Peso 18 de 16s 

s16P18 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP18&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP18&ListTopo$Gen=="16s")]), rep(175,1)))

s16P18 <- cbind(s16P18, "16s", "18")
colnames(s16P18) <- c("RFasim","nter", "Gen", "Peso")

head(s16P18)
str(s16P18)

# Peso 36 de 16s 

s16P36 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP36&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP36&ListTopo$Gen=="16s")]), rep(175,3)))

s16P36 <- cbind(s16P36, "16s", "36")
colnames(s16P36) <- c("RFasim","nter", "Gen", "Peso")

head(s16P36)
str(s16P36)

# Peso 72 de 16s 

s16P72 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP72&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP72&ListTopo$Gen=="16s")]), rep(175,3)))

s16P72 <- cbind(s16P72, "16s", "72")
colnames(s16P72) <- c("RFasim","nter", "Gen", "Peso")

head(s16P72)
str(s16P72)


# Peso Iguales del gen coI

coIIg <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolIguales&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolIguales&ListTopo$Gen=="coI")]), rep(51,5)))

coIIg <- cbind(coIIg, "coI", "iguales")
colnames(coIIg) <- c("RFasim","nter", "Gen", "Peso")

head(coIIg)
str(coIIg)

# Peso 2 Gen coI

coIP2 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP2&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP2&ListTopo$Gen=="coI")]), rep(51,1)))

coIP2 <- cbind(coIP2, "coI", "2")
colnames(coIP2) <- c("RFasim","nter", "Gen", "Peso")

head(coIP2)
str(coIP2)

# Peso 5 Gen coI

coIP5 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP5&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP5&ListTopo$Gen=="coI")]), rep(51,1)))

coIP5 <- cbind(coIP5, "coI", "5")
colnames(coIP5) <- c("RFasim","nter", "Gen", "Peso")

head(coIP5)
str(coIP5)

# Peso 9 Gen coI

coIP9 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP9&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP9&ListTopo$Gen=="coI")]), rep(51,1)))

coIP9 <- cbind(coIP9, "coI", "9")
colnames(coIP9) <- c("RFasim","nter", "Gen", "Peso")

head(coIP9)
str(coIP9)

# Peso 18 Gen coI

coIP18 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP18&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP18&ListTopo$Gen=="coI")]), rep(51,1)))

coIP18 <- cbind(coIP18, "coI", "18")
colnames(coIP18) <- c("RFasim","nter", "Gen", "Peso")

head(coIP18)
str(coIP18)

# Peso 36 Gen coI

coIP36 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP36&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP36&ListTopo$Gen=="coI")]), rep(51,1)))

coIP36 <- cbind(coIP36, "coI", "36")
colnames(coIP36) <- c("RFasim","nter", "Gen", "Peso")

head(coIP36)
str(coIP36)

# Peso 72 Gen coI

coIP72 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP72&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP72&ListTopo$Gen=="coI")]), rep(51,1)))

coIP72 <- cbind(coIP72, "coI", "72")
colnames(coIP72) <- c("RFasim","nter", "Gen", "Peso")

head(coIP72)
str(coIP72)

# Peso Iguales del gen 12coI

s12coIIg <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolIguales&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolIguales&ListTopo$Gen=="12coI")]), rep(166,40)))

s12coIIg <- cbind(s12coIIg, "s12coI", "iguales")
colnames(s12coIIg) <- c("RFasim","nter", "Gen", "Peso")

head(s12coIIg)
str(s12coIIg)

# Peso 2 Gen s12coI

s12coIP2 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP2&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP2&ListTopo$Gen=="12coI")]), rep(166,1)))

s12coIP2 <- cbind(s12coIP2, "s12coI", "2")
colnames(s12coIP2) <- c("RFasim","nter", "Gen", "Peso")

head(s12coIP2)
str(s12coIP2)

# Peso 5 Gen s12coI

s12coIP5 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP5&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP5&ListTopo$Gen=="12coI")]), rep(166,1)))

s12coIP5 <- cbind(s12coIP5, "s12coI", "5")
colnames(s12coIP5) <- c("RFasim","nter", "Gen", "Peso")

head(s12coIP5)
str(s12coIP5)

# Peso 9 Gen s12coI

s12coIP9 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP9&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP9&ListTopo$Gen=="12coI")]), rep(166,1)))

s12coIP9 <- cbind(s12coIP9, "s12coI", "9")
colnames(s12coIP9) <- c("RFasim","nter", "Gen", "Peso")

head(s12coIP9)
str(s12coIP9)

# Peso 18 Gen s12coI

s12coIP18 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP18&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP18&ListTopo$Gen=="12coI")]), rep(166,1)))

s12coIP18 <- cbind(s12coIP18, "s12coI", "18")
colnames(s12coIP18) <- c("RFasim","nter", "Gen", "Peso")

head(s12coIP18)
str(s12coIP18)

# Peso 36 Gen s12coI

s12coIP36 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP36&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP36&ListTopo$Gen=="12coI")]), rep(166,1)))

s12coIP36 <- cbind(s12coIP36, "s12coI", "36")
colnames(s12coIP36) <- c("RFasim","nter", "Gen", "Peso")

head(s12coIP36)
str(s12coIP36)

# Peso 72 Gen s12coI

s12coIP72 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP72&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP72&ListTopo$Gen=="12coI")]), rep(166,1)))

s12coIP72 <- cbind(s12coIP72, "s12coI", "72")
colnames(s12coIP72) <- c("RFasim","nter", "Gen", "Peso")

head(s12coIP72)
str(s12coIP72)

# Peso Iguales del gen 16coI

s16coIIg <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolIguales&ListTopo$Gen=="ETde3"), 
                                                        which(ArbolIguales&ListTopo$Gen=="16coI")]), rep(179,50)))

s16coIIg <- cbind(s16coIIg, "s16coI", "iguales")
colnames(s16coIIg) <- c("RFasim","nter", "Gen", "Peso")

head(s16coIIg)
str(s16coIIg)

# Peso 2 Gen s16coI

s16coIP2 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP2&ListTopo$Gen=="ETde3"), 
                                                        which(ArbolP2&ListTopo$Gen=="16coI")]), rep(179,1)))

s16coIP2 <- cbind(s16coIP2, "s16coI", "2")
colnames(s16coIP2) <- c("RFasim","nter", "Gen", "Peso")

head(s16coIP2)
str(s16coIP2)

# Peso 5 Gen s16coI

s16coIP5 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP5&ListTopo$Gen=="ETde3"), 
                                                        which(ArbolP5&ListTopo$Gen=="16coI")]), rep(179,1)))

s16coIP5 <- cbind(s16coIP5, "s16coI", "5")
colnames(s16coIP5) <- c("RFasim","nter", "Gen", "Peso")

head(s16coIP5)
str(s16coIP5)

# Peso 9 Gen s16coI

s16coIP9 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP9&ListTopo$Gen=="ETde3"), 
                                                        which(ArbolP9&ListTopo$Gen=="16coI")]), rep(179,1)))

s16coIP9 <- cbind(s16coIP9, "s16coI", "9")
colnames(s16coIP9) <- c("RFasim","nter", "Gen", "Peso")

head(s16coIP9)
str(s16coIP9)

# Peso 18 Gen s16coI

s16coIP18 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP18&ListTopo$Gen=="ETde3"), 
                                                         which(ArbolP18&ListTopo$Gen=="16coI")]), rep(179,1)))

s16coIP18 <- cbind(s16coIP18, "s16coI", "18")
colnames(s16coIP18) <- c("RFasim","nter", "Gen", "Peso")

head(s16coIP18)
str(s16coIP18)

# Peso 36 Gen s16coI

s16coIP36 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP36&ListTopo$Gen=="ETde3"), 
                                                         which(ArbolP36&ListTopo$Gen=="16coI")]), rep(179,1)))

s16coIP36 <- cbind(s16coIP36, "s16coI", "36")
colnames(s16coIP36) <- c("RFasim","nter", "Gen", "Peso")

head(s16coIP36)
str(s16coIP36)

# Peso 72 Gen s16coI

s16coIP72 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP72&ListTopo$Gen=="ETde3"), 
                                                         which(ArbolP72&ListTopo$Gen=="16coI")]), rep(179,1)))

s16coIP72 <- cbind(s16coIP72, "s16coI", "72")
colnames(s16coIP72) <- c("RFasim","nter", "Gen", "Peso")

head(s16coIP72)
str(s16coIP72)


# Peso Iguales del gen 1216 ####

s1216Ig <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolIguales&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolIguales&ListTopo$Gen=="1216")]), rep(179,5)))

s1216Ig <- cbind(s1216Ig, "1216", "iguales")
colnames(s1216Ig) <- c("RFasim","nter", "Gen", "Peso")

head(s1216Ig)
str(s1216Ig)

# Peso 2 de 1216

s1216P2 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP2&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP2&ListTopo$Gen=="1216")]), rep(179,1)))

s1216P2 <- cbind(s1216P2, "1216", "2")
colnames(s1216P2) <- c("RFasim","nter", "Gen", "Peso")

head(s1216P2)
str(s1216P2)

# Peso 5 de 1216 

s1216P5 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP5&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP5&ListTopo$Gen=="1216")]), rep(179,3)))

s1216P5 <- cbind(s1216P5, "1216", "5")
colnames(s1216P5) <- c("RFasim","nter", "Gen", "Peso")

head(s1216P5)
str(s1216P5)

# Peso 9 de 1216 

s1216P9 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP9&ListTopo$Gen=="ETde3"), 
                                                     which(ArbolP9&ListTopo$Gen=="1216")]), rep(179,1)))

s1216P9 <- cbind(s1216P9, "1216", "9")
colnames(s1216P9) <- c("RFasim","nter", "Gen", "Peso")

head(s1216P9)
str(s1216P9)

# Peso 18 de 1216 

s1216P18 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP18&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP18&ListTopo$Gen=="1216")]), rep(179,1)))

s1216P18 <- cbind(s1216P18, "1216", "18")
colnames(s1216P18) <- c("RFasim","nter", "Gen", "Peso")

head(s1216P18)
str(s1216P18)

# Peso 36 de 1216 

s1216P36 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP36&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP36&ListTopo$Gen=="1216")]), rep(179,3)))

s1216P36 <- cbind(s1216P36, "1216", "36")
colnames(s1216P36) <- c("RFasim","nter", "Gen", "Peso")

head(s1216P36)
str(s1216P36)

# Peso 72 de 1216 

s1216P72 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP72&ListTopo$Gen=="ETde3"), 
                                                      which(ArbolP72&ListTopo$Gen=="1216")]), rep(179,3)))

s1216P72 <- cbind(s1216P72, "1216", "72")
colnames(s1216P72) <- c("RFasim","nter", "Gen", "Peso")

head(s1216P72)
str(s1216P72)

#### ENTRE GENES por pesos ####
##   DE 12_16coI ####

s12_16coIIg <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolIguales&ListTopo$Gen=="16coI"), which(ArbolIguales&ListTopo$Gen=="12s")])))

s12_16coIIg <- cbind(s12_16coIIg, 0, as.data.frame("12_16coI"), "iguales")
colnames(s12_16coIIg) <- c("RFasim", "nter","Gen", "Peso")

head(s12_16coIIg)
str(s12_16coIIg)

# Peso 2 de 12_16coI

s12_16coIP2 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP2&ListTopo$Gen=="16coI"), 
                                                      which(ArbolP2&ListTopo$Gen=="12s")])))

s12_16coIP2 <- cbind(s12_16coIP2, 0,"12_16coI", "2")
colnames(s12_16coIP2) <- c("RFasim","nter", "Gen", "Peso")

head(s12_16coIP2)
str(s12_16coIP2)

# Peso 5 de 12_16coI 

s12_16coIP5 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP5&ListTopo$Gen=="16coI"), 
                                                      which(ArbolP5&ListTopo$Gen=="12s")])))

s12_16coIP5 <- cbind(s12_16coIP5, 0,"12_16coI", "5")
colnames(s12_16coIP5) <- c("RFasim","nter", "Gen", "Peso")

head(s12_16coIP5)
str(s12_16coIP5)

# Peso 9 de 12_16coI 

s12_16coIP9 <- cbind(as.numeric(RF_asim[which(ArbolP9&ListTopo$Gen=="16coI"), 
                                        which(ArbolP9&ListTopo$Gen=="12s")]),
                     
                     0, as.data.frame("12_16coI"), "9")

colnames(s12_16coIP9) <- c("RFasim","nter", "Gen", "Peso")

head(s12_16coIP9)
str(s12_16coIP9)

# Peso 18 de 12_16coI

s12_16coIP18 <- cbind(as.numeric(RF_asim[which(ArbolP18&ListTopo$Gen=="16coI"), 
                                         which(ArbolP18&ListTopo$Gen=="12s")]),
                      0, as.data.frame("12_16coI"), "18")

colnames(s12_16coIP18) <- c("RFasim","nter", "Gen", "Peso")

head(s12_16coIP18)
str(s12_16coIP18)

# Peso 36 de 12_16coI

s12_16coIP36 <- cbind(as.numeric(RF_asim[which(ArbolP36&ListTopo$Gen=="16coI"), 
                                         which(ArbolP36&ListTopo$Gen=="12s")]),
                      0, as.data.frame("12_16coI"), "36")

colnames(s12_16coIP36) <- c("RFasim","nter", "Gen", "Peso")

head(s12_16coIP36)
str(s12_16coIP36)

# Peso 72 de 12_16coI 

s12_16coIP72 <- cbind(as.numeric(RF_asim[which(ArbolP72&ListTopo$Gen=="16coI"), 
                                         which(ArbolP72&ListTopo$Gen=="12s")]), 
                      0, as.data.frame("12_16coI"), "72")

colnames(s12_16coIP72) <- c("RFasim","nter", "Gen", "Peso")

head(s12_16coIP72)
str(s12_16coIP72)


##   DE 16_12coI ####

s16_12coIIg <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolIguales&ListTopo$Gen=="16coI"), which(ArbolIguales&ListTopo$Gen=="12s")])))

s16_12coIIg <- cbind(s16_12coIIg, 0, as.data.frame("16_12coI"), "iguales")
colnames(s16_12coIIg) <- c("RFasim", "nter","Gen", "Peso")

head(s16_12coIIg)
str(s16_12coIIg)

# Peso 2 de 16_12coI

s16_12coIP2 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP2&ListTopo$Gen=="12coI"), 
                                                  which(ArbolP2&ListTopo$Gen=="16s")])))

s16_12coIP2 <- cbind(s16_12coIP2, 0,"16_12coI", "2")
colnames(s16_12coIP2) <- c("RFasim","nter", "Gen", "Peso")

head(s16_12coIP2)
str(s16_12coIP2)

# Peso 5 de 16_12coI 

s16_12coIP5 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP5&ListTopo$Gen=="12coI"), 
                                                  which(ArbolP5&ListTopo$Gen=="16s")])))

s16_12coIP5 <- cbind(s16_12coIP5, 0,"16_12coI", "5")
colnames(s16_12coIP5) <- c("RFasim","nter", "Gen", "Peso")

head(s16_12coIP5)
str(s16_12coIP5)

# Peso 9 de 16_12coI 

s16_12coIP9 <- cbind(as.numeric(RF_asim[which(ArbolP9&ListTopo$Gen=="12coI"), 
                                                  which(ArbolP9&ListTopo$Gen=="16s")]),
                     
                     0, as.data.frame("16_12coI"), "9")

colnames(s16_12coIP9) <- c("RFasim","nter", "Gen", "Peso")

head(s16_12coIP9)
str(s16_12coIP9)

# Peso 18 de 16_12coI

s16_12coIP18 <- cbind(as.numeric(RF_asim[which(ArbolP18&ListTopo$Gen=="12coI"), 
                                                   which(ArbolP18&ListTopo$Gen=="16s")]),
                   0, as.data.frame("16_12coI"), "18")

colnames(s16_12coIP18) <- c("RFasim","nter", "Gen", "Peso")

head(s16_12coIP18)
str(s16_12coIP18)

# Peso 36 de 16_12coI

s16_12coIP36 <- cbind(as.numeric(RF_asim[which(ArbolP36&ListTopo$Gen=="12coI"), 
                                                   which(ArbolP36&ListTopo$Gen=="16s")]),
                      0, as.data.frame("16_12coI"), "36")

colnames(s16_12coIP36) <- c("RFasim","nter", "Gen", "Peso")

head(s16_12coIP36)
str(s16_12coIP36)

# Peso 72 de 16_12coI 

s16_12coIP72 <- cbind(as.numeric(RF_asim[which(ArbolP72&ListTopo$Gen=="12coI"), 
                                                   which(ArbolP72&ListTopo$Gen=="16s")]), 
                      0, as.data.frame("16_12coI"), "72")

colnames(s16_12coIP72) <- c("RFasim","nter", "Gen", "Peso")

head(s16_12coIP72)
str(s16_12coIP72)

##   DE coI_1216 ####

coI_1216Ig <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolIguales&ListTopo$Gen=="1216"), which(ArbolIguales&ListTopo$Gen=="coI")])))

coI_1216Ig <- cbind(coI_1216Ig, 0, as.data.frame("coI_1216"), "iguales")
colnames(coI_1216Ig) <- c("RFasim", "nter","Gen", "Peso")

head(coI_1216Ig)
str(coI_1216Ig)

# Peso 2 de coI_1216

coI_1216P2 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP2&ListTopo$Gen=="1216"), 
                                                      which(ArbolP2&ListTopo$Gen=="coI")])))

coI_1216P2 <- cbind(coI_1216P2, 0,"coI_1216", "2")
colnames(coI_1216P2) <- c("RFasim","nter", "Gen", "Peso")

head(coI_1216P2)
str(coI_1216P2)

# Peso 5 de coI_1216 

coI_1216P5 <- as.data.frame(cbind(as.numeric(RF_asim[which(ArbolP5&ListTopo$Gen=="1216"), 
                                                      which(ArbolP5&ListTopo$Gen=="coI")])))

coI_1216P5 <- cbind(coI_1216P5, 0,"coI_1216", "5")
colnames(coI_1216P5) <- c("RFasim","nter", "Gen", "Peso")

head(coI_1216P5)
str(coI_1216P5)

# Peso 9 de coI_1216 

coI_1216P9 <- cbind(as.numeric(RF_asim[which(ArbolP9&ListTopo$Gen=="1216"), 
                                        which(ArbolP9&ListTopo$Gen=="coI")]),
                     
                     0, as.data.frame("coI_1216"), "9")

colnames(coI_1216P9) <- c("RFasim","nter", "Gen", "Peso")

head(coI_1216P9)
str(coI_1216P9)

# Peso 18 de coI_1216

coI_1216P18 <- cbind(as.numeric(RF_asim[which(ArbolP18&ListTopo$Gen=="1216"), 
                                         which(ArbolP18&ListTopo$Gen=="coI")]),
                      0, as.data.frame("coI_1216"), "18")

colnames(coI_1216P18) <- c("RFasim","nter", "Gen", "Peso")

head(coI_1216P18)
str(coI_1216P18)

# Peso 36 de coI_1216

coI_1216P36 <- cbind(as.numeric(RF_asim[which(ArbolP36&ListTopo$Gen=="1216"), 
                                         which(ArbolP36&ListTopo$Gen=="coI")]),
                      0, as.data.frame("coI_1216"), "36")

colnames(coI_1216P36) <- c("RFasim","nter", "Gen", "Peso")

head(coI_1216P36)
str(coI_1216P36)

# Peso 72 de coI_1216 

coI_1216P72 <- cbind(as.numeric(RF_asim[which(ArbolP72&ListTopo$Gen=="1216"), 
                                         which(ArbolP72&ListTopo$Gen=="coI")]), 
                      0, as.data.frame("coI_1216"), "72")

colnames(coI_1216P72) <- c("RFasim","nter", "Gen", "Peso")

head(coI_1216P72)
str(coI_1216P72)


# TABLA A MANO ####

setET <- rbind(s12Ig,
                s12P2,
                s12P5,
                s12P9,
                s12P18,
                s12P36,
                s12P72,
                s16Ig,
                s16P2,
                s16P5,
                s16P9,
                s16P18,
                s16P36,
                s16P72,
                coIIg,
                coIP2,
                coIP5,
                coIP9,
                coIP18,
                coIP36,
                coIP72,
                s12coIIg,
                s12coIP2,
                s12coIP5,
                s12coIP9,
                s12coIP18,
                s12coIP36,
                s12coIP72,
                s16coIIg,
                s16coIP2,
                s16coIP5,
                s16coIP9,
                s16coIP18,
                s16coIP36,
                s12coIP72,
                s1216Ig,
                s1216P2,
                s1216P5,
                s1216P9,
                s1216P18,
                s1216P36,
                s1216P72,
               s12_16coIIg,
               s12_16coIP2,
               s12_16coIP5,
               s12_16coIP9,
               s12_16coIP18,
               s12_16coIP36,
               s12_16coIP72,
               s16_12coIIg,
               s16_12coIP2,
               s16_12coIP5,
               s16_12coIP9,
               s16_12coIP18,
               s16_12coIP36,
               s16_12coIP72,
							 coI_1216Ig,
							 coI_1216P2,
							 coI_1216P5,
							 coI_1216P9,
							 coI_1216P18,
							 coI_1216P36,
							 coI_1216P72
                )
          
head(setET)
str(setET)
tail(setET)
write.csv(x = setET, file = "tablaETbien")

# Correlacion

cor.test(setET$nter, setET$RFasim)
#-0,815


#if (!is.null(dev.list())){
#  dev.off()}

#write.table(setET, file = "setET")


### GRAFICA A MANO ####

library("ggplot2")

color <- ggplot(setET, aes(x= setET$nter, y= setET$RFasim, group=setET$Peso, shape=setET$Gen)) + 
  geom_point(size=2, aes(colour=setET$Peso))  + 
  geom_smooth(se = FALSE, method = "lm", fullrange=TRUE) + 
  scale_shape_manual(values=c(16, 19, 4, 8, 7, 17)) + 
  theme_light() + 
  scale_x_continuous(breaks = pretty(setET$nter, n = 5)) + 
  scale_y_continuous(breaks = pretty(setET$RFasim, n = 9)) + 
  ggtitle("Gráfica Esquemas vs ET en terminos de Numero de Terminales)") + 
  xlab("Numero de terminales") + 
  ylab(" Distancia de RF") + 
  labs(colour = "Valores de K", shape= "Genes") 

##ggsave(filename="myPlot.pdf", plot=final)

### ILD 
# Así se creo la tabla de ILD
# LEN VS RF
# LEN VS ILD
# ILD VS RF

setET_IL <- read.csv("setET_IL", sep = ",", dec = ".")
  

unique(setET_IL$Gen)


AMANO_LenvsRF <- ggplot(setET_IL, aes(setET_IL$length, setET_IL$RFasim, group=setET_IL$Peso, shape=setET_IL$Gen)) + 
  geom_point(size=2, aes(colour=setET_IL$Peso))  + 
  geom_smooth(aes(colour=setET_IL$Peso), se = FALSE, method = "lm", fullrange=TRUE) +
  scale_shape_manual(values=c(5, 19, 4, 8, 7, 17, 3)) +
  theme_light() +
  scale_x_continuous(breaks = pretty(setET_IL$length, n = 9)) +
  scale_y_continuous(breaks = pretty(setET_IL$RFasim, n = 9)) +
  ggtitle("Ingroup: RF vs longitud de pb") + 
  xlab("longitud de pb") + 
  ylab("Distancia de RF") + 
  labs(colour = "Valores de K", shape= "Genes") 

setwd("../../../Fig_test/Figs/")
##ggsave("AMANO_LenvsRF.jpg", AMANO_LenvsRF)

AMANO_LenvsILD <- ggplot(setET_IL, aes(setET_IL$length, setET_IL$stand, group=setET_IL$Peso, shape=setET_IL$Gen)) + 
  geom_point(size=2, aes(colour=setET_IL$Peso))  + 
  geom_smooth(aes(colour=setET_IL$Peso), se = FALSE, method = "lm", fullrange=TRUE) +
  scale_shape_manual(values=c(5, 19, 4, 8, 7, 17, 3)) + 
  theme_light() + 
  scale_x_continuous(breaks = pretty(setET_IL$length, n = 9)) + 
  scale_y_continuous(breaks = pretty(setET_IL$stand, n = 9)) + 
  ggtitle("Ingroup: ILD vs longitud de pb)") + 
  xlab("longitud de pb") + 
  ylab("ILD") + 
  labs(colour = "Valores de K", shape= "Genes") 

##ggsave("AMANO_LenvsILD.jpg", AMANO_LenvsILD)

AMANO_RFvsILD <- ggplot(setET_IL, aes(setET_IL$stand, setET_IL$RFasim, group=setET_IL$Peso, shape=setET_IL$Gen)) + 
  geom_point(size=2, aes(colour=setET_IL$Peso)) + 
  geom_smooth(aes(colour=setET_IL$Peso), se = FALSE, method = "lm", fullrange=TRUE) +
  scale_shape_manual(values=c(5, 19, 4, 8, 7, 17, 3)) + 
  theme_light() + 
  scale_x_continuous(breaks = pretty(setET_IL$stand, n = 9)) + 
  ggtitle("Ingroup: RF vs ILD") + 
  xlab("longitud de pb") + 
  ylab("ILD") + 
  labs(colour = "Valores de K", shape= "Genes") 
#scale_y_continuous(breaks = pretty(setET_IL$RFasim, n = 2))

##ggsave("aMANO_RFvsILD.jpg", AMANO_RFvsILD)

cor.test(setET_IL$stand, setET_IL$RFasim)
#0.1940



### El for  ####

peso <- c("iguales", 2, 5, 9, 18, 36, 72)
InGroTABLA <- matrix(nrow = 0, ncol = 5)

for (g in 1:length(ListTopo$Gen)) {
  #
  for (p in peso) {
  
  mintab <- cbind(as.data.frame(as.character(ListTopo$Gen[g])),
                                as.character(ListTopo[g,2]),
                                as.numeric(RF_asim[which(ListTopo$peso==p&ListTopo$tipo=="arbol"&ListTopo$Gen==(as.character(ListTopo$Gen[g]))), which(ListTopo$peso==p&ListTopo$tipo=="arbol"&ListTopo$Gen=="ETde3")]),
                                as.numeric(ListTopo[g,8]), 
                                rep(ListTopo[which(ListTopo$peso==p&ListTopo$tipo=="arbol"&ListTopo$Gen==as.character(ListTopo$Gen[g])),5,1],5))

  
  InGroTABLA <- rbind(InGroTABLA, mintab)
  
  }
}


colnames(InGroTABLA) <- c("Gen","peso","RF", "length","nterm")
head(InGroTABLA)
tail(InGroTABLA)
str(InGroTABLA)
warning()

#write.table(InGroTABLA, file= "InGrotabla_de_RF", sep=",")
#Agrego a mano el ILD 
setwd("../../src/R/Par/")
InGroTABLAc <- read.csv("tablaETbienRF_ILD.csv")

head(InGroTABLAc)
tail(InGroTABLAc)
cor.test(InGroTABLAc$lenTot, InGroTABLAc$RFasim)
#-0.348
#scale_shape_manual(values=1:nlevels(InGroTABLA$Gen))
#theme_light() 

## Graficas.  con el for####

library("ggplot2")

InGroTABLA_RF <- rbind(InGroTABLAc[(InGroTABLAc$Gen=="12s"),],
                       InGroTABLAc[(InGroTABLAc$Gen=="16s"),],
                       InGroTABLAc[(InGroTABLAc$Gen=="coI"),],
                       InGroTABLAc[(InGroTABLAc$Gen=="12s"),],
                       InGroTABLAc[(InGroTABLAc$Gen=="12coI"),],
                       InGroTABLAc[(InGroTABLAc$Gen=="16coI"),],
                       InGroTABLAc[(InGroTABLAc$Gen=="1216"),],
                       InGroTABLAc[(InGroTABLAc$Gen=="12-16coI"),],
                       InGroTABLAc[(InGroTABLAc$Gen=="16-12coI"),],
                       InGroTABLAc[(InGroTABLAc$Gen=="coI-1216"),]
)

cor.test(InGroTABLA_RF$lenTot, InGroTABLA_RF$RFasim)

Ingroup_RFvsLen <- ggplot(InGroTABLA_RF, aes(InGroTABLA_RF$lenTot, InGroTABLA_RF$RFasim, group=InGroTABLA_RF$Peso, shape=InGroTABLA_RF$Gen)) + 
  geom_point(size=2, aes(colour=InGroTABLA_RF$Peso))  + 
  scale_shape_manual(values=1:nlevels(InGroTABLA_RF$Gen))+ 
  geom_smooth(aes(colour=InGroTABLA_RF$Peso),se = FALSE, method = "lm", fullrange=TRUE) +
  theme_minimal() + 
  scale_x_continuous(breaks = pretty(InGroTABLA_RF$lenTot, n = 9)) + 
  ggtitle("Ingroup: RF vs Número de pares de bases") + 
  xlab("Número de pb") + ylab(" Distancia de RF") + 
  labs(colour = "Valores de K", shape= "Genes") 

#setwd("/media/aleu/AleDocs/CONGRESO/Datos_Rivera&Daza/Fig_test/Figs/")
#ggsave("ingroupRF_Len.jpg", Ingroup_RFvsLen)

head(InGroTABLAc)
tail(InGroTABLAc)

InGroTABLA_ILRF <- rbind(InGroTABLAc[(InGroTABLAc$Gen=="12coI"),],
                         InGroTABLAc[(InGroTABLAc$Gen=="16coI"),],
                         InGroTABLAc[(InGroTABLAc$Gen=="1216"),])
#
                         InGroTABLAc[(InGroTABLAc$Gen=="12-16coI"),]
                         InGroTABLAc[(InGroTABLAc$Gen=="16-12coI"),]
                         InGroTABLAc[(InGroTABLAc$Gen=="coI-1216"),]



RFvsILD <- ggplot(InGroTABLA_ILRF, aes(InGroTABLA_ILRF$stan, InGroTABLA_ILRF$RFasim, group=InGroTABLA_ILRF$Peso, shape=InGroTABLA_ILRF$Gen)) + 
  geom_point(size=2, aes(colour=InGroTABLA_ILRF$Peso)) + 
  scale_shape_manual(values=1:nlevels(InGroTABLA_ILRF$Gen)) + 
  geom_smooth(aes(colour=InGroTABLA_ILRF$Peso), se = FALSE, method = "lm") + 
  theme_minimal() +  
  ggtitle("Ingroup: RF vs ILD") + 
  xlab("ILD") + 
  ylab(" Distancia de RF") + 
  labs(colour = "Valores de K", shape= "Genes") 

cor.test(InGroTABLA_ILRF$stan, InGroTABLA_ILRF$RF)

#0,197
#setwd("/media/aleu/AleDocs/CONGRESO/Datos_Rivera&Daza/Fig_test/Figs/")
ggsave("RFvsILD.jpg", RFvsILD)

ILDvsLen <- ggplot(InGroTABLA, aes(InGroTABLA$length, InGroTABLA$stan, group=InGroTABLA$peso, shape=InGroTABLA$Gen)) + 
  geom_point(size=2, aes(colour=InGroTABLA$peso)) + 
  scale_shape_manual(values=1:nlevels(InGroTABLA$Gen)) + 
  geom_smooth(aes(colour=InGroTABLA$peso), se = FALSE, method = "lm") + 
  theme_minimal() + 
  ggtitle("Ingroup: ILD vs longitud de pb") + 
  xlab("longitud de pb") + ylab(" ILD") + 
  labs(colour = "Valores de K", shape= "Genes") 

#cor.test(InGroTABLA$length, InGroTABLA$stan)
#0.287

#ggsave("ILDvsLen.jpg", ILDvsLen)
