#_____________________________________________________________________________________________________________________________
#
# Analisis relacionn buffers conservacion cacas con exito de genotipado
# paper de Valentina
#
# by Ines Sanchez-Donoso (2025)
#
# Es esta versión trabajo con réplica de PCR como unidad. Y tengo el coverage de cada réplica.
#
# RStudio 2024.12.0+467 "Kousa Dogwood" Release (cf37a3e5488c937207f992226d255be71f5e3f41, 2024-12-11) for windows
#_____________________________________________________________________________________________________________________________


# 12 cacas * 2 metodos conservacion * 7 replicas PCR * 32 micros

#Esta base de datos no incluye estos loci:
#- 2 del cromosoma sexual: Y650-79.3 y Y990-35
#- 1 del mitocondrial: Mtdll1/2 



#PACKAGES ----

library(ggplot2)
#install.packages("car")
library(car)  #Anova()  #VIF
#install.packages("patchwork")
library(patchwork) # To display 2 charts together

library(nlme) #lme

#install.packages("glmm")
library(glmm) #glmm

library(lme4) #glmm  lmer()
#install.packages("nlme")
library(nlme) #glmm  lme()

#install.packages("reshape2")
library(reshape2)

#install.packages("glmmTMB")
library(glmmTMB) 
citation("glmmTMB")

#install.packages("paletteer")
library(paletteer)

library(gridExtra)
library(grid)


#DATA ----

#setwd("C:/Users/Ines/Dropbox/A_Research/5_Varios/2024_Valentina") #portátil
#setwd("I:/Dropbox/A_Research/5_Varios/2024_Valentina")  #despacho

#despacho ord nuevo
setwd('C:/Users/Usuario/Dropbox/A_Research/5_Varios/2024_Valentina') 


d <- read.table("Base_Replica_Coverage_Genot.txt", header = T)
str(d)


d$caca <- as.factor(d$caca)
d$ID_Field <- as.factor(d$ID_Field)
d$ID_Lab <- as.factor(d$ID_Lab)
d$replica <- as.factor(d$replica)
d$caca_replica <- as.factor(d$caca_replica)
d$metodo <- as.factor(d$metodo)
d$Locus <- as.factor(d$Locus)
d$amplificado <- as.factor(d$amplificado)
d$genotipoCorrecto <- as.factor(d$genotipoCorrecto)
d$tipoincorrecto <- as.factor(d$tipoincorrecto)
d$introducidodeTablaS2 <- as.factor(d$introducidodeTablaS2)

str(d)

#reescalo coverage:
d$s.coverage <- scale(d$coverage)


#para Fig4:
dfig4 <- read.table("LocusGenotCorrecto.txt", header = T)
str(dfig4)

Nlocigenot_caca1Et <- sum(subset(dfig4, dfig4$caca=="caca1" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca1NAP <- sum(subset(dfig4, dfig4$caca=="caca1" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca2Et <- sum(subset(dfig4, dfig4$caca=="caca2" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca2NAP <- sum(subset(dfig4, dfig4$caca=="caca2" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca3Et <- sum(subset(dfig4, dfig4$caca=="caca3" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca3NAP <- sum(subset(dfig4, dfig4$caca=="caca3" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca4Et <- sum(subset(dfig4, dfig4$caca=="caca4" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca4NAP <- sum(subset(dfig4, dfig4$caca=="caca4" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca5Et <- sum(subset(dfig4, dfig4$caca=="caca5" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca5NAP <- sum(subset(dfig4, dfig4$caca=="caca5" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca6Et <- sum(subset(dfig4, dfig4$caca=="caca6" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca6NAP <- sum(subset(dfig4, dfig4$caca=="caca6" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca7Et <- sum(subset(dfig4, dfig4$caca=="caca7" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca7NAP <- sum(subset(dfig4, dfig4$caca=="caca7" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca8Et <- sum(subset(dfig4, dfig4$caca=="caca8" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca8NAP <- sum(subset(dfig4, dfig4$caca=="caca8" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca9Et <- sum(subset(dfig4, dfig4$caca=="caca9" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca9NAP <- sum(subset(dfig4, dfig4$caca=="caca9" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca10Et <- sum(subset(dfig4, dfig4$caca=="caca10" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca10NAP <- sum(subset(dfig4, dfig4$caca=="caca10" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca11Et <- sum(subset(dfig4, dfig4$caca=="caca11" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca11NAP <- sum(subset(dfig4, dfig4$caca=="caca11" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)
Nlocigenot_caca12Et <- sum(subset(dfig4, dfig4$caca=="caca12" & dfig4$metodo=="Ethanol")$locusgenotcorrecto_print)
Nlocigenot_caca12NAP <- sum(subset(dfig4, dfig4$caca=="caca12" & dfig4$metodo=="NAP")$locusgenotcorrecto_print)



dtf <- data.frame(
  NAP = c(dNAP_NPCRsNoAmplif, dNAP_NPCRsCorrectAmplif, dNAP_NPCRs_false_allele, dNAP_NPCRs_dropout, dNAP_NPCRs_ambiguo, dNAP_NPCRs_CoverageInsuf, dNAP_NPCRs_Unscored),
  Etanol = c(dEt_NPCRsNoAmplif, dEt_NPCRsCorrectAmplif, dEt_NPCRs_false_allele, dEt_NPCRs_dropout, dEt_NPCRs_ambiguo, dEt_NPCRs_CoverageInsuf, dEt_NPCRs_Unscored),
  row.names = c("7_NPCRsNoAmplif", "6_NPCRsCorrectAmplif", "5_Nfalse_alleles", "4_Ndropouts", "3_Nambiguos", "2_NCoverage_insuficiente", "1_NUnscored")
)

dtfsumas$category <- row.names(dtfsumas)
longdtfsumas <- melt(dtfsumas, id.vars = "category")

ggplot(longdtfsumas, aes(variable, value, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=value),size=2.5, position=position_stack(0.5)) +
  labs(title= "d", x="method",y= "NAmplifications") +
  theme_bw()




#EXPLORACION ----

#NAP
dNAP <- subset(d, d$metodo=="NAP")  #N PCRs NAP 2688 ok
nrow(dNAP) 
dNAP_NPCRsNoAmplif <- nrow(subset(dNAP, dNAP$amplificado=="no")) #N PCRs NAP no amplificación  1278
dNAP_NPCRsSiAmplif <- nrow(subset(dNAP, dNAP$amplificado=="si")) #N PCRs NAP si amplificación  1410

#De las si amplificación NAP, cuántas no genotipadas:
dNAP_NPCRsNogenotip <- nrow(subset(dNAP, dNAP$amplificado=="si" & dNAP$algungenotipo==0)) #744
  #por coverage insuficiente (0/0):
  dNAP_NPCRsCoverageInsu <-  nrow(subset(dNAP, dNAP$tipoincorrecto=="CoverageInsuficiente")) #732
  #por ruido o statterbands (Unscored/Unscored):
  dNAP_NPCRsUnscored <-  nrow(subset(dNAP, dNAP$genotipo=="Unscored/Unscored")) #12

#De las si amplificación NAP, cuántas si genotipadas:
  #Genotipo ambiguo:
  dNAP_NPCRsAmbiguo <-  nrow(subset(dNAP, dNAP$tipoincorrecto=="ambiguo")) #20
  #Genotipo mal:
    #por drop out:
    dNAP_NPCRsDropout <-  nrow(subset(dNAP, dNAP$tipoincorrecto=="dropout_allele")) #98
    #por false allele:
    dNAP_NPCRsFalso <-  nrow(subset(dNAP, dNAP$tipoincorrecto=="false_allele")) #26
  #Genotipo correcto:
  dNAP_NPCRsCorrecto <-  nrow(subset(dNAP, dNAP$genotipoCorrecto01==1)) #522
    
 
#Ethanol

dEt <- subset(d, d$metodo=="Ethanol")  #N PCRs Ethanol 2688 ok
nrow(dEt) 
dEt_NPCRsNoAmplif <- nrow(subset(dEt, dEt$amplificado=="no")) #N PCRs Etanol no amplificación  1030
dEt_NPCRsSiAmplif <- nrow(subset(dEt, dEt$amplificado=="si")) #N PCRs Etanol si amplificación  1658
  
#De las si amplificación Etanol, cuántas no genotipadas:
dEt_NPCRsNogenotip <- nrow(subset(dEt, dEt$amplificado=="si" & dEt$algungenotipo==0)) #766
  #por coverage insuficiente (0/0):
  dEt_NPCRsCoverageInsu <-  nrow(subset(dEt, dEt$tipoincorrecto=="CoverageInsuficiente")) #748
  #por ruido o statterbands (Unscored/Unscored):
  dEt_NPCRsUnscored <-  nrow(subset(dEt, dEt$genotipo=="Unscored/Unscored")) #18
  
#De las si amplificación Etanol, cuántas si genotipadas con
  #Genotipo ambiguo:
  dEt_NPCRsAmbiguo <-  nrow(subset(dEt, dEt$tipoincorrecto=="ambiguo")) #57
  #Genotipo mal:
    #por drop out:
    dEt_NPCRsDropout <-  nrow(subset(dEt, dEt$tipoincorrecto=="dropout_allele")) #94
    #por false allele:
    dEt_NPCRsFalso <-  nrow(subset(dEt, dEt$tipoincorrecto=="false_allele")) #22
  #Genotipo correcto:
  dEt_NPCRsCorrecto <-  nrow(subset(dEt, dEt$genotipoCorrecto01==1)) #719



dtfsumas <- data.frame(
  NAP = c(dNAP_NPCRsNoAmplif, dNAP_NPCRsCoverageInsu, dNAP_NPCRsUnscored, dNAP_NPCRsAmbiguo, dNAP_NPCRsDropout, dNAP_NPCRsFalso, dNAP_NPCRsCorrecto),
  Etanol = c(dEt_NPCRsNoAmplif, dEt_NPCRsCoverageInsu, dEt_NPCRsUnscored, dEt_NPCRsAmbiguo, dEt_NPCRsDropout, dEt_NPCRsFalso, dEt_NPCRsCorrecto),
  row.names = c("7_NPCRsNoAmplif", "6_NCoverage_insuficiente", "5_NUnscored", "4_Nambiguos", "3_Ndropouts", "2_Nfalse_alleles", "1_NPCRsCorrectGenotype")
)

dtfsumas$category <- row.names(dtfsumas)
longdtfsumas <- melt(dtfsumas, id.vars = "category")

ggplot(longdtfsumas, aes(variable, value, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=value),size=2.5, position=position_stack(0.5)) +
  labs(title= "d", x="method",y= "NAmplifications") +
  theme_bw()



#Figure 3 ----
# Use in a ggplot2 chart:
scale_colour_paletteer_d("rcartocolor::Prism")
scale_fill_paletteer_d("rcartocolor::Prism")


# Definir el orden deseado de los niveles de SAMPLE
longdtfsumas$category <- factor(longdtfsumas$category, levels = c("7_NPCRsNoAmplif", 
                                                                  "6_NCoverage_insuficiente", 
                                                                  "5_NUnscored", 
                                                                  "4_Nambiguos", 
                                                                  "3_Ndropouts", 
                                                                  "2_Nfalse_alleles", 
                                                                  "1_NPCRsCorrectGenotype"))

Etanol <-  ggplot(subset(longdtfsumas, longdtfsumas$variable =="Etanol")) +
  geom_col(aes(x = value, y = category, fill = category), width = 0.6) +
  
  # Etiquetas de los valores al final de las barras
  geom_text(aes(x = value, y = category, label = value),
            position = position_nudge(x = - 80),  # Mueve el texto un poco hacia afuera
            hjust = 0,  # Alineación hacia la derecha
            color = "black", size = 3) +
  scale_x_reverse(limits = c(1300, 0)) +
  labs(title = 'Ethanol', x = 'Number of PCRs') +
    # Personalización de los colores de las categorías
  scale_fill_manual(values = c('#dad7cd', "#E17C05FF","#994E95FF", "#666666FF", '#669bbc', '#CC503EFF', "#73AF48FF")) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),  # Opcional: Oculta el título del eje Y
        axis.text.y = element_text(color = "white"),   # Opcional: Oculta los textos del eje Y
        legend.position = "")

  
NAP <-  ggplot(subset(longdtfsumas, longdtfsumas$variable =="NAP")) +
  geom_col(aes(x = value, y = category, fill = category), width = 0.6) +
  # Etiquetas de los valores al final de las barras
  geom_text(aes(x = value, y = category, label = value),
            position = position_nudge(x = 20),  # Mueve el texto un poco hacia afuera
            hjust = 0,  # Alineación hacia la derecha
            color = "black", size = 3) +
  labs(title = 'NAP', x = 'Number of PCRs') +
  xlim(0, 1300) +
  # Personalización de los colores de las categorías
  scale_fill_manual(values = c('#dad7cd', "#E17C05FF","#994E95FF", "#666666FF", '#669bbc', '#CC503EFF', "#73AF48FF")) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), # Opcional: Oculta los textos del eje Y
        legend.position = "") 

##plot ----
Etanol 
NAP

#Graficar en dos columnas una enfrente de la otra
#Genotipo_clasificado <- grid.arrange(Etanol, NAP, ncol=2)



#Figure 2 ----
#showtext::showtext_auto()

#Relevel caca
dNAP$caca <- factor(dNAP$caca, levels = c("caca1", "caca2", "caca3", "caca4","caca5", "caca6","caca7", "caca8", "caca9", "caca10", "caca11", "caca12"))
dEt$caca <- factor(dEt$caca, levels = c("caca1", "caca2", "caca3", "caca4","caca5", "caca6","caca7", "caca8", "caca9", "caca10", "caca11", "caca12"))


ggplot(dEt, aes(x = coverage, y = Locus, fill = caca)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c('#000000' ,'gold', 'mediumpurple3', '#f58231', 'khaki1', '#46f0f0', 'mediumblue', '#bcf60c','#3cb44b', "darkorchid4",'#aaffc3','#f032e6' )) +
  ggtitle("Ethanol") +
  scale_x_reverse(limits = c(153000, 0)) +
  xlab("") + # Las facetas comparten el mismo eje Y
  theme_classic() +
  theme(axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size=10),
        plot.margin = unit(c(1,-1,1,0), "mm"),
        plot.title = element_text(hjust = 0),
        legend.position = "left")


ggplot(dNAP, aes(x = coverage, y = Locus, fill = caca)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c('#000000' ,'gold', 'mediumpurple3', '#f58231', 'khaki1', '#46f0f0', 'mediumblue', '#bcf60c','#3cb44b', "darkorchid4",'#aaffc3','#f032e6' )) +
  ggtitle("NAP") +
  xlim(0, 153000) +
  xlab("") + # Las facetas comparten el mismo eje Y
  theme_classic() +
  theme(axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size=10),
        plot.margin = unit(c(1,-1,1,0), "mm"),
        plot.title = element_text(hjust = 0),
        legend.position = "left")


#Mean Coverage de Amplificado
dmean <- data.frame(
  Etanol = mean(subset(dEt, dEt$amplificado=="si")$coverage),
  NAP = mean(subset(dNAP, dNAP$amplificado=="si")$coverage),
  row.names = c("meanCoverage"))

dmean$category <- row.names(dmean)
longdmean <- melt(dmean, id.vars = "category")

ggplot(longdmean, aes(variable, value, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=value),size=2.5, position=position_stack(0.5)) +
  labs(title= "Amplificado Si", x="method",y= "Mean Coverage") +
  theme_bw()


#Mean Coverage de Genotipado correcto
dmean2 <- data.frame(
  Etanol = mean(subset(dEt, dEt$genotipoCorrecto=="si")$coverage),
  NAP = mean(subset(dNAP, dNAP$genotipoCorrecto=="si")$coverage),
  row.names = c("meanCoverage"))

dmean2$category <- row.names(dmean2)
longdmean2 <- melt(dmean2, id.vars = "category")

ggplot(longdmean2, aes(variable, value, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=value),size=2.5, position=position_stack(0.5)) +
  labs(title= "Genotipado Correcto", x="method",y= "Mean Coverage") +
  theme_bw()


#Amplificado correcto/Amplificado incorrecto
dmean3 <- data.frame(
  Etanol = c(mean(subset(dEt, dEt$genotipoCorrecto=="si")$coverage), mean(subset(dEt, dEt$genotipoCorrecto=="no")$coverage)),
  NAP = c(mean(subset(dNAP, dNAP$genotipoCorrecto=="si")$coverage), mean(subset(dNAP, dNAP$genotipoCorrecto=="no")$coverage)),
  row.names = c("meanCoveragegenotipoCorrecto", "meanCoverageAmplifNOCorrecto")
)

dmean3$category <- row.names(dmean3)
longdmean3 <- melt(dmean3, id.vars = "category")

ggplot(longdmean3, aes(variable, value, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=value),size=2.5, position=position_stack(0.5)) +
  labs(title= "d", x="method",y= "Mean Coverage") +
  theme_bw()


#No entiendo por qué en los siguientes tests me sale que el coverage afecta a que 
#se ha amplificado bien o no
#pero no a que se haya amplificado algo o no




#TESTS ----


##AS - Amplification success ----

m <- glmmTMB(amplificado01 ~ metodo + s.coverage + (1|caca) + (1|Locus) +(1|caca:Locus),
             data = d, family = nbinom1) #problema convergencia
qqnorm(residuals(m)); qqline(residuals(m))  #regular ajuste
summary(m)


m <- glmmTMB(amplificado01 ~ metodo + s.coverage + (1|caca) + (1|Locus) +(1|caca:Locus),   ###ESTE ----
             data = d, family = nbinom2) #problema convergencia 
qqnorm(residuals(m)); qqline(residuals(m))  #regular ajuste
summary(m)
Anova(m)
str(d)




##GS - Genotyping success ----


m <- glmmTMB(genotipoCorrecto01 ~ metodo + s.coverage + (1|caca) + (1|Locus) +(1|caca:Locus), 
             data = d, family = nbinom1) #problema convergencia
qqnorm(residuals(m)); qqline(residuals(m))  #regular ajuste
summary(m)

#ó:
m <- glmmTMB(genotipoCorrecto01 ~ metodo + s.coverage + (1|caca) + (1|Locus) +(1|caca:Locus), 
             data = d, family = nbinom2) ###ESTE ----
qqnorm(residuals(m)); qqline(residuals(m))  #regular ajuste
summary(m)  
Anova(m)

#los ajustes de los dos modelos (con nbinom1 o nbinom2) son iguales.








#Coverage está muy relacionado con el método
  wilcox.test(coverage ~ as.factor(metodo), data=d)
  
  dEtanol <- subset(d, d$metodo =="Ethanol")
  sum(dEtanol$coverage)
  
  dNAP <- subset(d, d$metodo =="NAP")
  sum(dNAP$coverage)
  
  dtf <- data.frame(
    Coverage = c(sum(dNAP$coverage), sum(dEtanol$coverage)), 
    metodo = c("NAP", "Etanol"))
  
  ggplot(data=dtf)+
    geom_bar(aes(x=metodo, y=Coverage), stat="identity", width = 0.5)

totalcoverage <-  sum(dEtanol$coverage) + sum(dNAP$coverage)

#por lo tanto, no sé si tiene mucho sentido que meta coverage en el modelo.

#Puedo tener dos VR diferentes:
  #1- Amplificado: si/no
  #2- Amplificado correcto: si/no (de los que sí han amplifaicado, cuáles son correctos)
  
#Metiendo coverage:
#y con VR Amplificado correcto: 
#Después de muchas pruebas veo que el mejor modelo es
  #con glmmTMB(), ya que me permite zero-inflation
  #con binomial negativa
  #reescalando coverage
  #con locus y caca como factores random cruzados, y su interacción.
  
  
 
  
#Pruebas anteriores ----
  
  
  
m <- glmer(amplificado ~ metodo + coverage + (1|Locus) + (1|caca), data = d, family = binomial) #no convergencia
plot(m)
qqnorm(residuals(m)); qqline(residuals(m))  #Mal ajuste
summary(m)



m <- glmer(amplificado ~ metodo + s.coverage + (1|Locus) + (1|caca), data = d, family = binomial)
plot(m)
qqnorm(residuals(m)); qqline(residuals(m))  #Mal ajuste
summary(m)


m <- glmmTMB(amplificado ~ metodo + coverage + (1|Locus) + (1|caca), #problema convergencia
                         data=d, family=binomial)

m <- glmmTMB(amplificado ~ metodo + s.coverage + (1|Locus) + (1|caca), #problema convergencia
             data=d, family=binomial)

#puede que est? dando fallo de convergencia porque hay muchos 0 ("no") en "amplificado"??





#VR Amplificado Correcto si/no
  m <- glmer(genotipoCorrecto ~ metodo + coverage + (1|Locus) + (1|caca), 
             data = subset(d, d$amplificado=="si"), family = binomial)  #pide reescalar
  plot(m)
  qqnorm(residuals(m)); qqline(residuals(m))  #Mal ajuste
  summary(m)
  
    #con coverage reescalado:
    m <- glmer(genotipoCorrecto ~ metodo + s.coverage + (1|Locus) + (1|caca), 
               data = subset(d, d$amplificado=="si"), family = binomial) #lo traga bien
    plot(m)
    qqnorm(residuals(m)); qqline(residuals(m))  #Ajuste regular
    summary(m)

  #pruebo con glmmTMB para models with various extensions, including zero-inflation
  m <- glmmTMB(genotipoCorrecto ~ metodo + coverage + (1|Locus) + (1|caca), 
             data = subset(d, d$amplificado=="si"), family = binomial) #lo traga bien pero
  qqnorm(residuals(m)); qqline(residuals(m))  #Muy mal ajuste
  summary(m)

    #reescalo coverage:
    m <- glmmTMB(genotipoCorrecto ~ metodo + s.coverage + (1|Locus) + (1|caca), 
                 data = subset(d, d$amplificado=="si"), family = binomial) #lo traga bien pero
    qqnorm(residuals(m)); qqline(residuals(m))  #Muy mal ajuste
    summary(m)

    
  #pruebo con glmmTMB y family = nbinom2, pongo VR como 0 y 1 (lo necesita):
  m <- glmmTMB(genotipoCorrecto01 ~ metodo + coverage + (1|Locus) + (1|caca), 
             data = subset(d, d$amplificado=="si"), family = nbinom2) #lo traga bien
  qqnorm(residuals(m)); qqline(residuals(m)) #mal ajuste
  summary(m)


  #m <- glmmTMB(as.factor(genotipoCorrecto) ~ 
  #               as.factor(metodo) + coverage + (1|fcaca:fLocus), 
  #             data = subset(d, d$amplificado=="si"), family = binomial)
  #qqnorm(residuals(m)); qqline(residuals(m))  #MUY MAL ajuste
  #summary(m)



#locus no est? nested en caca, son crossed 

#(1|block)+(1|year)  #crossed
#(1|year/block) #block nested en year
#(1|year)+(1|block)+(1|year:block) #year, block and their interaction
 str(d)


m <- glmmTMB(genotipoCorrecto01 ~ metodo + s.coverage + (1|caca) + (1|Locus) +(1|caca:Locus), 
             data = subset(d, d$amplificado=="si"), family = nbinom1) #se lo traga
qqnorm(residuals(m)); qqline(residuals(m))  #regular ajuste
summary(m)


m <- glmmTMB(genotipoCorrecto01 ~ metodo + s.coverage + (1|caca) + (1|Locus) +(1|caca:Locus), 
             data = subset(d, d$amplificado=="si"), family = nbinom2)  ##ESTE
qqnorm(residuals(m)); qqline(residuals(m)) #no tan mal
summary(m)




#por lo tanto, no creo que est? muy bien incluir coverage en el modelo, no?
#Eliminando coverage del modelo:

m <- glmmTMB(genotipoCorrecto01 ~ metodo + (1|caca) + (1|Locus), 
             data = subset(d, d$amplificado=="si"), family = nbinom2)  #problema de convergencia
qqnorm(residuals(m)); qqline(residuals(m))
summary(m)
#problema de convergencia

m <- glmmTMB(genotipoCorrecto01 ~ metodo + (1|caca) + (1|Locus) +(1|caca:Locus), 
             data = subset(d, d$amplificado=="si"), family = nbinom2) ## o ESTE
qqnorm(residuals(m)); qqline(residuals(m)) #no tan mal
summary(m)
#NO problema de convergencia !!!

m <- glmmTMB(genotipoCorrecto01 ~ metodo + (1|caca) + (1|Locus), 
             data = subset(d, d$amplificado=="si"), family = nbinom1)  #problema de convergencia
qqnorm(residuals(m)); qqline(residuals(m))
summary(m)


m <- glmmTMB(genotipoCorrecto01 ~ metodo + (1|caca) + (1|Locus) +(1|caca:Locus), 
             data = subset(d, d$amplificado=="si"), family = nbinom1) #problema de convergencia
qqnorm(residuals(m)); qqline(residuals(m)) #no tan mal
summary(m)


#Por probar, quito m?todo y dejo coverage, para ver c?mo de correlacionados est?n:
m <- glmmTMB(genotipoCorrecto01 ~ s.coverage + (1|fcaca) + (1|fLocus) +(1|fcaca:fLocus), 
             data = subset(d, d$amplificado=="si"), family = nbinom2)
qqnorm(residuals(m)); qqline(residuals(m)) #no tan mal
summary(m)
#Problema de convergencia

m <- glmmTMB(genotipoCorrecto01 ~ s.coverage + (1|fcaca) + (1|fLocus) +(1|fcaca:fLocus), 
             data = subset(d, d$amplificado=="si"), family = nbinom1)
qqnorm(residuals(m)); qqline(residuals(m)) #no tan mal
summary(m)
#Problema de convergencia


#Tiene sentido biol?gico poner la interacci?n de caca y locus?
#Hay locus que funcionan mejor con algunas cacas que con otras? Supongo que puede ser.







#Con glmer si le pongo negativa binomial me sale que tengo que darle un valor de theta (creo), 
   #y yo ese valor no s? cu?l debe ser.
#Con glmer me sale que debo reescalar variables, por lo tanto, usar glmmTMB, que no tiene ese problema
#nbinom2 (y nbinom1) con locus anidado dentro de caca se ajusta bien, pero tiene problema de convergencia, que debo ver c?mo resolver
#Negativas bionomiales necesitan variable "y" num?rica.
#escalando coverage me desaparece el problema del warninig alertando de falsa convergencia

#https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf
#https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html
#Models with convergence problems should be excluded from further consideration, in general.
#In some cases, extreme eigenvalues may be caused by having predictor variables that are on very different scales: try rescaling, and centering, continuous predictors in the model.


#https://search.r-project.org/CRAN/refmans/glmmTMB/html/nbinom2.html
#family glmmTMB
#nbinom2: Negative binomial distribution: quadratic parameterization (Hardin & Hilbe 2007). 
#nbinom1: Negative binomial distribution: linear parameterization (Hardin & Hilbe 2007)

#https://www.nature.com/articles/nmeth.3137
#Con nested locus dentro de cada me sale bien, pero no son nested, son crossed.
#Si los hago crossed o sigo los hago crossed + interacci?n me salen problemas de convergencia.