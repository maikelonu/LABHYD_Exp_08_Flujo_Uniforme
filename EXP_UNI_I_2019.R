# //////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# https://www.tec.ac.cr
# Session: FLUJO UNIFORME

# M.Sc. Eng. Maikel Mendez M
# Water Resources + GIS + DataScience
# Instituto Tecnologico de Costa Rica
# https://www.tec.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://twitter.com/MaikelMendezM
# https://github.com/maikelonu
# //////////////////////////////////////////////////////////////////////////////////

# INFO:
# Analisis grafico avanzado
# ggplot2
# lattice
# Normalizacion y homogenizacion de variables
# Exportaci?n ASCII"
# //////////////////////////////////////////////////////////////////////////////////

# Scientific notation is suppress
options(scipen = 0)

# Workspace is cleared
rm(list = ls())

# Working directory is selected
# setwd("/media/maikel/Trabajo/R_ITC/R_LABHYD/EXP_UNI")
setwd("C:/DATOS/R_ITC/R_LABHYD/EXP_UNI")

# CRAN libraries are loaded
require(Agreement)
require(DescTools)
require(effects)
require(ggplot2)
require(MASS)
require(nls2)
require(nlstools)
require(pastecs)
require(reshape)
require(visreg)

# /////////////////////////////////////////////////////////////
# BLOCK: Custom function, round data.frame to specif digits
# /////////////////////////////////////////////////////////////
round_df <- function(df, digits) {
  options(scipen = 0)
  options(scipen = -3)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return(df)
}

# ////////////////////////////////////////////////////////
# BLOCK: Declarations
# ////////////////////////////////////////////////////////
base_m <- 0.086
visco <- 1e-06

# ////////////////////////////////////////////////////////
# BLOCK: Data input
# ////////////////////////////////////////////////////////

# Input data is loaded and a data.frame is created
df.base <- read.table("uniforme.txt", header = TRUE)

# Desc {DescTools} function is requested
Desc(df.base, plotit = TRUE)

# names {base} function is requested
names(df.base)

# a new slope factor variable (as %) is created
df.base$slope_perc <- as.factor(df.base$slope_m_m * 100)

# flow units are converted from m3/h to m3/s
df.base$q_m3_s <- (df.base$q_m3_h) / 3600

# water depth units are converted from cm to m
df.base$y_m <- (df.base$y_cm) / 100

# hydraulic area is calculated
df.base$area <- (df.base$y_m) * base_m

# hydraulic perimeter is calculated
df.base$perimeter <- ((df.base$y_m) * 2) + base_m

# hydraulic radius is calculated
df.base$radius <- (df.base$area / df.base$perimeter)

# square root of hydraulic radius is calculated
df.base$radius.root <- (df.base$radius) ^ 0.5

# water velocity is calculated
df.base$vel <- (df.base$q_m3_s / df.base$area)

# Chezy experimental coefficient is calculated
df.base$c.exp <- df.base$vel / ((df.base$radius * df.base$slope_m_m) ^ 0.5)

# n of Manning is calculated
df.base$n.Manning <- (df.base$radius ^ (1/6)) / df.base$c.exp

# m of Kutter is calculated
df.base$m.Kutter <- (df.base$radius.root * 100 / df.base$c.exp) - df.base$radius.root

# s of Bazin is calculated
df.base$s.Bazin <- ((87 / df.base$c.exp) - 1) * df.base$radius.root

# Froude number is calculated
df.base$Froude <- df.base$vel / ((df.base$area * 9.81 / base_m) ^ 0.5)

# Froude number is rounded to 2 decimals
# df.base$Froude <- round(df.base$Froude, 2)

# coefficients are normalized based on their mean
df.base$n.Manning.norm <- df.base$n.Manning / median(df.base$n.Manning)
df.base$m.Kutter.norm <- df.base$m.Kutter / median(df.base$m.Kutter)
df.base$s.Bazin.norm <- df.base$s.Bazin / median(df.base$s.Bazin)

# stat.desc {pastecs} function is requested and rounded to 3 digits  
df.base.desc <- round(stat.desc(df.base[, 18:20]),3)

# A ggplot object is created
# As velocity is too large we round to 1 decimal just for display
df.base$tempV <- round((df.base$vel), 2)

fg01 <- ggplot(aes(y = y_cm,x = slope_perc),data=df.base) +
  geom_boxplot(aes(colour = slope_perc),size = 0.75) +
  geom_point(aes(colour = slope_perc),size = 3.5) +
  geom_text(aes(label = q_m3_h),size = 6.0,angle = 45,
            hjust = 0.015,vjust = 1.0,parse = FALSE) +
  geom_text(aes(label = tempV),size = 4.0,angle = 0, 
            hjust = 1.0,vjust = 0.0,parse = FALSE, color = "blue") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Boxplot - Distribucion de tirante por pendiente. Etiquetado por Caudal y Velocidad") +
  xlab("Pendiente (%)") +
  ylab("Tirante (cm)") +
  theme_bw(base_size = 18.0)

# A ggplot object is requested
fg01

# We get rid of tempV
df.base$tempV <- NULL

# A ggplot object is created
# As Froude # is too large we round to 1 decimal just for display
df.base$tempF <- round((df.base$Froude), 2)

# As velocity is too large we round to 1 decimal just for display
df.base$tempV <- round((df.base$vel), 2)

fg02 <- ggplot(aes(y = y_cm,x = slope_perc),data=df.base) +
  geom_boxplot(aes(colour = slope_perc),size = 0.75) +
  geom_point(aes(colour = slope_perc),size = 3.5) +
  geom_text(aes(label = tempF),size = 6.0,angle = 45,
            hjust = 0.015,vjust = 1.0,parse = FALSE) +
  geom_text(aes(label = tempV),size = 4.0,angle = 0, 
            hjust = 1.0,vjust = 0.0,parse = FALSE, color = "blue") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Boxplot - Distribucion de tirante por pendiente. Etiquetado por # Froude y Velocidad") +
  xlab("Pendiente (%)") +
  ylab("Tirante (cm)") +
  theme_bw(base_size = 18.0)

# A ggplot object is requested
fg02

# We get rid of tempF
df.base$tempF <- NULL

# We get rid of tempV
df.base$tempV <- NULL

# A ggplot object is created
fg03 <- ggplot(aes(y = n.Manning,x = slope_perc),data=df.base) +
  geom_boxplot(aes(colour = slope_perc),size = 0.75,outlier.size = 0.01) +
  geom_point(aes(colour = slope_perc),data=df.base,size = 3.0) +
  geom_hline(data=df.base,yintercept = 0.0,size = 0.75,linetype = 2,alpha = 0.60) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  geom_text(aes(x = slope_perc,y = n.Manning,label = q_m3_h),data=df.base,
            size = 6.0,angle = 45,hjust = 0.0,vjust = 1.0,parse = FALSE) +
  ggtitle("Boxplot - Distribucion del coeficiente de Manning") +
  xlab("Pendiente (%)") +
  ylab("Coeficiente (-)") +
  theme_bw(base_size = 18.0)

# A ggplot object is requested
fg03

# A ggplot object is created
fg04 <- ggplot(aes(y = m.Kutter,x = slope_perc),data=df.base) +
  geom_boxplot(aes(colour = slope_perc),size = 0.75,outlier.size = 0.01) +
  geom_point(aes(colour = slope_perc),data=df.base,size = 3.0) +
  geom_hline(data=df.base,yintercept = 0.0,size = 0.75,linetype = 2,alpha = 0.60) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Boxplot - Distribucion del coeficiente de Kutter") +
  xlab("Pendiente (%)") +
  ylab("Coeficiente (-)") +
  theme_bw(base_size = 18.0)

# A ggplot object is requested
fg04

# A ggplot object is created
fg05 <- ggplot(aes(y = s.Bazin,x = slope_perc),data=df.base) +
  geom_boxplot(aes(colour = slope_perc),size = 0.75,outlier.size = 0.01) +
  geom_point(aes(colour = slope_perc),data=df.base,size = 3.0) +
  geom_hline(data=df.base,yintercept = 0.0,size = 0.75,linetype = 2,alpha = 0.60) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Boxplot - Distribucion del coeficiente de Bazin") +
  xlab("Pendiente (%)") +
  ylab("Coeficiente (-)") +
  theme_bw(base_size = 18.0)

# A ggplot object is requested
fg05

# A normalized-coefficients data.frame is created
df.norm <- df.base[c("slope_perc",
                     "n.Manning.norm",
                     "m.Kutter.norm",
                     "s.Bazin.norm")]

# melt {reshape} function is requested to convert data 
# from "wide" to "long" format
df.norm <- melt(df.norm)

# A ggplot object is created
fg06 <- ggplot(aes(y = value,x = slope_perc),data=df.norm) +
  geom_boxplot(aes(colour = variable), size = 0.75,outlier.size = 3.5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  geom_hline(data=df.base,yintercept = 0.0,size = 0.75,linetype = 2,alpha = 0.60) +
  ggtitle("Boxplot - Comparacion de coeficientes normalizados (respecto de la mediana)") +
  xlab("Pendiente (%)") +
  ylab("Coeficiente normalizado (-)") +
  theme_bw(base_size = 18.0) 

# A ggplot object is requested
fg06

# round_df function is applied to relevant data.frames
df.output <- round_df(df=df.base, digits=3)
df.base.desc <- round_df(df=df.base.desc, digits=3)

# Objects to export:
# df.output, df.base.desc
# fg01, fg02, fg03, fg06
write.csv(df.output, file = "df.output.csv")
write.csv(df.base.desc, file = "df.base.desc.csv")

# /////////////////////////////////////////////////////////////
# END OF SCRIPT
# /////////////////////////////////////////////////////////////
