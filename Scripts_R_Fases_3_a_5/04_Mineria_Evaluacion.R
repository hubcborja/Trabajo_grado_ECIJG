#Objetivo: Estimación de modelos, pruebas de robustez y ranking predictivo (Objetivos 3 y 4 del KDD).
#Contenido:
#     TERGM (Minería): Estimación por bloques anidados (Modelo 1 Base, Modelo 2 Preferencial, Modelo 3 Triádico) 
#           usando el paquete btergm y el bootstrap temporal (500 réplicas).
#     Robustez: Reestimación de los modelos cambiando las ventanas temporales (2015-2024, 2020-2024) y variando 
#           el umbral de binarización (50M y 500M).
#     Evaluación (GOF): Simulación de las 100 redes y generación de los gráficos de bondad de ajuste 
#           (grados, socios compartidos, distancias geodésicas, curvas ROC y PR).
#     Predicción e Inversión del Modelo: El script del modelo logit temporal para predecir la creación del CDT 
#           (con rezagos $t-1$ y filtrando el riesgo, como corregimos recientemente).
#     Salida (Output): El ranking final (Top 15) de países candidatos a negociar un tratado con Colombia, 
#           listo para entregar a la DIAN.

library(tidyverse)
library(network)
library(ergm)
library(Rglpk)
library(igraph)
library(ergm.count)
library(ergm.rank)
library(latentnet)
library(parallel)
library(texreg)

#===========================================
# Importacion y preparación para red desde TERGM y VERGM
#===========================================

baseERGM <- read.csv("bases_datos/base_comercio_ERGM.csv")

#===========================================
# TERGM
#===========================================


#===========================================
# Robustez
#===========================================

# VERGM


#===========================================
# Evaluación
#===========================================


#===========================================
# Salida
#===========================================