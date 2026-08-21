# Liberias
#install.packages("tidyverse")
#install.packages("igraph")
#install.packages("ergm.count")
#install.packages("ergm.rank")
#install.packages("ergm.ego")
#install.packages("latentnet")

library(tidyverse)
library(network)
library(ergm)
library(Rglpk)
library(igraph)
library(ergm.count)
library(ergm.rank)
library(latentnet)
library(parallel)

nucleos <- detectCores() - 1

# Importación de la base
baseERGM <- read.csv("bases_datos/base_comercio_ERGM.csv")

# Datos para 2014 y dejar las primeras columnas según los códigos ISO de primeras
baseERGM24 <- baseERGM[baseERGM$year==2024, ]

baseERGM24 <- baseERGM24 %>%
  mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))
# Selección de las variables de interés para no tener correlación alta
baseERGM24_limpia <- baseERGM24 %>%
  select(
    iso3_o,
    iso3_d,
    ln_peso_comercio_total,
    ln_dist,
    ln_pib_per_capita_dolares_constantes_o,
    ln_pib_per_capita_dolares_constantes_d,
    CDT_vigente,
    fta_wto,
    comlang_off,
    comrelig,
    entry_cost_o,
    entry_cost_d
  )
baseERGM24_limpia$ln_peso_comercio_total <- round(baseERGM24_limpia$ln_peso_comercio_total)
# Crear la red (Con dirección, sin self-loops, matriz tipo desde lista de aristas)
network_comtrade <- as.network(baseERGM24_limpia, directed=TRUE, matrix.type = "edgelist", 
                               ignore.eval = FALSE, loops=FALSE)
# Estadística de la red
summary(network_comtrade)
network.vertex.names(network_comtrade)

colnames(baseERGM24)



#graficación de los nodos
plot(network_comtrade, main = "peso_comercio_total", label = network.vertex.names(network_comtrade))
plot(network_comtrade, edge.lwd = baseERGM24_limpia$ln_peso_comercio_total, main = "peso_comercio_total")


#ERGM base con solo los convenios y las estadísticas endógenas
ergm_base <- network_comtrade ~ sum + mutual(form="min") + 
  edgecov(network_comtrade, "CDT_vigente") 
#transitiveweights("min", "max", "min") + 

ergm.fit_base <- ergm(formula = ergm_base, 
                      response='ln_peso_comercio_total', 
                      reference = ~Poisson,
                      control = control.ergm(parallel = nucleos, 
                                             parallel.type = "PSOCK",
                                             MCMC.prop = ~triadic + sparse,   
                                             MCMC.samplesize = 8000,          
                                             MCMC.interval = 4096))
summary(ergm.fit_base)
mcmc.diagnostics(ergm.fit_base)
# Bondad de ajuste
gof1 <- gof(ergm.fit_base)
par(mfrow=c(2,2))
plot(gof1)
gof1


# Segundo ERGM base con convenios, distancia entre paises y las estadísticas endógenas
ergm_base2 <- network_comtrade ~ sum + mutual(form="min") + 
  edgecov("CDT_vigente") +
  edgecov("ln_dist") 

ergm.fit.base2 <- ergm(formula = ergm_base2, 
                       response = 'ln_peso_comercio_total',
                       reference = ~Poisson,
                       control = control.ergm(parallel = nucleos, 
                                              parallel.type = "PSOCK",
                                              MCMC.prop = ~triadic + sparse,   
                                              MCMC.samplesize = 8000,          
                                              MCMC.interval = 4096))
summary(ergm.fit_base2)
mcmc.diagnostics(ergm.fit_base2)
# Bondad de ajuste
gof2 <- gof(ergm.fit_base2)
par(mfrow=c(2,2))
plot(gof2)
gof2

# Tercer ERGM base con convenios, distancia entre paises, diferencia entre pib, pertenencia OMC
ergm_base3 <- network_comtrade ~ sum + mutual(form="min") + 
  edgecov("CDT_vigente") +
  edgecov("ln_dist") +
  edgecov("fta_wto") 

ergm.fit.base3 <- ergm(formula = ergm_base3, 
                       response='ln_peso_comercio_total', 
                       reference = ~Poisson,
                       control = control.ergm(parallel = nucleos, 
                                              parallel.type = "PSOCK",
                                              MCMC.prop = ~triadic + sparse,   
                                              MCMC.samplesize = 8000,          
                                              MCMC.interval = 4096))

summary(ergm.fit_base3)
mcmc.diagnostics(ergm.fit_base3)
# Bondad de ajuste
gof3 <- gof(ergm.fit_base3)
par(mfrow=c(2,2))
plot(gof3)
gof3

# Cuarto ERGM con todas las variables exógenas y las endógenas
ergm_todas <- network_comtrade ~ sum + mutual(form="min") + 
  edgecov("CDT_vigente") + 
  edgecov("fta_wto") +   
  edgecov("ln_dist") +
  edgecov("ln_suma_PIB_capita")
  edgecov("comlang_off") 

ergm.fit.todas <- ergm(formula = ergm_todas, 
                       response='ln_peso_comercio_total', 
                       reference = ~Poisson,
                       control = control.ergm(parallel = nucleos, 
                                              parallel.type = "PSOCK",
                                              MCMC.prop = ~triadic + sparse,   
                                              MCMC.samplesize = 8000,          
                                              MCMC.interval = 4096))

summary(ergm.fit.todas)
mcmc.diagnostics(ergm.fit.todas)
# Bondad de ajuste
gof4 <- gof(ergm.fit.todas)
par(mfrow=c(2,2))
plot(gof4)
gof4