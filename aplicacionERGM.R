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

# Importación de la base
baseERGM <- read.csv("bases_datos/base_comercio_ERGM.csv")

# Datos para 2014 y dejar las primeras columnas según los códigos ISO de primeras
baseERGM24 <- baseERGM[baseERGM$year==2024, ]
baseERGM24 <- baseERGM24[, -c(1,4,5,7,8)]

baseERGM24 <- baseERGM24 %>%
  mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

# Crear la red (Con dirección, sin self-loops, matriz tipo desde lista de aristas)
network_comtrade <- as.network(baseERGM24[, c(1,2,3)], directed=TRUE, matrix.type = "adjacency", 
                               ignore.eval = FALSE, loops=FALSE)
# Estadística de la red
summary(network_comtrade)
network.vertex.names(network_comtrade)

colnames(baseERGM24)

#A <- igraph::as_adjacency_matrix(graph = network_comtrade)
#isSymmetric(as.matrix(A))

#Atributos de las aristas
network::set.edge.attribute(x = network_comtrade, attrname = "CDT_vigente",    value = baseERGM24$CDT_vigente)
network::set.edge.attribute(x = network_comtrade, attrname = "norm_ln_dist",  value = baseERGM24$norm_ln_dist)
network::set.edge.attribute(x = network_comtrade, attrname = "contig",    value = baseERGM24$contig)
network::set.edge.attribute(x = network_comtrade, attrname = "comlang_off", value = baseERGM24$comlang_off)
network::set.edge.attribute(x = network_comtrade, attrname = "comrelig",    value = baseERGM24$comrelig)
network::set.edge.attribute(x = network_comtrade, attrname = "fta_wto",  value = baseERGM24$fta_wto)
network::set.vertex.attribute(x = network_comtrade, attrname = "norm_ln_pib_per_capita_dolares_constantes_o",    value = baseERGM24$norm_ln_pib_per_capita_dolares_constantes_o)
network::set.vertex.attribute(x = network_comtrade, attrname = "norm_ln_pib_per_capita_dolares_constantes_d", value = baseERGM24$norm_ln_pib_per_capita_dolares_constantes_d)
network::set.vertex.attribute(x = network_comtrade, attrname = "norm_entry_cost_o",    value = baseERGM24$norm_entry_cost_o)
network::set.vertex.attribute(x = network_comtrade, attrname = "norm_entry_cost_d",  value = baseERGM24$norm_entry_cost_d)

peso_comercio <- baseERGM24[, c(3)]
peso_comercio

#graficación de los nodos
plot(network_comtrade, main = "peso_comercio_total", label = network.vertex.names(network_comtrade))
plot(network_comtrade, edge.cex = peso_comercio, main = "peso_comercio_total")


ergm_paper <- network_comtrade ~ sum + 
  edgecov(network_comtrade, "CDT_vigente") +
  edgecov(network_comtrade, "norm_ln_dist") +
  absdiff("norm_ln_pib_per_capita_dolares_constantes_o") +
  mutual(form="min") + 
  transitiveweights("min", "max", "min")

ergm.fit.paper <- ergm(formula = ergm_paper, response='peso_comercio_total', reference = ~Poisson)
summary(ergm.fit.paper)

#Establecimiento del ERGM (endogenidad y valores exógenos)
ergm_base <- network_comtrade ~ sum +
  edgecov(network_comtrade, "CDT_vigente") + 
  edgecov(network_comtrade, "fta_wto") +    
  edgecov(network_comtrade, "contig") +
  edgecov(network_comtrade, "norm_ln_dist") +
  edgecov(network_comtrade, "comlang_off") +
  edgecov(network_comtrade, "comrelig") +
  nodeocov("norm_ln_pib_per_capita_dolares_constantes_o") + 
  nodeicov("norm_ln_pib_per_capita_dolares_constantes_d") +
  nodeocov("norm_entry_cost_o") +
  nodeicov("norm_entry_cost_d")

ergm.fit.base <- ergm(formula = ergm_base, response='peso_comercio_total', reference = ~Poisson)
summary(ergm.fit.base)


