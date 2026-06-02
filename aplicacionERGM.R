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

baseERGM <- read.csv("bases_datos/base_comercio_ERGM.csv")


baseERGM24 <- baseERGM[baseERGM$year==2024, ]
baseERGM24 <- baseERGM24[, -c(1,4,5,7,8)]
#com_network <- baseERGM24[, c("iso3_o", "iso3_d", "peso_comercio_total")]
#com_network <- com_network[com_network$iso3_o!=com_network$iso3_d,]

network_comtrade <- as.network(baseERGM24[, 1:3], directed=TRUE, matrix.type = "edgelist", ignore.eval = FALSE,
                               loops=FALSE)
summary(network_comtrade)
network.vertex.names(network_comtrade)

colnames(baseERGM24)

A <- igraph::as_adjacency_matrix(graph = network_comtrade)
isSymmetric(as.matrix(A))

network::set.vertex.attribute(x = network_comtrade, attrname = "CDT_vigente",    value = baseERGM24$CDT_vigente)
network::set.vertex.attribute(x = network_comtrade, attrname = "norm_ln_dist",  value = baseERGM24$norm_ln_dist)
network::set.vertex.attribute(x = network_comtrade, attrname = "contig",    value = baseERGM24$contig)
network::set.vertex.attribute(x = network_comtrade, attrname = "comlang_off", value = baseERGM24$comlang_off)
network::set.vertex.attribute(x = network_comtrade, attrname = "comrelig",    value = baseERGM24$comrelig)
network::set.vertex.attribute(x = network_comtrade, attrname = "fta_wto",  value = baseERGM24$fta_wto)
network::set.vertex.attribute(x = network_comtrade, attrname = "norm_ln_pib_per_capita_dolares_constantes_o",    value = baseERGM24$norm_ln_pib_per_capita_dolares_constantes_o)
network::set.vertex.attribute(x = network_comtrade, attrname = "norm_ln_pib_per_capita_dolares_constantes_d", value = baseERGM24$norm_ln_pib_per_capita_dolares_constantes_d)
network::set.vertex.attribute(x = network_comtrade, attrname = "norm_entry_cost_o",    value = baseERGM24$norm_entry_cost_o)
network::set.vertex.attribute(x = network_comtrade, attrname = "norm_entry_cost_d",  value = baseERGM24$norm_entry_cost_d)

peso_comercio <- baseERGM24 %v% 'peso_comercio_total'
peso_comercio

plot(network_comtrade, main = "peso_comercio_total", label = network.vertex.names(network_comtrade))
plot(network_comtrade, vertex.cex = peso_comercio, main = "peso_comercio_total")

attr_comercio <- network_comtrade %v% network_comtrade$peso_comercio_total


ergm_aristas_triangulos <- network_comtrade~edges+ triangles
summary(ergm_aristas_triangulos)
ergm.fit.01 <- ergm(formula = ergm_aristas_triangulos)
summary(ergm.fit.01)



ergm_aristas_triangulos <- network_comtrade~edges+ match("CDT_vigente")+ 
  match("contig")+
  nodemain("norm_ln_dist")+
  nodemain("comlang_off")+
  nodemain("comrelig")+
  match("fta_wto")+
  nodemain("norm_ln_pib_per_capita_dolares_constantes_o")+
  nodemain("norm_ln_pib_per_capita_dolares_constantes_d")+
  nodemain("norm_entry_cost_o")+
  nodemain("norm_entry_cost_d")

summary(ergm_aristas_triangulos)
ergm.fit.01 <- ergm(formula = ergm_aristas_triangulos)
summary(ergm.fit.01)

network_comtrade <- as.network(baseERGM24[, 1:2], 
                               directed = TRUE, 
                               matrix.type = "edgelist")


set.edge.value(network_comtrade, "comercio_valor", baseERGM24[, 3])


summary(network_comtrade, verbose = FALSE)
nombres_nodos <- network.vertex.names(network_comtrade)
nombres_nodos

atributos_nodos <- baseERGM24 %>%
  select(iso3_o, norm_ln_pib_per_capita_dolares_constantes_o, CDT_vigente) %>%
  distinct(iso3_o, .keep_all = TRUE)
atributos_nodos
# Mapeamos los atributos a la red[cite: 1, 2]
set.vertex.attribute(network_comtrade, "pib_capita", atributos_nodos$norm_ln_pib_per_capita_dolares_constantes_o)
set.vertex.attribute(network_comtrade, "CDT", atributos_nodos$CDT_vigente)


modelo_test <- ergm(network_comtrade ~ sum + mutual, 
                    response = "comercio_valor", 
                    reference = ~Poisson)

summary(modelo_test)
