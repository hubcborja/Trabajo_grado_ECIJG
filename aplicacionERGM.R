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
library(btergm)
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
    ln_suma_PIB_capita,
    CDT_vigente,
    fta_wto,
    comlang_off,
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

resumen_medias <- baseERGM24_limpia %>%
  group_by(CDT_factor) %>%
  summarise(
    N = n(),
    Media = mean(ln_peso_comercio_total, na.rm = TRUE),
    Mediana = median(ln_peso_comercio_total, na.rm = TRUE),
    Desv_Est = sd(ln_peso_comercio_total, na.rm = TRUE)
  )

print(resumen_medias)

# 2. Prueba t de Welch para diferencia de medias (Bilateral)
prueba_t <- t.test(ln_peso_comercio_total ~ CDT_factor, 
                   data = baseERGM24_limpia,
                   var.equal = FALSE) # var.equal = FALSE aplica la corrección de Welch

print(prueba_t)

# 3. Prueba de robustez no paramétrica (Mann-Whitney / Wilcoxon)
# Útil si la distribución presenta sesgo o colas pesadas
prueba_wilcox <- wilcox.test(ln_peso_comercio_total ~ CDT_factor, 
                             data = baseERGM24_limpia)

print(prueba_wilcox)

# Histograma del volumen de comercio
ggplot(baseERGM24_limpia, aes(x = ln_peso_comercio_total)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribución del Logaritmo del Comercio Total",
       x = "Log(Peso Comercio Total)",
       y = "Frecuencia")
# Convertir las variables dummy a factor para la gráfica
baseERGM24_limpia <- baseERGM24_limpia %>%
  mutate(CDT_factor = ifelse(CDT_vigente == 1, "Con CDT", "Sin CDT"),
         FTA_factor = ifelse(fta_wto == 1, "Con TLC/OMC", "Sin TLC/OMC"))

# Boxplot para Convenios de Doble Tributación
ggplot(baseERGM24_limpia, aes(x = CDT_factor, y = ln_peso_comercio_total, fill = CDT_factor)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Volumen de Comercio según Vigencia de CDT",
       x = "Estado del Convenio",
       y = "Log(Peso Comercio Total)") +
  theme(legend.position = "none")

# Boxplot para Tratados de Libre Comercio
ggplot(baseERGM24_limpia, aes(x = FTA_factor, y = ln_peso_comercio_total, fill = FTA_factor)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Volumen de Comercio según Tratados (TLC/OMC)",
       x = "Estado del Tratado",
       y = "Log(Peso Comercio Total)") +
  theme(legend.position = "none")

ggplot(baseERGM24_limpia, aes(x = ln_dist, y = ln_peso_comercio_total)) +
  geom_point(alpha = 0.3, color = "darkred") +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Efecto Gravitacional: Distancia vs Comercio",
       x = "Log(Distancia)",
       y = "Log(Peso Comercio Total)")

# Crear grafo dirigido a partir de la lista de aristas
g_comtrade <- graph_from_data_frame(d = baseERGM24_limpia, directed = TRUE)

# Calcular el grado de entrada (quién importa de más socios) y salida (quién exporta a más socios)
in_degree <- degree(g_comtrade, mode = "in")
out_degree <- degree(g_comtrade, mode = "out")

# Top 10 países con más socios de exportación
head(sort(out_degree, decreasing = TRUE), 10)

# Detectar comunidades (clusters) para ver cómo se agrupan los países
# Nota: cluster_walktrap funciona bien con redes densas y dirigidas
comunidades <- cluster_walktrap(g_comtrade)
modularity(comunidades)

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
  edgecov(network_comtrade, "CDT_vigente") +
  edgecov(network_comtrade, "ln_dist") 

ergm.fit.base2 <- ergm(formula = ergm_base2, 
                       response = 'ln_peso_comercio_total',
                       reference = ~Poisson,
                       control = control.ergm(parallel = nucleos, 
                                              parallel.type = "PSOCK",
                                              MCMC.prop = ~triadic + sparse,   
                                              MCMC.samplesize = 8000,          
                                              MCMC.interval = 4096))
summary(ergm.fit.base2)
mcmc.diagnostics(ergm.fit.base2)
# Bondad de ajuste
gof2 <- gof(ergm.fit.base2)
par(mfrow=c(2,2))
plot(gof2)
gof2

# Tercer ERGM base con convenios, distancia entre paises, diferencia entre pib, pertenencia OMC
ergm_base3 <- network_comtrade ~ sum + mutual(form="min") + 
  edgecov(network_comtrade, "CDT_vigente") +
  edgecov(network_comtrade, "ln_dist") +
  edgecov(network_comtrade, "fta_wto") 

ergm.fit_base3 <- ergm(formula = ergm_base3, 
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
  edgecov(network_comtrade, "CDT_vigente") + 
  edgecov(network_comtrade, "fta_wto") +   
  edgecov(network_comtrade, "ln_dist") +
  edgecov(network_comtrade, "ln_suma_PIB_capita")+
  edgecov(network_comtrade, "comlang_off")

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

######################################
# TERGM
# 1. Definir el universo de los 217 países (Ordenados alfabéticamente para evitar cruces)
paises_217 <- sort(unique(c(baseERGM24_limpia$iso3_o, baseERGM24_limpia$iso3_d)))

# 2. Crear el esqueleto masivo (21 años x 217 orígenes x 217 destinos = 984,312 filas)
esqueleto_masivo <- expand_grid(
  year = 2004:2024,
  iso3_o = paises_217,
  iso3_d = paises_217
) %>% 
  filter(iso3_o != iso3_d) # Quitamos la diagonal principal

# 3. Rescatar variables estáticas (distancia y lenguaje) para que no se pierdan
vars_estaticas <- baseERGM %>%
  select(iso3_o, iso3_d, ln_dist, comlang_off) %>%
  filter(!is.na(ln_dist)) %>%
  distinct(iso3_o, iso3_d, .keep_all = TRUE)

# 4. Cruzar el esqueleto masivo con tu base histórica
baseERGM_masiva <- esqueleto_masivo %>%
  left_join(baseERGM %>% select(-ln_dist, -comlang_off), 
            by = c("year", "iso3_o", "iso3_d")) %>%
  left_join(vars_estaticas, by = c("iso3_o", "iso3_d"))

# 5. Rellenar los vacíos (El paso crítico)
baseERGM_masiva <- baseERGM_masiva %>%
  mutate(
    # Si Comtrade no reportó, el comercio es 0
    peso_comercio_total = replace_na(peso_comercio_total, 0),
    # Si no hay registro, asumimos que no hay convenio ni idioma común
    CDT_vigente = replace_na(CDT_vigente, 0),
    fta_wto = replace_na(fta_wto, 0),
    comlang_off = replace_na(comlang_off, 0)
  ) %>%
  # Imputar NAs en variables continuas con la mediana de ese año
  group_by(year) %>%
  mutate(
    ln_suma_PIB_capita = ifelse(is.na(ln_suma_PIB_capita), median(ln_suma_PIB_capita, na.rm = TRUE), ln_suma_PIB_capita),
    ln_dist = ifelse(is.na(ln_dist), median(ln_dist, na.rm = TRUE), ln_dist)
  ) %>%
  ungroup()

# Validar que el panel quedó perfecto (Cada año DEBE tener 46,872 filas exactas)
table(baseERGM_masiva$year)

# 6. Construir las listas temporales para el btergm
anios_panel <- 2004:2024 
lista_redes <- list()
cov_CDT <- list()
cov_fta <- list()
cov_dist <- list()
cov_PIB <- list()
cov_lang <- list()

for (i in seq_along(anios_panel)) {
  anio_actual <- anios_panel[i]
  
  base_temp <- baseERGM_masiva %>% filter(year == anio_actual)
  base_temp$lazo_comercial <- ifelse(base_temp$peso_comercio_total > 0, 1, 0)
  
  # Red
  net_df <- base_temp %>% select(iso3_o, iso3_d, lazo_comercial)
  lista_redes[[i]] <- as.network(net_df, directed=TRUE, matrix.type="edgelist", loops=FALSE)
  
  # Matrices 217x217
  mat_cdt <- as.matrix(pivot_wider(base_temp %>% select(iso3_o, iso3_d, CDT_vigente), names_from = iso3_d, values_from = CDT_vigente)[,-1])
  mat_fta <- as.matrix(pivot_wider(base_temp %>% select(iso3_o, iso3_d, fta_wto), names_from = iso3_d, values_from = fta_wto)[,-1])
  mat_dist <- as.matrix(pivot_wider(base_temp %>% select(iso3_o, iso3_d, ln_dist), names_from = iso3_d, values_from = ln_dist)[,-1])
  mat_pib <- as.matrix(pivot_wider(base_temp %>% select(iso3_o, iso3_d, ln_suma_PIB_capita), names_from = iso3_d, values_from = ln_suma_PIB_capita)[,-1])
  mat_lang <- as.matrix(pivot_wider(base_temp %>% select(iso3_o, iso3_d, comlang_off), names_from = iso3_d, values_from = comlang_off)[,-1])
  
  # Limpieza y Diagonal
  mat_cdt[is.na(mat_cdt) | is.infinite(mat_cdt)] <- 0; diag(mat_cdt) <- 0
  mat_fta[is.na(mat_fta) | is.infinite(mat_fta)] <- 0; diag(mat_fta) <- 0
  mat_lang[is.na(mat_lang) | is.infinite(mat_lang)] <- 0; diag(mat_lang) <- 0
  mat_dist[is.na(mat_dist) | is.infinite(mat_dist)] <- median(mat_dist, na.rm = TRUE); diag(mat_dist) <- 0
  mat_pib[is.na(mat_pib) | is.infinite(mat_pib)] <- median(mat_pib, na.rm = TRUE); diag(mat_pib) <- 0
  
  # Nombrar dimensiones
  dimnames(mat_cdt) <- list(paises_217, paises_217)
  dimnames(mat_fta) <- list(paises_217, paises_217)
  dimnames(mat_dist) <- list(paises_217, paises_217)
  dimnames(mat_pib) <- list(paises_217, paises_217)
  dimnames(mat_lang) <- list(paises_217, paises_217)
  
  cov_CDT[[i]] <- mat_cdt; cov_fta[[i]] <- mat_fta; cov_dist[[i]] <- mat_dist
  cov_PIB[[i]] <- mat_pib; cov_lang[[i]] <- mat_lang
  
  cat("Año", anio_actual, "procesado (Matriz 217x217).\n")
}

# 7. Modelo final
tergm_todas <- lista_redes ~ edges + 
  edgecov(cov_CDT) + edgecov(cov_fta) + edgecov(cov_dist) +
  edgecov(cov_PIB) + edgecov(cov_lang) + memory(type = "autoregression")

fit_btergm <- btergm(
  formula = tergm_todas,
  R = 100,                         
  parallel = "snow",               
  ncpus = nucleos,                 
  verbose = TRUE                   
)

summary(fit_btergm)

