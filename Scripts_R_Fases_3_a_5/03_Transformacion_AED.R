#Objetivo: Proyección a grafos, exploración topológica y estadística descriptiva.
#Contenido:
#     Carga del CSV base completa y conversión a objetos network iterativos para cada año.
#     Binarización de la red.
#     Análisis Exploratorio
#         Pruebas de hipótesis volumen comercial
#         Generación del gráfico de la distribución logarítmica del comercio
#     Topología de Red: Extracción de la tabla de evolución descriptiva anual 
#         (Nodos, Aristas, Densidad, Agrupamiento, Reciprocidad, Distancia Media).
#     Salida (Output): Tablas descriptivas y la visualización de la evolución de la red (los grafos por bloques económicos).
#install.packages("countrycode")

library(tidyverse)
library(igraph)
library(graphlayouts)
library(countrycode)

#===========================================
# Carga del CSV base completa y conversión a objetos network iterativos para cada año.
#===========================================

baseERGM <- read.csv("bases_datos/base_comercio_ERGM.csv")


# Extraer vector 'reporteros' (Países con serie completa)
reporteros <- baseERGM %>%
  group_by(iso3_o) %>%
  summarise(n = n_distinct(year)) %>%
  filter(n == 21) %>%
  pull(iso3_o) %>%
  sort()

reporteros
# 2. Carga de datos
baseERGM <- read.csv("bases_datos/base_comercio_ERGM.csv")

# Extraer vector 'reporteros' (Países con serie completa)
reporteros <- baseERGM %>%
  group_by(iso3_o) %>%
  summarise(n = n_distinct(year)) %>%
  filter(n == 21) %>%
  pull(iso3_o) %>%
  sort()

# Datos para corte de EDA 2024
baseERGM24 <- baseERGM %>% filter(year == 2024) %>%
  mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

baseERGM24_limpia <- baseERGM24 %>%
  select(iso3_o, iso3_d, ln_peso_comercio_total, ln_dist, ln_gdp_prod, CDT_vigente, fta_wto, comlang_off) %>%
  mutate(
    CDT_factor = ifelse(CDT_vigente == 1, "Con CDT", "Sin CDT"),
    FTA_factor = ifelse(fta_wto == 1, "Con TLC/OMC", "Sin TLC/OMC")
  )

# 3. Análisis Exploratorio Estadístico (Pruebas de Hipótesis)
cat("\n--- Resumen de Medias por CDT ---\n")
print(baseERGM24_limpia %>% group_by(CDT_factor) %>%
        summarise(N = n(), Media = mean(ln_peso_comercio_total, na.rm = TRUE)))

cat("\n--- Prueba T de Welch ---\n")
print(t.test(ln_peso_comercio_total ~ CDT_factor, data = baseERGM24_limpia, var.equal = FALSE))

cat("\n--- Prueba de Wilcoxon ---\n")
print(wilcox.test(ln_peso_comercio_total ~ CDT_factor, data = baseERGM24_limpia))

# 4. Visualizaciones (GGPlot)
p1 <- ggplot(baseERGM24_limpia, aes(x = ln_peso_comercio_total)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  theme_minimal() + labs(title = "Distribución del Logaritmo del Comercio Total", x = "Log(Comercio)", y = "Frecuencia")

p2 <- ggplot(baseERGM24_limpia, aes(x = CDT_factor, y = ln_peso_comercio_total, fill = CDT_factor)) +
  geom_boxplot(alpha = 0.6) + theme_minimal() + theme(legend.position = "none") +
  labs(title = "Comercio vs CDT", x = "", y = "Log(Comercio)")

p3 <- ggplot(baseERGM24_limpia, aes(x = ln_dist, y = ln_peso_comercio_total)) +
  geom_point(alpha = 0.3, color = "darkred") + geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() + labs(title = "Gravedad: Distancia vs Comercio")

# 5. Cálculo de Evolución Topológica (Tabla de la tesis)
df_redes <- baseERGM %>% filter(peso_comercio_total > 0)

calcular_metricas_red <- function(df_year) {
  g <- graph_from_data_frame(df_year[, c("iso3_o", "iso3_d", "peso_comercio_total")], directed = TRUE)
  data.frame(
    Year = unique(df_year$year),
    Nodos = vcount(g),
    Aristas = ecount(g),
    Densidad = round(edge_density(g), 3),
    Agrupamiento = round(transitivity(g, type = "global"), 3),
    Reciprocidad = round(reciprocity(g), 3),
    Distancia_Media = round(mean_distance(g, directed = TRUE, weights = NA), 3)
  )
}

tabla_topologica <- df_redes %>% group_by(year) %>% group_split() %>% map_dfr(calcular_metricas_red)
write.csv(tabla_topologica, "tabla_topologica_evolucion.csv", row.names = FALSE)

# 6. Red Ego de Colombia (2024 - Top 15 Socios)
socios_top15 <- df_redes %>% filter(year == 2024 & (iso3_o == "COL" | iso3_d == "COL")) %>%
  arrange(desc(peso_comercio_total)) %>% slice(1:15)
codigos_top15 <- setdiff(unique(c(socios_top15$iso3_o, socios_top15$iso3_d)), "COL")

df_ego_col <- df_redes %>% filter(year == 2024 & 
                                    (iso3_o == "COL" | iso3_d == "COL" | (iso3_o %in% codigos_top15 & iso3_d %in% codigos_top15)))

ego_col <- graph_from_data_frame(df_ego_col[, c("iso3_o", "iso3_d")], directed = TRUE)

V(ego_col)$color <- ifelse(V(ego_col)$name == "COL", "tomato", "gray85")
E(ego_col)$color <- ifelse(E(ego_col)[.from("COL") | .to("COL")], "dodgerblue3", adjustcolor("orange", alpha.f = 0.15))

l <- layout_in_circle(ego_col)
l[which(V(ego_col)$name == "COL"), ] <- c(0, 0) # COL al centro

png("ego_network_colombia.png", width = 2000, height = 2000, res = 300)
plot(ego_col, layout = l, vertex.size = ifelse(V(ego_col)$name == "COL", 20, 7), 
     vertex.label.color = "black", edge.arrow.size = 0.1, main = "Red Ego Comercial de Colombia (2024)")
dev.off()


# ==============================================================================
# 7. EVOLUCIÓN ESTRUCTURAL DE LA RED (Grilla Multipanel de la Tesis)
# ==============================================================================

# A. Definición de la grilla de años y grupos económicos
anios_panel <- 2004:2024
anios_fig <- c(2004, 2014, 2020, 2024) 

can   <- c("BOL", "ECU", "PER", "COL")         
asean <- c("BRN", "KHM", "IDN", "LAO", "MYS", "MMR", "PHL", "SGP", "THA", "VNM")
ocde  <- c("AUS","AUT","BEL","CAN","CHL","CZE","DNK","EST","FIN","FRA",
           "DEU","GRC","HUN","ISL","IRL","ISR","ITA","JPN","KOR","LVA",
           "LTU","LUX","MEX","NLD","NZL","NOR","POL","PRT","SVK","SVN",
           "ESP","SWE","CHE","TUR","GBR","USA","CRI")
mercosur <- c("ARG", "BRA", "PRY", "URY")

grupo_de <- function(iso) {
  dplyr::case_when(
    iso == "COL"      ~ "Colombia",
    iso %in% c("CHN","IND")   ~ "Asia emergente",
    iso %in% can      ~ "CAN",
    iso %in% asean    ~ "ASEAN",
    iso %in% ocde     ~ "OCDE",
    iso %in% mercosur ~ "Mercosur",
    TRUE              ~ "Otros"
  )
}

# Verificación de clasificación
cat("\n--- Clasificación de Países por Bloque Económico ---\n")
print(table(grupo_de(reporteros)))

# B. Construcción del Panel Base Balanceado para la Matriz
base_panel <- expand_grid(year = 2004:2024, iso3_o = reporteros, iso3_d = reporteros) %>%
  filter(iso3_o != iso3_d) %>%
  left_join(baseERGM, by = c("year", "iso3_o", "iso3_d")) %>%
  mutate(peso_comercio_total = replace_na(peso_comercio_total, 0))

# C. Justificación de Umbrales (Barrido de densidades)
cat("\n--- Evaluación de Densidades por Umbral de Binarización ---\n")
densidades_umbral <- base_panel %>% group_by(year) %>%
  summarise(d_50M  = mean(peso_comercio_total > 5e7),
            d_100M = mean(peso_comercio_total > 1e8),
            d_500M = mean(peso_comercio_total > 5e8),
            d_1B   = mean(peso_comercio_total > 1e9),
            d_5B   = mean(peso_comercio_total > 5e9))
print(densidades_umbral, n = 21)

# Función para convertir dataframe a matriz de adyacencia
hacer_matriz <- function(df, var) {
  m <- df %>% select(iso3_o, iso3_d, valor = all_of(var)) %>% pivot_wider(names_from = iso3_d, values_from = valor) %>% arrange(iso3_o)
  mat <- as.matrix(m[,-1]); rownames(mat) <- m$iso3_o; diag(mat) <- 0
  return(mat)
}

# Construir la lista de redes binarias (umbral de 100M USD)
UMBRAL <- 1e8
lista_redes <- list()

for (i in seq_along(anios_panel)) {
  bt <- base_panel %>% filter(year == anios_panel[i]) %>% mutate(lazo = ifelse(peso_comercio_total > UMBRAL, 1, 0))
  lista_redes[[i]] <- hacer_matriz(bt, "lazo")
}

# D. Configuración Visual Constante (Layout Fijo)
# Usar la matriz del último año de la grilla (2024) como referencia estructural
g_ref <- graph_from_adjacency_matrix(lista_redes[[which(anios_panel == max(anios_fig))]], mode = "directed")
layout_fijo <- layout_with_stress(g_ref)
rownames(layout_fijo) <- V(g_ref)$name

paleta <- c("Colombia" = "#D64545", 
            "CAN" = "#E89C3F", 
            "Mercosur" = "#E8C33F",
            "Asia emergente" = "#690101", 
            "ASEAN" = "#4A9D5B", 
            "OCDE" = "#4A7BA7", 
            "Otros" = "#B0B0B0")

# E. Función Generadora de Panel
dibujar_panel <- function(anio) {
  i <- which(anios_panel == anio)
  A <- lista_redes[[i]]
  g <- graph_from_adjacency_matrix(A, mode = "directed")
  
  grado <- degree(g, mode = "all")
  V(g)$size <- 3 + 12 * (grado / max(grado))
  
  pesos <- base_panel %>% filter(year == anio) %>% select(iso3_o, iso3_d, peso_comercio_total)
  el <- as_edgelist(g)
  w <- log1p(pesos$peso_comercio_total[match(paste(el[,1], el[,2]), paste(pesos$iso3_o, pesos$iso3_d))])
  E(g)$width <- 0.2 + 1.8 * (w - min(w, na.rm=TRUE)) / diff(range(w, na.rm=TRUE))
  
  es_col <- V(g)$name == "COL"
  V(g)$label <- V(g)$name
  V(g)$label.cex <- ifelse(es_col, 0.8, 0.7)
  V(g)$color <- paleta[grupo_de(V(g)$name)]
  V(g)$frame.color <- ifelse(es_col, "black", "transparent")
  V(g)$label[es_col] <- "COL"
  V(g)$size <- 3 + 8 * (grado / max(grado))  
  
  plot(g,
       layout = layout_fijo[V(g)$name, ],
       vertex.label.dist = 0.8,    
       vertex.label.degree = pi/2,   
       vertex.label.color = "black",
       vertex.label.family = "sans",
       edge.color = adjustcolor("#8FA8BE", alpha.f = 0.15),
       edge.arrow.size = 0,
       edge.curved = 0.1,
       margin = c(0.05, 0.05, 0.05, 0.05),
       main = paste0("(", letters[which(anios_fig == anio)], ") ", anio),
       cex.main = 1.6)
}

# F. Impresión de la Cuadrícula Gráfica
png("figura_redes_comercio.png", width = 3000, height = 3200, res = 300)
layout(matrix(c(1,2, 3,4, 5,5), nrow = 3, byrow = TRUE), heights = c(5, 5, 0.7))
par(mar = c(1, 1, 2.5, 1))
for (a in anios_fig) dibujar_panel(a)
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = names(paleta), pt.bg = paleta,
       pch = 21, pt.cex = 2.2, horiz = TRUE, bty = "n", cex = 1.3)
dev.off()

# G. Documentación de Discontinuidad
solicitados <- sort(unique(baseERGM$iso3_o))
excluidos <- setdiff(solicitados, reporteros)
cat("\n--- Países Excluidos del Panel Estructural ---\n")
print(excluidos)

cat("\n--- Cobertura Temporal de Casos Incompletos ---\n")
baseERGM %>% filter(iso3_o %in% c("ARE","BRN","VNM")) %>% group_by(iso3_o) %>%
  summarise(faltante = paste(setdiff(2004:2024, unique(year)), collapse=", ")) %>% print()