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
options(scipen = 999)
library(tidyverse)
library(igraph)
library(graphlayouts)
library(countrycode)

#===========================================
# Carga del CSV base completa 
#===========================================

baseERGM <- read.csv("bases_datos/base_comercio_ERGM.csv")

# 100 millones USD
UMBRAL_TERGM <- 1e8 
anios_panel <- 2004:2024

#===========================================
# verificaciones base de datos
#===========================================

# conteo por año para verificar uniformidad
baseERGM %>% count(year)

# logaritmo de exportaciones
baseERGM$ln_X <- log1p(baseERGM$X)

# Extraer vector 'reporteros' (Países con serie completa)
reporteros <- baseERGM %>%
  group_by(iso3_o) %>%
  summarise(n = n_distinct(year)) %>%
  filter(n == 21) %>%
  pull(iso3_o) %>%
  sort()

reporteros

# Datos para corte de EDA 2024
baseERGM24 <- baseERGM %>% filter(year == 2024) %>%
  mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

baseERGM24_limpia <- baseERGM24 %>%
  mutate(
    CDT_factor = ifelse(CDT_vigente == 1, "Con CDT", "Sin CDT"),
    FTA_factor = ifelse(fta_wto == 1, "Con TLC/OMC", "Sin TLC/OMC")
  )

cat("\n--- Cuantiles de Exportaciones (X) en 2024 ---\n")

print(summary(baseERGM24_limpia$X))
print(quantile(baseERGM24_limpia$X, probs = c(0.25, 0.45, 0.50, 0.54, 0.75, 0.90)))

cat("Porcentaje de díadas que superan 100M USD en 2024 (Densidad):", mean(baseERGM24_limpia$X > 1e8), "\n")

cat("\n--- Resumen de Medias y Medianas de ln_X por CDT (2024) ---\n")
print(
  baseERGM24_limpia %>%
    group_by(CDT_factor) %>%
    summarise(
      N = n(),
      Media_lnX = mean(ln_X, na.rm = TRUE),
      DE_lnX = sd(ln_X, na.rm = TRUE),
      Mediana_lnX = median(ln_X, na.rm = TRUE),
      IQR_lnX = IQR(ln_X, na.rm = TRUE)
    )
)

cat("\n--- Prueba T de Welch (sobre ln_X) ---\n")
print(t.test(ln_X ~ CDT_factor, data = baseERGM24_limpia, var.equal = FALSE))

cat("\n--- Prueba de Wilcoxon (sobre ln_X) ---\n")
print(wilcox.test(ln_X ~ CDT_factor, data = baseERGM24_limpia))


# 4. Visualizaciones (GGPlot)
p1 <- ggplot(baseERGM24_limpia, aes(x = ln_X)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = log1p(UMBRAL_TERGM), color = "red", linetype = "dashed", linewidth = 1) +
  theme_minimal() +
  labs(x = "Log(Exportaciones + 1)", y = "Frecuencia")
print(p1)

p2 <- ggplot(baseERGM24_limpia, aes(x = CDT_factor, y = ln_X, fill = CDT_factor)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "", y = "Log(Exportaciones + 1)")
print(p2)

p3 <- ggplot(baseERGM24_limpia, aes(x = ln_dist, y = ln_X)) +
  geom_point(alpha = 0.3, color = "darkred") +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  theme_minimal() +
  labs(x = "Log(Distancia)", y = "Log(Exportaciones + 1)")
print(p3)

# 5. Cálculo de Evolución Topológica 
umbrales <- c("50M" = 5e7, "100M" = 1e8, "500M" = 5e8, "1000M" = 1e9, "5000M" = 5e9)

barrido_X <- sapply(umbrales, function(u) {
  sapply(anios_panel, function(y) {
    df_t <- baseERGM %>% filter(year == y)
    mean(df_t$X > u)
  })
})
rownames(barrido_X) <- anios_panel
cat("\n--- Barrido de Densidades sobre Exportaciones (X) ---\n")
print(round(barrido_X, 3))

# Función de matriz indexada
hacer_matriz <- function(df, var, nodos = reporteros) {
  mat <- matrix(0, length(nodos), length(nodos), dimnames = list(nodos, nodos))
  mat[cbind(match(df$iso3_o, nodos), match(df$iso3_d, nodos))] <- df[[var]]
  diag(mat) <- 0
  stopifnot(!anyNA(mat))
  return(mat)
}

tabla_topologica <- map_dfr(anios_panel, function(y) {
  df_t <- baseERGM %>%
    filter(year == y) %>%
    mutate(lazo_binario = ifelse(X > UMBRAL_TERGM, 1, 0))
  
  A <- hacer_matriz(df_t, "lazo_binario")
  g <- graph_from_adjacency_matrix(A, mode = "directed")
  
  tibble(
    Year = y,
    Nodos = vcount(g),
    Aristas = ecount(g),
    Densidad = round(edge_density(g), 3),
    Agrupamiento = round(transitivity(g, type = "global"), 3),
    Reciprocidad = round(reciprocity(g), 3),
    Distancia_Media = round(mean_distance(g, directed = TRUE, weights = NA), 3)
  )
})

cat("\n--- Tabla Topológica (Red Dirigida X > 100M USD) ---\n")
print(as.data.frame(tabla_topologica))
write.csv(tabla_topologica, "tabla_topologica_evolucion.csv", row.names = FALSE)

# 6. Red Ego de Colombia (2024 - Top 15 Socios)
df_redes_100M <- baseERGM %>% filter(year == 2024 & X > UMBRAL_TERGM)

# 1. Identificar exactamente los 15 socios principales sumando el flujo bilateral X con COL
socios_top15 <- df_redes_100M %>%
  filter(iso3_o == "COL" | iso3_d == "COL") %>%
  mutate(socio = ifelse(iso3_o == "COL", iso3_d, iso3_o)) %>%
  group_by(socio) %>%
  summarise(flujo_max = max(X), .groups = "drop") %>%
  arrange(desc(flujo_max)) %>%
  slice(1:15) %>%
  pull(socio)

nodos_ego <- c("COL", socios_top15)

# 2. Filtrar las aristas donde AMBOS nodos pertenezcan a la red ego (COL + Top 15)
df_ego_col <- df_redes_100M %>%
  filter(iso3_o %in% nodos_ego & iso3_d %in% nodos_ego)

ego_col <- graph_from_data_frame(df_ego_col[, c("iso3_o", "iso3_d")], directed = TRUE, vertices = nodos_ego)

# 3. Asignar colores sin usar ifelse() sobre secuencias de igraph (evita Invalid indexing)
V(ego_col)$color <- ifelse(V(ego_col)$name == "COL", "tomato", "gray85")

# Primero pintamos todas las aristas entre terceros (naranja tenue) y luego sobrescribimos las de COL (azul)
E(ego_col)$color <- adjustcolor("orange", alpha.f = 0.25)
E(ego_col)[.from("COL") | .to("COL")]$color <- "dodgerblue3"

# 4. Layout circular con Colombia en el centro exacto (0, 0)
l <- layout_in_circle(ego_col)
idx_col <- which(V(ego_col)$name == "COL")
nodos_resto <- setdiff(seq_len(vcount(ego_col)), idx_col)
angulos <- seq(0, 2 * pi, length.out = length(nodos_resto) + 1)[-(length(nodos_resto) + 1)]
l[nodos_resto, ] <- cbind(cos(angulos), sin(angulos))
l[idx_col, ] <- c(0, 0)

png("ego_network_colombia.png", width = 2000, height = 2000, res = 300)
plot(
  ego_col,
  layout = l,
  vertex.size = ifelse(V(ego_col)$name == "COL", 20, 9),
  vertex.label.color = "black",
  vertex.label.cex = 0.85,
  edge.arrow.size = 0.35,
  main = "Red Ego de Exportaciones de Colombia (2024)"
)
dev.off()
# ==============================================================================
# 7. EVOLUCIÓN ESTRUCTURAL DE LA RED 
# ==============================================================================
install.packages("ggraph")
library(ggraph) # o graphlayouts para layout_with_stress

# A. Definición de la grilla de años y grupos económicos (59 países activos)
anios_panel <- 2004:2024
anios_fig   <- c(2004, 2014, 2020, 2024) 

can      <- c("BOL", "ECU", "PER", "COL")         
asean    <- c("BRN", "KHM", "IDN", "MYS", "PHL", "SGP", "THA", "VNM") # Sin LAO ni MMR
ocde     <- c("AUS","AUT","BEL","CAN","CHL","CZE","DNK","EST","FIN","FRA",
              "DEU","GRC","HUN","ISL","IRL","ISR","ITA","JPN","KOR","LVA",
              "LTU","LUX","MEX","NLD","NZL","NOR","POL","PRT","SVK","SVN",
              "ESP","SWE","CHE","TUR","GBR","USA","CRI")
mercosur <- c("ARG", "BRA", "PRY", "URY")

grupo_de <- function(iso) {
  dplyr::case_when(
    iso == "COL"             ~ "Colombia",
    iso %in% c("CHN", "IND") ~ "Asia emergente",
    iso %in% can             ~ "CAN",
    iso %in% asean           ~ "ASEAN",
    iso %in% ocde            ~ "OCDE",
    iso %in% mercosur        ~ "Mercosur",
    TRUE                     ~ "Otros" # Agrupa a ARE, GUY, PAN, SUR
  )
}

# Verificación de clasificación (suma exacta de 59 economías)
cat("\n--- Clasificación de Países por Bloque Económico ---\n")
print(table(grupo_de(reporteros)))

# B. Panel Base Balanceado (baseERGM ya viene balanceada desde Python con 71.862 filas)
base_panel <- baseERGM %>%
  mutate(
    X = replace_na(X, 0),
    peso_comercio_total = replace_na(peso_comercio_total, 0)
  )

stopifnot(nrow(base_panel) == 59 * 58 * 21)

# Construir la lista de redes binarias dirigidas (umbral de 100M USD sobre X)
lista_redes <- list()

for (i in seq_along(anios_panel)) {
  bt <- base_panel %>%
    filter(year == anios_panel[i]) %>%
    mutate(lazo = ifelse(X > UMBRAL_TERGM, 1, 0))
  lista_redes[[i]] <- hacer_matriz(bt, "lazo")
}

# D. Configuración Visual Constante (Layout Fijo con semilla reproducible)
set.seed(2026)
g_ref <- graph_from_adjacency_matrix(lista_redes[[which(anios_panel == max(anios_fig))]], mode = "directed")
layout_fijo <- graphlayouts::layout_with_stress(g_ref)
rownames(layout_fijo) <- V(g_ref)$name

paleta <- c(
  "Colombia"       = "#D64545", 
  "CAN"            = "#E89C3F", 
  "Mercosur"       = "#E8C33F",
  "Asia emergente" = "#690101", 
  "ASEAN"          = "#4A9D5B", 
  "OCDE"           = "#4A7BA7", 
  "Otros"          = "#B0B0B0"
)

# E. Función Generadora de Panel (Ponderada por Exportaciones X)
dibujar_panel <- function(anio) {
  i <- which(anios_panel == anio)
  A <- lista_redes[[i]]
  g <- graph_from_adjacency_matrix(A, mode = "directed")
  
  grado <- degree(g, mode = "all")
  V(g)$size <- 3 + 8 * (grado / max(grado))
  
  pesos <- base_panel %>% filter(year == anio) %>% select(iso3_o, iso3_d, X)
  el <- as_edgelist(g)
  w <- log1p(pesos$X[match(paste(el[, 1], el[, 2]), paste(pesos$iso3_o, pesos$iso3_d))])
  E(g)$width <- 0.2 + 1.8 * (w - min(w, na.rm = TRUE)) / diff(range(w, na.rm = TRUE))
  
  es_col <- V(g)$name == "COL"
  V(g)$label <- V(g)$name
  V(g)$label.cex <- ifelse(es_col, 0.85, 0.7)
  V(g)$color <- paleta[grupo_de(V(g)$name)]
  V(g)$frame.color <- ifelse(es_col, "black", "transparent")
  
  plot(
    g,
    layout = layout_fijo[V(g)$name, ],
    vertex.label.dist = 0.8,    
    vertex.label.degree = pi / 2,   
    vertex.label.color = "black",
    vertex.label.family = "sans",
    edge.color = adjustcolor("#8FA8BE", alpha.f = 0.18),
    edge.arrow.size = 0.15,
    edge.curved = 0.1,
    margin = c(0.05, 0.05, 0.05, 0.05),
    main = paste0("(", letters[which(anios_fig == anio)], ") ", anio),
    cex.main = 1.6
  )
}

# F. Impresión de la Cuadrícula Gráfica
png("figura_redes_comercio.png", width = 3000, height = 3200, res = 300)
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(5, 5, 0.7))
par(mar = c(1, 1, 2.5, 1))
for (a in anios_fig) dibujar_panel(a)
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center", legend = names(paleta), pt.bg = paleta,
  pch = 21, pt.cex = 2.2, horiz = TRUE, bty = "n", cex = 1.3
)
dev.off()

# G. Verificación de Completitud (59 activos vs. 5 excluidos en preprocesamiento)
excluidos_preprocesamiento <- c("CUB", "LAO", "MMR", "TLS", "VEN")
cat("\n--- Países Excluidos por Discontinuidad Estructural (>5 años) ---\n")
print(excluidos_preprocesamiento)

cat("\n--- Verificación de Cobertura Completa en Casos Imputados (t-1) ---\n")
base_panel %>%
  filter(iso3_o %in% c("ARE", "BRN", "VNM")) %>%
  group_by(iso3_o) %>%
  summarise(
    anios_activos = n_distinct(year),
    faltante = paste(setdiff(2004:2024, unique(year)), collapse = ", ")
  ) %>%
  print()
