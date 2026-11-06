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
# Carga del CSV base completa y conversión a objetos network iterativos para cada año.
#===========================================

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
    ln_gdp_prod,
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
network.edgecount(network_comtrade)

colnames(baseERGM24)

baseERGM24_limpia <- baseERGM24_limpia %>%
  mutate(CDT_factor = ifelse(CDT_vigente == 1, "Con CDT", "Sin CDT"),
         FTA_factor = ifelse(fta_wto == 1, "Con TLC/OMC", "Sin TLC/OMC"))


#===========================================
# Binarización de la red.
#===========================================


#===========================================
# Análisis Exploratorio
#===========================================

resumen_medias <- baseERGM24_limpia %>%
  group_by(CDT_factor) %>%
  summarise(
    N = n(),
    Media = mean(ln_peso_comercio_total, na.rm = TRUE),
    Mediana = median(ln_peso_comercio_total, na.rm = TRUE),
    Desv_Est = sd(ln_peso_comercio_total, na.rm = TRUE)
  )

print(resumen_medias)

# Prueba t de Welch para diferencia de medias (Bilateral)
prueba_t <- t.test(ln_peso_comercio_total ~ CDT_factor, 
                   data = baseERGM24_limpia,
                   var.equal = FALSE) 

print(prueba_t)

# Prueba de robustez no paramétrica (Mann-Whitney / Wilcoxon)
prueba_wilcox <- wilcox.test(ln_peso_comercio_total ~ CDT_factor, 
                             data = baseERGM24_limpia)

print(prueba_wilcox)

# ── 2. Red 27×27: solo díadas entre reporteros ──
reporteros <- sort(unique(baseERGM24$iso3_o))  


#===========================================
# Topología de Red
#===========================================
# -----------------------------------------------------------------------------
# Años a comparar 
# -----------------------------------------------------------------------------
anios_fig <- c(2004, 2014, 2020, 2024) 


can   <- c("BOL", "ECU", "PER", "COL")         
asean <- c("BRN", "KHM", "IDN", "LAO", "MYS", "MMR",
           "PHL", "SGP", "THA", "VNM")
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

# Verificar que ningún país quede mal clasificado
table(grupo_de(reporteros))
setdiff(reporteros, c(can, asean, ocde, mercosur, "COL"))

# Layout FIJO calculado sobre la red más densa: así los países no se mueven
# entre paneles y la comparación visual es válida.
g_ref <- graph_from_adjacency_matrix(
  as.matrix(lista_redes[[which(anios_panel == max(anios_fig))]]),
  mode = "directed")
layout_fijo <- layout_with_stress(g_ref)
rownames(layout_fijo) <- V(g_ref)$name
# -----------------------------------------------------------------------------
#  Función que dibuja un panel
# -----------------------------------------------------------------------------
paleta <- c("Colombia" = "#D64545",
            "CAN"      = "#E89C3F",
            "Mercosur" = "#E8C33F",
            "Asia emergente" = "#690101",
            "ASEAN"    = "#4A9D5B",
            "OCDE"     = "#4A7BA7",
            "Otros"    = "#B0B0B0")

dibujar_panel <- function(anio, etiquetar_top = 15) {
  i <- which(anios_panel == anio)
  A <- as.matrix(lista_redes[[i]])
  g <- graph_from_adjacency_matrix(A, mode = "directed")
  
  # Tamaño por grado
  grado <- degree(g, mode = "all")
  V(g)$size <- 3 + 12 * (grado / max(grado))
  
  # Grosor por volumen
  pesos <- base_panel %>% filter(year == anio) %>%
    select(iso3_o, iso3_d, peso_comercio_total)
  el <- as_edgelist(g)
  w <- log1p(pesos$peso_comercio_total[
    match(paste(el[,1], el[,2]), paste(pesos$iso3_o, pesos$iso3_d))])
  E(g)$width <- 0.2 + 1.8 * (w - min(w, na.rm=TRUE)) / diff(range(w, na.rm=TRUE))
  
  # Etiquetas de los más centrales
  es_col <- V(g)$name == "COL"
  V(g)$label <- V(g)$name
  V(g)$label.cex <- 0.7
  V(g)$label.cex[es_col] <- 0.8
  
  # Color por bloque 
  grupos <- grupo_de(V(g)$name)
  V(g)$color <- paleta[grupos]
  
  # Colombia resaltada
  V(g)$frame.color <- ifelse(es_col, "black", "transparent")
  V(g)$label[es_col] <- "COL"
  V(g)$size <- 3 + 8 * (grado / max(grado))  
  
  plot(g,
       layout = layout_fijo[V(g)$name, ],
       vertex.label.dist  = 0.8,    
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


# Figura comparativa de tres paneles
---
png("figura_redes_comercio.png", width = 3000, height = 3200, res = 300)

layout(matrix(c(1,2,
                3,4,
                5,5), nrow = 3, byrow = TRUE),
       heights = c(5, 5, 0.7))
par(mar = c(1, 1, 2.5, 1))
for (a in anios_fig) dibujar_panel(a)
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = names(paleta), pt.bg = paleta,
       pch = 21, pt.cex = 2.2, horiz = TRUE, bty = "n", cex = 1.3)
dev.off()


solicitados <- sort(unique(baseERGM$iso3_o))
excluidos <- setdiff(solicitados, reporteros)
excluidos

# Y cuántos años tiene cada uno, para documentar la cobertura:
baseERGM %>% group_by(iso3_o) %>%
  summarise(anios = n_distinct(year)) %>%
  filter(anios < 21) %>% arrange(anios)
baseERGM %>% filter(iso3_o %in% c("ARE","BRN","VNM")) %>%
  group_by(iso3_o) %>%
  summarise(faltante = paste(setdiff(2004:2024, unique(year)), collapse=", "))

# Tabla descriptiva que acompaña la figura 

tabla_descriptiva <- lapply(seq_along(anios_panel), function(i) {
  g <- graph_from_adjacency_matrix(as.matrix(lista_redes[[i]]), mode = "directed")
  data.frame(
    anio        = anios_panel[i],
    nodos       = vcount(g),
    aristas     = ecount(g),
    densidad    = round(edge_density(g), 3),
    agrupamiento= round(transitivity(g, type = "global"), 3),
    reciprocidad= round(reciprocity(g), 3),
    dist_media  = round(mean_distance(g, directed = TRUE), 3)
  )
}) %>% bind_rows()

print(tabla_descriptiva, row.names = FALSE)
write.csv(tabla_descriptiva, "tabla_descriptiva_red.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Evolución de la centralidad de países clave (complemento narrativo)
# -----------------------------------------------------------------------------
paises_foco <- c("COL", "USA", "CHN", "BRA", "ESP", "MEX", "PER", "CHL")

evolucion <- lapply(seq_along(anios_panel), function(i) {
  g <- graph_from_adjacency_matrix(as.matrix(lista_redes[[i]]), mode = "directed")
  data.frame(anio = anios_panel[i],
             pais = paises_foco,
             grado_salida = degree(g, v = paises_foco, mode = "out"),
             grado_entrada = degree(g, v = paises_foco, mode = "in"),
             intermediacion = round(betweenness(g, v = paises_foco, directed = TRUE), 1))
}) %>% bind_rows()

write.csv(evolucion, "evolucion_centralidad.csv", row.names = FALSE)

library(countrycode) 
#install.packages("countrycode")
tabla_paises <- data.frame(
  iso = reporteros,
  pais = countrycode(reporteros, "iso3c", "country.name.es")
)
print(tabla_paises, row.names = FALSE)

n <- nrow(tabla_paises); mitad <- ceiling(n/2)
izq <- tabla_paises[1:mitad, ]
der <- tabla_paises[(mitad+1):n, ]
if (nrow(der) < nrow(izq)) der <- rbind(der, data.frame(iso="", pais=""))
cat(paste0(izq$iso, " & ", izq$pais, " & ", der$iso, " & ", der$pais, " \\\\"),
    sep = "\n")