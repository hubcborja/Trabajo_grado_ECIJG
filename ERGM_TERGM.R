
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

nucleos <- detectCores() - 1

hacer_matriz <- function(df, variable, relleno = 0) {
  m <- df %>% select(iso3_o, iso3_d, valor = all_of(variable)) %>%
    pivot_wider(names_from = iso3_d, values_from = valor) %>% arrange(iso3_o)
  filas <- m$iso3_o; m <- as.matrix(m[, -1]); rownames(m) <- filas
  m <- m[reporteros, reporteros]
  m[is.na(m) | is.infinite(m)] <- relleno
  diag(m) <- 0
  m
}

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



# Al inicio, una sola vez para todo el script
rep_2024 <- sort(unique(baseERGM24$iso3_o))
rep_tergm <- baseERGM %>% group_by(iso3_o) %>%
  summarise(n = n_distinct(year)) %>% filter(n == 21) %>% pull(iso3_o) %>% sort()

length(rep_2024); length(rep_tergm)
setdiff(rep_2024, rep_tergm)

reporteros <- baseERGM %>%
  group_by(iso3_o) %>% summarise(n_anios = n_distinct(year)) %>%
  filter(n_anios == 21) %>% pull(iso3_o) %>% sort()

base27 <- baseERGM24 %>%
  filter(iso3_o %in% reporteros, iso3_d %in% reporteros) %>%   # <- también filtra el origen
  mutate(ln_peso = round(ln_peso_comercio_total))

nrow(base27)                          # ¿cuántas filas hay?
length(reporteros) * (length(reporteros) - 1)   # ¿cuántas debería haber?
sum(mat_cdt)                          # compara con el conteo del panel del TERGM
sum(cov_CDT[[21]])  

# Matriz valuada 27×27 con dimnames (adyacencia, no edgelist)

b24 <- base_panel %>% filter(year == 2024) %>%
  mutate(ln_peso = round(ln_peso_comercio_total))

ggplot(b24, aes(x = ln_peso_comercio_total)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(x = "Log(Peso Comercio Total)",
       y = "Frecuencia")

mat_peso <- hacer_matriz(b24, "ln_peso")
mat_cdt  <- hacer_matriz(b24, "CDT_vigente")
mat_fta  <- hacer_matriz(b24, "fta_wto")
mat_dist <- hacer_matriz(b24, "ln_dist")
mat_pib  <- hacer_matriz(b24, "ln_gdp_prod")
mat_lang <- hacer_matriz(b24, "comlang_off")

red27 <- as.network(mat_peso, directed = TRUE, matrix.type = "adjacency",
                    ignore.eval = FALSE, names.eval = "ln_peso", loops = FALSE)

# Verificaciones
sum(mat_cdt) == sum(cov_CDT[[21]])    # debe dar TRUE
dim(mat_peso)     

# ── 3. Modelo por bloques — 
# Bloque A: estructura endógena pura
f_A <- red27 ~ sum + mutual(form = "min") +
  nodeicovar(transform = "sqrt") +
  nodeocovar(transform = "sqrt")

fit_A <- ergm(f_A, response = "ln_peso", reference = ~Poisson,
              control = control.ergm(MCMLE.termination = "Hummel",
                                     MCMLE.maxit = 60,
                                     MCMC.samplesize = 16000,
                                     MCMC.interval = 2048))

# Bloque B: + dependencia triádica (completa el bloque endógeno)
f_B <- update(f_A, . ~ . + transitiveweights("min","max","min"))

fit_B <- ergm(f_B, response = "ln_peso", reference = ~Poisson,
              control = control.ergm(MCMLE.termination = "Hummel",
                                     MCMLE.maxit = 60,
                                     MCMC.samplesize = 8000,
                                     MCMC.interval = 2048))

# Bloque C: + controles exógenos (gravedad e instituciones)
f_C <- update(f_B, . ~ . +
                edgecov(mat_pib,  attrname = "PIB")  +
                edgecov(mat_dist, attrname = "DIST") +
                edgecov(mat_lang, attrname = "LANG") +
                edgecov(mat_fta,  attrname = "FTA")  +
                edgecov(mat_cdt,  attrname = "CDT"))

fit_C <- ergm(f_C, response = "ln_peso", reference = ~Poisson,
              control = control.ergm(MCMLE.termination = "Hummel",
                                     MCMLE.maxit = 100,
                                     MCMC.samplesize = 8000,
                                     MCMC.interval = 2048))

screenreg(list(fit_A, fit_B, fit_C),
          custom.model.names = c("Endógeno", "+ Triádico", "+ Controles"))


######################################
# TERGM

library(btergm)
library(fastglm)

# Si al limpiar aparecen díadas duplicadas, consolidar sumando el comercio:
dups <- baseERGM %>% count(year, iso3_o, iso3_d) %>% filter(n > 1)
if (nrow(dups) > 0) {
  baseERGM <- baseERGM %>%
    group_by(year, iso3_o, iso3_d) %>%
    summarise(peso_comercio_total = sum(peso_comercio_total, na.rm = TRUE),
              across(everything(), first), .groups = "drop")
}


reporteros <- baseERGM %>%
  group_by(iso3_o) %>% summarise(n_anios = n_distinct(year)) %>%
  filter(n_anios == 21) %>% pull(iso3_o) %>% sort()

anios_panel <- 2004:2024
cat("Reporteros con serie completa:", length(reporteros), "\n")

# ── 2. Esqueleto balanceado 27x27 x 21 años 
esqueleto <- expand_grid(year = anios_panel,
                         iso3_o = reporteros,
                         iso3_d = reporteros) %>%
  filter(iso3_o != iso3_d)

vars_estaticas <- baseERGM %>%
  select(iso3_o, iso3_d, ln_dist, comlang_off) %>%
  filter(!is.na(ln_dist)) %>%
  distinct(iso3_o, iso3_d, .keep_all = TRUE)

base_panel <- esqueleto %>%
  left_join(baseERGM %>% select(-ln_dist, -comlang_off),
            by = c("year", "iso3_o", "iso3_d")) %>%
  left_join(vars_estaticas, by = c("iso3_o", "iso3_d")) %>%
  mutate(peso_comercio_total = replace_na(peso_comercio_total, 0),
         CDT_vigente = replace_na(CDT_vigente, 0),
         fta_wto     = replace_na(fta_wto, 0),
         comlang_off = replace_na(comlang_off, 0)) %>%
  group_by(year) %>%
  mutate(ln_gdp_prod = ifelse(is.na(ln_gdp_prod),
                                     median(ln_gdp_prod, na.rm = TRUE),
                                     ln_gdp_prod),
         ln_dist = ifelse(is.na(ln_dist), median(ln_dist, na.rm = TRUE), ln_dist)) %>%
  ungroup()

stopifnot(all(table(base_panel$year) == length(reporteros)*(length(reporteros)-1)))




base_panel %>% group_by(year) %>%
  summarise(d_50M  = mean(peso_comercio_total > 5e7),
            d_100M = mean(peso_comercio_total > 1e8),
            d_500M = mean(peso_comercio_total > 5e8),
            d_1B   = mean(peso_comercio_total > 1e9),
            d_5B   = mean(peso_comercio_total > 5e9)) %>% print(n = 21)

UMBRAL <- 1e8

# ── 4. Construcción de listas temporales ─────────────────────────────────────


lista_redes <- list(); cov_CDT <- list(); cov_fta <- list(); cov_PIB <- list()

for (i in seq_along(anios_panel)) {
  bt <- base_panel %>% filter(year == anios_panel[i]) %>%
    mutate(lazo = ifelse(peso_comercio_total > UMBRAL, 1, 0))
  
  lista_redes[[i]] <- as.network(hacer_matriz(bt, "lazo"), directed = TRUE,
                                 matrix.type = "adjacency", loops = FALSE)
  cov_CDT[[i]] <- hacer_matriz(bt, "CDT_vigente")
  cov_fta[[i]] <- hacer_matriz(bt, "fta_wto")
  cov_PIB[[i]] <- hacer_matriz(bt, "ln_gdp_prod",
                               relleno = median(bt$ln_gdp_prod, na.rm = TRUE))
  cat("Año", anios_panel[i], "- densidad:",
      round(network.density(lista_redes[[i]]), 3), "\n")
}

b24 <- base_panel %>% filter(year == 2024)
cov_dist <- hacer_matriz(b24, "ln_dist", relleno = median(b24$ln_dist, na.rm = TRUE))
cov_lang <- hacer_matriz(b24, "comlang_off")

# ── 5. Modelos por bloques (espejo del ERGM: A -> B -> C) ────────────────────
# Bloque A: endógenas básicas + covariables + memoria

f_A <- lista_redes ~ edges + mutual +
  edgecov(cov_CDT) + edgecov(cov_fta) +
  edgecov(cov_dist) + edgecov(cov_PIB) + edgecov(cov_lang) +
  memory(type = "stability", lag = 1) +
  delrecip(mutuality = FALSE, lag = 1)

fit_A <- btergm(f_A, R = 500, parallel = "snow", ncpus = 4,
                usefastglm = TRUE, verbose = TRUE)
summary(fit_A)

# Bloque B: + attachment preferencial (espejo de nodeicovar/nodeocovar)
f_B <- lista_redes ~ edges + mutual +
  gwidegree(0.5, fixed = TRUE) + gwodegree(0.5, fixed = TRUE) +
  edgecov(cov_CDT) + edgecov(cov_fta) +
  edgecov(cov_dist) + edgecov(cov_PIB) + edgecov(cov_lang) +
  memory(type = "stability", lag = 1) +
  delrecip(mutuality = FALSE, lag = 1)

fit_B <- btergm(f_B, R = 500, parallel = "snow", ncpus = 4,
                usefastglm = TRUE, verbose = TRUE)
summary(fit_B)

# Bloque C: + triádicos (espejo de transitiveweights/cyclicalweights)
f_C <- lista_redes ~ edges + mutual +
  gwesp(0.69, fixed = TRUE) + ctriple +
  gwidegree(0.5, fixed = TRUE) + gwodegree(0.5, fixed = TRUE) +
  edgecov(cov_CDT) + edgecov(cov_fta) +
  edgecov(cov_dist) + edgecov(cov_PIB) + edgecov(cov_lang) +
  memory(type = "stability", lag = 1) +
  delrecip(mutuality = FALSE, lag = 1)

fit_C <- btergm(f_C, R = 500, parallel = "snow", ncpus = 4,
                usefastglm = TRUE, verbose = TRUE)
summary(fit_C)

# ── 6. Corrida final y GOF (solo con la especificación ya cerrada) ───────────
fit_final <- btergm(f_C, R = 500, parallel = "snow", ncpus = 4,
                    usefastglm = TRUE, verbose = TRUE)
summary(fit_final)

gof_final <- gof(fit_final, nsim = 100,
                 statistics = c(esp, dsp, geodesic, ideg, odeg,
                                triad.directed, rocpr))
plot(gof_final)

screenreg(list(fit_A, fit_B, fit_C))


idx <- which(anios_panel >= 2015)
redes_10a <- lista_redes[idx]
cdt_10a   <- cov_CDT[idx]
fta_10a   <- cov_fta[idx]
pib_10a   <- cov_PIB[idx]
# cov_dist y cov_lang son matrices únicas: btergm las recicla, no se recortan

f_10a <- redes_10a ~ edges + mutual +
  gwesp(0.69, fixed = TRUE) + ctriple +
  gwidegree(0.5, fixed = TRUE) + gwodegree(0.5, fixed = TRUE) +
  edgecov(cdt_10a) + edgecov(fta_10a) +
  edgecov(cov_dist) + edgecov(pib_10a) + edgecov(cov_lang) +
  memory(type = "stability", lag = 1) +
  delrecip(mutuality = FALSE, lag = 1)

fit_10a <- btergm(f_10a, R = 500, parallel = "snow", ncpus = 4,
                  usefastglm = TRUE, verbose = TRUE)

summary(fit_10a)

idx <- which(anios_panel >= 2020)
redes_5a <- lista_redes[idx]
cdt_5a   <- cov_CDT[idx]
fta_5a   <- cov_fta[idx]
pib_5a   <- cov_PIB[idx]
# cov_dist y cov_lang son matrices únicas: btergm las recicla, no se recortan

f_5a <- redes_5a ~ edges + mutual +
  gwesp(0.69, fixed = TRUE) + ctriple +
  gwidegree(0.5, fixed = TRUE) + gwodegree(0.5, fixed = TRUE) +
  edgecov(cdt_5a) + edgecov(fta_5a) +
  edgecov(cov_dist) + edgecov(pib_5a) + edgecov(cov_lang) +
  memory(type = "stability", lag = 1) +
  delrecip(mutuality = FALSE, lag = 1)

fit_5a <- btergm(f_5a, R = 500, parallel = "snow", ncpus = 4,
                 usefastglm = TRUE, verbose = TRUE)

summary(fit_5a)

screenreg(list(fit_final, fit_10a, fit_5a),
          custom.model.names = c("2004-2024", "2015-2024", "2020-2024"))


# ── 7. Robustez por umbral de binarización ──
redes_con_umbral <- function(umbral) {
  lapply(seq_along(anios_panel), function(i) {
    bt <- base_panel %>% filter(year == anios_panel[i]) %>%
      mutate(lazo = ifelse(peso_comercio_total > umbral, 1, 0))
    as.network(hacer_matriz(bt, "lazo"), directed = TRUE,
               matrix.type = "adjacency", loops = FALSE)
  })
}

redes_50M  <- redes_con_umbral(5e7)
redes_500M <- redes_con_umbral(5e8)

sapply(redes_50M,  network.density)   # esperado 0.53 – 0.69
sapply(redes_500M, network.density)   # esperado 0.26 – 0.42

correr_umbral <- function(redes, R = 500) {
  f <- redes ~ edges + mutual +
    gwesp(0.69, fixed = TRUE) + ctriple +
    gwidegree(0.5, fixed = TRUE) + gwodegree(0.5, fixed = TRUE) +
    edgecov(cov_CDT) + edgecov(cov_fta) +
    edgecov(cov_dist) + edgecov(cov_PIB) + edgecov(cov_lang) +
    memory(type = "stability", lag = 1) +
    delrecip(mutuality = FALSE, lag = 1)
  btergm(f, R = R, parallel = "snow", ncpus = 4,
         usefastglm = TRUE, verbose = FALSE)
}

fit_50M  <- correr_umbral(redes_50M)
fit_500M <- correr_umbral(redes_500M)

screenreg(list(fit_50M, fit_final, fit_500M),
          custom.model.names = c("50M", "100M (base)", "500M"))

# ¿Cuántas díadas tienen CDT=1 y quiénes son?
i <- 21  # año 2024
which(cov_CDT[[i]] == 1, arr.ind = TRUE) |> 
  apply(1, \(x) paste(reporteros[x[1]], "->", reporteros[x[2]])) |> head(30)
sum(cov_CDT[[i]])
