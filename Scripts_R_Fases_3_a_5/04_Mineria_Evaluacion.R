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

#install.packages("ergm")
#install.packages("btergm")
#install.packages("parallel")
install.packages("logistf")
install.packages("pROC")

library(tidyverse)
library(network)
library(ergm)
library(btergm)
library(texreg)
library(parallel)
library(logistf) 
library(pROC)

#===========================================
# Importacion y preparación para red desde TERGM y VERGM
#===========================================

baseERGM <- read.csv("bases_datos/base_comercio_ERGM.csv")

nucleos <- detectCores() - 1
UMBRAL_TERGM <- 1e8


reporteros <- baseERGM %>% group_by(iso3_o) %>% 
  summarise(n = n_distinct(year)) %>% filter(n == 21) %>% 
  pull(iso3_o) %>% sort()


base_panel <- expand_grid(year = 2004:2024, iso3_o = reporteros, iso3_d = reporteros) %>%
  filter(iso3_o != iso3_d) %>%
  left_join(baseERGM, by = c("year", "iso3_o", "iso3_d")) %>%
  group_by(iso3_o, iso3_d) %>%
  fill(ln_dist, comlang_off, .direction = "downup") %>% 
  ungroup() %>%
  mutate(
    peso_comercio_total = replace_na(peso_comercio_total, 0),
    CDT_vigente = replace_na(CDT_vigente, 0),
    fta_wto = replace_na(fta_wto, 0),
    lazo_binario = ifelse(peso_comercio_total > UMBRAL_TERGM, 1, 0),
    ln_gdp_prod = ifelse(is.na(ln_gdp_prod), median(ln_gdp_prod, na.rm=TRUE), ln_gdp_prod)
  )

# Función para asegurar alineación exacta de filas y columnas (por indexación)
hacer_matriz <- function(df, var, nodos = reporteros) {
  mat <- matrix(0, length(nodos), length(nodos), dimnames = list(nodos, nodos))
  mat[cbind(match(df$iso3_o, nodos), match(df$iso3_d, nodos))] <- df[[var]]
  diag(mat) <- 0
  stopifnot(!anyNA(mat))
  return(mat)
}

lista_redes <- list(); cov_CDT <- list(); cov_fta <- list(); cov_PIB <- list()
pa_in <- list(); pa_out <- list()
anios_panel <- 2004:2024
n_rep <- length(reporteros)

for (i in seq_along(anios_panel)) {
  df_t <- base_panel %>% filter(year == anios_panel[i])
  A <- hacer_matriz(df_t, "lazo_binario")
  lista_redes[[i]] <- as.network(A, directed=TRUE)
  
  cov_CDT[[i]] <- hacer_matriz(df_t, "CDT_vigente")
  cov_fta[[i]] <- hacer_matriz(df_t, "fta_wto")
  cov_PIB[[i]] <- hacer_matriz(df_t, "ln_gdp_prod")
  
  idx_rezago <- max(i - 1, 1)
  A_rezago <- hacer_matriz(base_panel %>% filter(year == anios_panel[idx_rezago]), "lazo_binario")
  pa_in[[i]]  <- matrix(sqrt(colSums(A_rezago)), n_rep, n_rep, byrow = TRUE)
  pa_out[[i]] <- matrix(sqrt(rowSums(A_rezago)), n_rep, n_rep, byrow = FALSE)
}

b24 <- base_panel %>% filter(year == 2024)
mat_dist <- hacer_matriz(b24, "ln_dist"); mat_lang <- hacer_matriz(b24, "comlang_off")

# ==============================================================================
# 3. OBJETIVO 3: MINERÍA DE DATOS (Contribución por Bloques: M1 -> M2 -> M3)
# ==============================================================================
redes_train <- lista_redes[1:20]
pa_in_train <- pa_in[1:20]; pa_out_train <- pa_out[1:20]
cdt_train <- cov_CDT[1:20]; fta_train <- cov_fta[1:20]; pib_train <- cov_PIB[1:20]
reps <- 500

f_M1 <- redes_train ~ edges + mutual + gwesp(0.69, fixed = TRUE) + ctriple +
  memory(type = "stability", lag = 1) + delrecip(mutuality = FALSE, lag = 1) +
  edgecov(pa_in_train) + edgecov(pa_out_train) 
fit_M1 <- btergm(f_M1, R = reps, parallel = "snow", ncpus = nucleos)

f_M2 <- update(f_M1, ~ . + edgecov(fta_train) + edgecov(mat_dist) + edgecov(pib_train) + edgecov(mat_lang))
fit_M2 <- btergm(f_M2, R = reps, parallel = "snow", ncpus = nucleos)

f_M3 <- update(f_M2, ~ . + edgecov(cdt_train))
fit_M3 <- btergm(f_M3, R = reps, parallel = "snow", ncpus = nucleos)

texreg(list(fit_M1, fit_M2, fit_M3), custom.model.names = c("M1: Endógeno", "M2: +Gravedad", "M3: +CDT"), file = "resultados_TERGM_bloques.tex")

# ==============================================================================
# 4. EVALUACIÓN FUERA DE MUESTRA (Prediciendo 2024)
# ==============================================================================
cat("\n--- Evaluando GOF Fuera de Muestra (Red 2024) ---\n")
f_M3_test <- lista_redes[20:21] ~ edges + mutual + gwesp(0.69, fixed = TRUE) + ctriple +
  memory(type = "stability", lag = 1) + delrecip(mutuality = FALSE, lag = 1) +
  edgecov(pa_in[20:21]) + edgecov(pa_out[20:21]) + edgecov(cov_fta[20:21]) + 
  edgecov(mat_dist) + edgecov(cov_PIB[20:21]) + edgecov(mat_lang) + edgecov(cov_CDT[20:21])

gof_M3 <- gof(fit_M3, nsim = 100, target = lista_redes[[21]],
              formula = f_M3_test, coef = coef(fit_M3), 
              statistics = c(rocpr, ideg, odeg, esp, geodesic))
png("gof_M3_prediccion2024.png", width=800, height=600); plot(gof_M3); dev.off()

# ==============================================================================
# 5. ROBUSTEZ I: VENTANAS TEMPORALES
# ==============================================================================
cat("\n--- Ejecutando Robustez: Ventanas Temporales ---\n")

# Última década (2015-2023 train) -> Índices 12 a 20
idx_10a <- 12:20
f_10a <- lista_redes[idx_10a] ~ edges + mutual + gwesp(0.69, fixed = TRUE) + ctriple +
  memory(type = "stability", lag = 1) + delrecip(mutuality = FALSE, lag = 1) +
  edgecov(pa_in[idx_10a]) + edgecov(pa_out[idx_10a]) + edgecov(cov_fta[idx_10a]) + 
  edgecov(mat_dist) + edgecov(cov_PIB[idx_10a]) + edgecov(mat_lang) + edgecov(cov_CDT[idx_10a])
fit_10a <- btergm(f_10a, R = reps, parallel = "snow", ncpus = nucleos)

# Último quinquenio (2019-2023 train) -> Índices 16 a 20
idx_5a <- 16:20
f_5a <- lista_redes[idx_5a] ~ edges + mutual + gwesp(0.69, fixed = TRUE) + ctriple +
  memory(type = "stability", lag = 1) + delrecip(mutuality = FALSE, lag = 1) +
  edgecov(pa_in[idx_5a]) + edgecov(pa_out[idx_5a]) + edgecov(cov_fta[idx_5a]) + 
  edgecov(mat_dist) + edgecov(cov_PIB[idx_5a]) + edgecov(mat_lang) + edgecov(cov_CDT[idx_5a])
fit_5a <- btergm(f_5a, R = reps, parallel = "snow", ncpus = nucleos)

texreg(list(fit_M3, fit_10a, fit_5a), custom.model.names = c("2004-2023 (Base)", "2015-2023", "2019-2023"), file = "robustez_ventanas.tex")

# ==============================================================================
# 6. ROBUSTEZ II: UMBRALES DE BINARIZACIÓN
# ==============================================================================
cat("\n--- Ejecutando Robustez: Umbrales de Binarización ---\n")

correr_umbral <- function(umbral, R_sim = reps) {
  redes_u <- list(); pa_in_u <- list(); pa_out_u <- list()
  for (i in seq_along(anios_panel)) {
    df_t <- base_panel %>% filter(year == anios_panel[i]) %>% mutate(lazo = ifelse(peso_comercio_total > umbral, 1, 0))
    redes_u[[i]] <- as.network(hacer_matriz(df_t, "lazo"), directed = TRUE)
    
    idx_rezago <- max(i - 1, 1)
    df_rez <- base_panel %>% filter(year == anios_panel[idx_rezago]) %>% mutate(lazo = ifelse(peso_comercio_total > umbral, 1, 0))
    A_rez <- hacer_matriz(df_rez, "lazo")
    pa_in_u[[i]]  <- matrix(sqrt(colSums(A_rez)), n_rep, n_rep, byrow = TRUE)
    pa_out_u[[i]] <- matrix(sqrt(rowSums(A_rez)), n_rep, n_rep, byrow = FALSE)
  }
  
  f_u <- redes_u[1:20] ~ edges + mutual + gwesp(0.69, fixed = TRUE) + ctriple +
    memory(type = "stability", lag = 1) + delrecip(mutuality = FALSE, lag = 1) +
    edgecov(pa_in_u[1:20]) + edgecov(pa_out_u[1:20]) + edgecov(fta_train) +
    edgecov(mat_dist) + edgecov(pib_train) + edgecov(mat_lang) + edgecov(cdt_train)
  
  return(btergm(f_u, R = R_sim, parallel = "snow", ncpus = nucleos))
}

fit_50M  <- correr_umbral(5e7)
fit_500M <- correr_umbral(5e8)

texreg(list(fit_50M, fit_M3, fit_500M), custom.model.names = c("Umbral 50M", "Umbral 100M (Base)", "Umbral 500M"), file = "robustez_umbrales.tex")

# ==============================================================================
# 7. OBJETIVO 4: MODELO DE FORMACIÓN DE CDT (Logit de eventos raros)
# ==============================================================================
cat("\n--- Preparando Panel de Riesgo para Formación de CDT ---\n")

panel_cdt <- base_panel %>%
  mutate(a = pmin(iso3_o, iso3_d), b = pmax(iso3_o, iso3_d)) %>%
  group_by(year, a, b) %>%
  summarise(
    CDT = max(CDT_vigente), comercio = log1p(sum(peso_comercio_total)),
    lazo = max(lazo_binario), ln_dist = first(ln_dist),
    comlang_off = first(comlang_off), ln_gdp_prod = first(ln_gdp_prod), 
    .groups = "drop"
  ) %>%
  arrange(a, b, year) %>%
  group_by(a, b) %>%
  mutate(
    CDT_nuevo = ifelse(lag(CDT) == 0 & CDT == 1, 1, 0),
    CDT_t1 = lag(CDT), comercio_t1 = lag(comercio), pib_t1 = lag(ln_gdp_prod)
  ) %>%
  ungroup() %>%
  filter(year != 2017 & CDT_t1 == 0) %>% 
  drop_na(CDT_nuevo, comercio_t1)

train_cdt <- panel_cdt %>% filter(year <= 2018)
test_cdt  <- panel_cdt %>% filter(year >= 2019)

modelo_formacion_cdt <- logistf(CDT_nuevo ~ comercio_t1 + ln_dist + comlang_off + pib_t1, data = train_cdt)

cat("\n--- RESULTADO: ¿El comercio predice la creación del CDT? ---\n")
summary(modelo_formacion_cdt)

predicciones_cdt <- predict(modelo_formacion_cdt, newdata = test_cdt, type = "response")
roc_cdt <- roc(test_cdt$CDT_nuevo, predicciones_cdt)
cat("\nAUC-ROC para Predicción de Tratados Nuevos (2019-2024):", auc(roc_cdt), "\n")