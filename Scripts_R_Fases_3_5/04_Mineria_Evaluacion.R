# ==============================================================================
# 04_Mineria_Evaluacion.R
# Objetivo: Estimación de TERGM por bloques anidados (M1, M2, M3), evaluación 
#           predictiva fuera de muestra (2024), pruebas de robustez (ventanas, 
#           umbrales y especificación de grado) y modelo longitudinal inverso 
#           de formación de CDT (Logit penalizado de Firth).
# ==============================================================================

#install.packages("ergm")
#install.packages("btergm")
#install.packages("parallel")
#install.packages("logistf")
#install.packages("pROC")
#install.packages("PRROC")


library(tidyverse)
library(network)
library(ergm)
library(btergm)
library(texreg)
library(parallel)


# Crear directorios de salida si no existen
dir.create("Salidas/tablas", recursive = TRUE, showWarnings = FALSE)
dir.create("Salidas/figuras", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. IMPORTACIÓN Y PREPARACIÓN DE MATRICES RELACIONALES (59 PAÍSES)
# ==============================================================================
baseERGM <- read.csv("bases_datos/base_comercio_ERGM.csv")

nucleos <- max(1, detectCores() - 1)
UMBRAL_TERGM <- 1e8 # 100 millones de USD sobre Exportaciones (X)
anios_panel <- 2004:2024

reporteros <- baseERGM %>%
  group_by(iso3_o) %>%
  summarise(n = n_distinct(year)) %>%
  filter(n == 21) %>%
  pull(iso3_o) %>%
  sort()

n_rep <- length(reporteros)

# La base ya viene balanceada e imputada desde Python (71.862 filas exactas)
base_panel <- baseERGM %>%
  arrange(year, iso3_o, iso3_d) %>%
  mutate(
    X = replace_na(X, 0),
    M = replace_na(M, 0),
    peso_comercio_total = replace_na(peso_comercio_total, 0),
    CDT_vigente = replace_na(CDT_vigente, 0),
    fta_wto = replace_na(fta_wto, 0),
    comlang_off = replace_na(comlang_off, 0),
    lazo_binario = ifelse(X > UMBRAL_TERGM, 1, 0)
  )

# Función para asegurar alineación exacta de filas y columnas por indexación
hacer_matriz <- function(df, var, nodos = reporteros) {
  mat <- matrix(0, length(nodos), length(nodos), dimnames = list(nodos, nodos))
  mat[cbind(match(df$iso3_o, nodos), match(df$iso3_d, nodos))] <- df[[var]]
  diag(mat) <- 0
  stopifnot(!anyNA(mat))
  return(mat)
}

lista_redes <- list()
cov_CDT <- list()
cov_fta <- list()
cov_PIB <- list()
pa_in <- list()
pa_out <- list()

for (i in seq_along(anios_panel)) {
  df_t <- base_panel %>% filter(year == anios_panel[i])
  A <- hacer_matriz(df_t, "lazo_binario")
  lista_redes[[i]] <- as.network(A, directed = TRUE)
  
  cov_CDT[[i]] <- hacer_matriz(df_t, "CDT_vigente")
  cov_fta[[i]] <- hacer_matriz(df_t, "fta_wto")
  cov_PIB[[i]] <- hacer_matriz(df_t, "ln_gdp_prod")
  
  idx_rezago <- max(i - 1, 1)
  A_rezago <- hacer_matriz(base_panel %>% filter(year == anios_panel[idx_rezago]), "lazo_binario")
  pa_in[[i]]  <- matrix(sqrt(colSums(A_rezago)), n_rep, n_rep, byrow = TRUE, dimnames = list(reporteros, reporteros))
  pa_out[[i]] <- matrix(sqrt(rowSums(A_rezago)), n_rep, n_rep, byrow = FALSE, dimnames = list(reporteros, reporteros))
}

b24 <- base_panel %>% filter(year == 2024)
mat_dist <- hacer_matriz(b24, "ln_dist")
mat_lang <- hacer_matriz(b24, "comlang_off")

# ==============================================================================
# 2. PRUEBAS DE VERIFICACIÓN DE COHERENCIA METODOLÓGICA
# ==============================================================================
cat("\n--- Ejecutando Verificaciones de Coherencia ---\n")
stopifnot(length(reporteros) == 59)
stopifnot(nrow(base_panel) == 59 * 58 * 21) # 71.862 díadas dirigidas

cdt_check <- function(o, d, y) {
  base_panel$CDT_vigente[base_panel$iso3_o == o & base_panel$iso3_d == d & base_panel$year == y]
}
stopifnot(cdt_check("COL", "FRA", 2022) == 0, cdt_check("COL", "FRA", 2023) == 1) # Vigencia, no firma
stopifnot(cdt_check("BRA", "DEU", 2005) == 1, cdt_check("BRA", "DEU", 2006) == 0) # Terminación aplicada
stopifnot(cdt_check("CAN", "CHN", 2016) == cdt_check("CAN", "CHN", 2017))         # Sin salto por MLI

# Simetría del tratado bilateral y alineación de matrices
stopifnot(all(cov_CDT[[21]] == t(cov_CDT[[21]])))
stopifnot(all(rownames(cov_PIB[[1]]) == colnames(cov_PIB[[1]])))

cat("Evolución de convenios activos entre los 59 países (2004-2024):\n")
print(setNames(sapply(cov_CDT, function(m) sum(m[upper.tri(m)])), anios_panel))

cat("Densidad anual de la red dirigida de exportaciones (X > 100M USD):\n")
print(setNames(round(sapply(lista_redes, network.density), 3), anios_panel))

A_21 <- as.matrix(lista_redes[[21]])
cat("Proporción de díadas asimétricas en 2024 (verifica red dirigida):",
    round(mean(A_21[upper.tri(A_21)] != t(A_21)[upper.tri(A_21)]), 3), "\n")

# ==============================================================================
# 3. OBJETIVO 3: ESTIMACIÓN TERGM POR BLOQUES ANIDADOS (M1 -> M2 -> M3)
# ==============================================================================
cat("\n--- Estimando Bloques TERGM (2004-2023) ---\n")
redes_train  <- lista_redes[1:20]
pa_in_train  <- pa_in[1:20]
pa_out_train <- pa_out[1:20]
cdt_train    <- cov_CDT[1:20]
fta_train    <- cov_fta[1:20]
pib_train    <- cov_PIB[1:20]
reps <- 500

set.seed(2026)
f_M1 <- redes_train ~ edges + mutual + gwesp(0.69, fixed = TRUE) + ctriple +
  memory(type = "stability", lag = 1) + delrecip(mutuality = FALSE, lag = 1) +
  edgecov(pa_in_train) + edgecov(pa_out_train)
fit_M1 <- btergm(f_M1, R = reps, parallel = "snow", ncpus = nucleos)

set.seed(2026)
f_M2 <- update(f_M1, ~ . + edgecov(fta_train) + edgecov(mat_dist) + edgecov(pib_train) + edgecov(mat_lang))
fit_M2 <- btergm(f_M2, R = reps, parallel = "snow", ncpus = nucleos)

set.seed(2026)
f_M3 <- update(f_M2, ~ . + edgecov(cdt_train))
fit_M3 <- btergm(f_M3, R = reps, parallel = "snow", ncpus = nucleos)

# Nombres unificados con la Ecuación 3 y la Tabla 1 de la tesis
nombres_coef <- c(
  "edges"                        = "Edges (Densidad base)",
  "mutual"                       = "Mutual (Reciprocidad)",
  "gwesp.fixed.0.69"             = "GWESP (0,69)",
  "ctriple"                      = "CTriple (Ciclicidad)",
  "edgecov.memory[[i]]"          = "Memory (Estabilidad t-1)",
  "edgecov.delrecip[[i]]"        = "DelRecip (Reciprocidad t-1)",
  "edgecov.pa_in_train[[i]]"     = "Popularidad (PA_in)",
  "edgecov.pa_out_train[[i]]"    = "Actividad (PA_out)",
  "edgecov.fta_train[[i]]"       = "TLC / OMC",
  "edgecov.mat_dist"             = "ln(Distancia)",
  "edgecov.pib_train[[i]]"       = "ln(PIB_i * PIB_j)",
  "edgecov.mat_lang"             = "Idioma común",
  "edgecov.cdt_train[[i]]"       = "CDT vigente"
)

texreg(
  list(fit_M1, fit_M2, fit_M3),
  custom.model.names = c("M1: Endógeno", "M2: +Gravedad", "M3: +CDT"),
  custom.coef.map = as.list(nombres_coef),
  file = "Salidas/tablas/resultados_TERGM_bloques.tex"
)

# ==============================================================================
# 4. EVALUACIÓN FUERA DE MUESTRA EN 2024 (M1, M2 y M3 + Chequeo sin ARE/VNM)
# ==============================================================================
library(pROC)
library(PRROC)

# A. Construir matrices rezagadas manuales de 2023 (t=20) para predecir 2024 (t=21)
A23      <- as.matrix(lista_redes[[20]])
mem_24   <- ifelse(A23 == 1, 1, -1) # Memory (stability): 1 si había lazo, -1 si no había
delr_24  <- t(A23)                  # DelRecip: transpuesta de A23
diag(mem_24) <- 0; diag(delr_24) <- 0

net24   <- lista_redes[[21]]
pain_24 <- pa_in[[21]]
pout_24 <- pa_out[[21]]
fta_24  <- cov_fta[[21]]
pib_24  <- cov_PIB[[21]]
cdt_24  <- cov_CDT[[21]]

# Fórmulas equivalentes en ergm puro para 2024
f_ergm_M1 <- net24 ~ edges + mutual + gwesp(0.69, fixed = TRUE) + ctriple +
  edgecov(mem_24) + edgecov(delr_24) + edgecov(pain_24) + edgecov(pout_24)

f_ergm_M2 <- update(f_ergm_M1, ~ . + edgecov(fta_24) + edgecov(mat_dist) +
                      edgecov(pib_24) + edgecov(mat_lang))

f_ergm_M3 <- update(f_ergm_M2, ~ . + edgecov(cdt_24))

# B. Prueba diagnóstica de simulate() directo (para ver si ergm 4.12 simula sin error)
sim_test <- simulate(f_ergm_M1, nsim = 1, coef = unname(coef(fit_M1)))
print(sim_test)

# C. Cálculo exacto de AUC-ROC y AUC-PR fuera de muestra (2024) vía ergmMPLE
eval_mple_2024 <- function(formula_ergm, modelo_btergm, nombre, excluir_imputados = FALSE) {
  mple <- ergmMPLE(formula_ergm, output = "matrix")
  probs <- as.numeric(plogis(mple$predictor %*% unname(coef(modelo_btergm))))
  y_obs <- mple$response
  
  roc_obj <- roc(y_obs, probs, direction = "<", quiet = TRUE)
  pr_obj  <- pr.curve(scores.class0 = probs[y_obs == 1],
                      scores.class1 = probs[y_obs == 0])
  
  tibble(
    Modelo  = nombre,
    AUC_ROC = round(as.numeric(auc(roc_obj)), 4),
    AUC_PR  = round(pr_obj$auc.integral, 4)
  )
}

tabla_auc_tergm <- bind_rows(
  eval_mple_2024(f_ergm_M1, fit_M1, "M1: Endógeno"),
  eval_mple_2024(f_ergm_M2, fit_M2, "M2: +Gravedad"),
  eval_mple_2024(f_ergm_M3, fit_M3, "M3: +CDT")
)
cat("\n--- Desempeño Predictivo Fuera de Muestra en 2024 (vía MPLE) ---\n")
print(tabla_auc_tergm)
write.csv(tabla_auc_tergm, "Salidas/tablas/desempeno_bloques_2024.csv", row.names = FALSE)

# D. Chequeo de sensibilidad excluyendo ARE y VNM en 2024 (sobre las 100 redes simuladas con ergm)
set.seed(2026)
sim_M3_100 <- simulate(f_ergm_M3, nsim = 100, coef = unname(coef(fit_M3)))
prob_mat_M3 <- Reduce("+", lapply(sim_M3_100, as.matrix)) / 100

A_2024_full <- as.matrix(net24)
diag(A_2024_full) <- NA; diag(prob_mat_M3) <- NA
auc_sim_full <- auc(roc(as.vector(A_2024_full), as.vector(prob_mat_M3), direction = "<", quiet = TRUE))

idx_limpio <- which(!(reporteros %in% c("ARE", "VNM")))
auc_sim_limpio <- auc(roc(as.vector(A_2024_full[idx_limpio, idx_limpio]),
                          as.vector(prob_mat_M3[idx_limpio, idx_limpio]), direction = "<", quiet = TRUE))

cat("\nAUC-ROC 2024 por simulación (100 redes) - Completo:", round(as.numeric(auc_sim_full), 4),
    "| Sin ARE y VNM:", round(as.numeric(auc_sim_limpio), 4), "\n")

# E. Gráfico de Bondad de Ajuste (GOF estructural de M3)
gof_M3_in <- btergm::gof(fit_M3, nsim = 100, statistics = c(esp, dsp,geodesic, ideg, odeg, triad.directed, rocpr))
png("Salidas/figuras/gof_M3_prediccion2024.png", width = 2400, height = 1800, res = 300)
plot(gof_M3_in)
dev.off()
# ==============================================================================
# 5. ROBUSTEZ I: VENTANAS TEMPORALES (2015-2023 y 2019-2023)
# ==============================================================================
cat("\n--- Ejecutando Robustez I: Ventanas Temporales ---\n")

correr_ventana <- function(idx) {
  r_v    <- lista_redes[idx]
  pain_v <- pa_in[idx]
  pout_v <- pa_out[idx]
  fta_v  <- cov_fta[idx]
  pib_v  <- cov_PIB[idx]
  cdt_v  <- cov_CDT[idx]
  
  f_v <- r_v ~ edges + mutual + gwesp(0.69, fixed = TRUE) + ctriple +
    memory(type = "stability", lag = 1) + delrecip(mutuality = FALSE, lag = 1) +
    edgecov(pain_v) + edgecov(pout_v) + edgecov(fta_v) +
    edgecov(mat_dist) + edgecov(pib_v) + edgecov(mat_lang) + edgecov(cdt_v)
  
  set.seed(2026)
  btergm(f_v, R = reps, parallel = "snow", ncpus = nucleos)
}

fit_10a <- correr_ventana(12:20) # 2015-2023
fit_5a  <- correr_ventana(16:20) # 2019-2023

texreg(
  list(fit_M3, fit_10a, fit_5a),
  custom.model.names = c("2004-2023 (Base)", "2015-2023", "2019-2023"),
  file = "Salidas/tablas/robustez_ventanas.tex"
)

# ==============================================================================
# 6. ROBUSTEZ II: UMBRALES DE BINARIZACIÓN (50M y 500M sobre Exportaciones X)
# ==============================================================================
cat("\n--- Ejecutando Robustez II: Umbrales de Binarización ---\n")

correr_umbral <- function(umbral, R_sim = reps) {
  redes_u <- list(); pa_in_u <- list(); pa_out_u <- list()
  for (i in seq_along(anios_panel)) {
    df_t <- base_panel %>%
      filter(year == anios_panel[i]) %>%
      mutate(lazo = ifelse(X > umbral, 1, 0))
    redes_u[[i]] <- as.network(hacer_matriz(df_t, "lazo"), directed = TRUE)
    
    idx_rezago <- max(i - 1, 1)
    df_rez <- base_panel %>%
      filter(year == anios_panel[idx_rezago]) %>%
      mutate(lazo = ifelse(X > umbral, 1, 0))
    A_rez <- hacer_matriz(df_rez, "lazo")
    pa_in_u[[i]]  <- matrix(sqrt(colSums(A_rez)), n_rep, n_rep, byrow = TRUE, dimnames = list(reporteros, reporteros))
    pa_out_u[[i]] <- matrix(sqrt(rowSums(A_rez)), n_rep, n_rep, byrow = FALSE, dimnames = list(reporteros, reporteros))
  }
  
  r_u    <- redes_u[1:20]
  pain_u <- pa_in_u[1:20]
  pout_u <- pa_out_u[1:20]
  
  f_u <- r_u ~ edges + mutual + gwesp(0.69, fixed = TRUE) + ctriple +
    memory(type = "stability", lag = 1) + delrecip(mutuality = FALSE, lag = 1) +
    edgecov(pain_u) + edgecov(pout_u) + edgecov(fta_train) +
    edgecov(mat_dist) + edgecov(pib_train) + edgecov(mat_lang) + edgecov(cdt_train)
  
  set.seed(2026)
  return(btergm(f_u, R = R_sim, parallel = "snow", ncpus = nucleos))
}

fit_50M  <- correr_umbral(5e7)
fit_500M <- correr_umbral(5e8)

texreg(
  list(fit_50M, fit_M3, fit_500M),
  custom.model.names = c("Umbral 50M", "Umbral 100M (Base)", "Umbral 500M"),
  file = "Salidas/tablas/robustez_umbrales.tex"
)

# ==============================================================================
# 7. ROBUSTEZ III: ESPECIFICACIÓN ALTERNATIVA DE GRADO (GWIDEGREE / GWODEGREE)
# ==============================================================================
cat("\n--- Ejecutando Robustez III: gwidegree y gwodegree ---\n")
set.seed(2026)
f_M3_gwdeg <- redes_train ~ edges + mutual + gwesp(0.69, fixed = TRUE) + ctriple +
  gwidegree(0.5, fixed = TRUE) + gwodegree(0.5, fixed = TRUE) +
  memory(type = "stability", lag = 1) + delrecip(mutuality = FALSE, lag = 1) +
  edgecov(fta_train) + edgecov(mat_dist) + edgecov(pib_train) + edgecov(mat_lang) + edgecov(cdt_train)

fit_M3_gwdeg <- btergm(f_M3_gwdeg, R = reps, parallel = "snow", ncpus = nucleos)

texreg(
  list(fit_M3, fit_M3_gwdeg),
  custom.model.names = c("M3: Apego Preferencial (PA)", "M3: GW-Degree (0,5)"),
  file = "Salidas/tablas/robustez_grados.tex"
)

# ==============================================================================
# MODELO DE FORMACIÓN DE CDT (Logit Penalizado de Firth)
# ==============================================================================
library(logistf)

cat("\n--- Preparando Panel No Dirigido en Riesgo para Formación de CDT ---\n")

# Calcular grado en la red de convenios para controlar estructura previa de tratados
grados_cdt <- map_dfr(seq_along(anios_panel), function(i) {
  tibble(year = anios_panel[i], iso = reporteros, deg_cdt = rowSums(cov_CDT[[i]]))
})

panel_cdt <- base_panel %>%
  mutate(a = pmin(iso3_o, iso3_d), b = pmax(iso3_o, iso3_d)) %>%
  group_by(year, a, b) %>%
  summarise(
    CDT         = max(CDT_vigente),
    comercio    = log1p(sum(X)),          # Flujo bilateral sumado X_ab + X_ba
    lazo_activo = max(lazo_binario),      # Al menos un enlace comercial > 100M USD
    fta_wto     = max(fta_wto),
    ln_dist     = first(ln_dist),
    comlang_off = first(comlang_off),
    ln_gdp_prod = first(ln_gdp_prod),
    .groups     = "drop"
  ) %>%
  left_join(grados_cdt, by = c("year", "a" = "iso")) %>%
  left_join(grados_cdt, by = c("year", "b" = "iso"), suffix = c("_a", "_b")) %>%
  mutate(deg_cdt_sum = deg_cdt_a + deg_cdt_b) %>%
  arrange(a, b, year) %>%
  group_by(a, b) %>%
  mutate(
    CDT_nuevo   = ifelse(lag(CDT) == 0 & CDT == 1, 1, 0),
    CDT_t1      = lag(CDT),
    comercio_t1 = lag(comercio),
    lazo_t1     = lag(lazo_activo),
    pib_t1      = lag(ln_gdp_prod),
    fta_t1      = lag(fta_wto),
    deg_cdt_t1  = lag(deg_cdt_sum)
  ) %>%
  ungroup() %>%
  filter(CDT_t1 == 0) %>% # Díadas en riesgo de formar un nuevo tratado (sin excluir 2017)
  drop_na(CDT_nuevo, comercio_t1, pib_t1)

train_cdt <- panel_cdt %>% filter(year <= 2018)
test_cdt  <- panel_cdt %>% filter(year >= 2019)

cat("\nEventos de entrada en vigor (CDT_nuevo == 1) en Train (2005-2018) y Test (2019-2024):\n")
print(table(Train = train_cdt$CDT_nuevo))
print(table(Test  = test_cdt$CDT_nuevo))

# A. Modelo 1 (Referencia Gravitacional pura)
logit_grav <- logistf(
  CDT_nuevo ~ pib_t1 + ln_dist + comlang_off + fta_t1,
  data = train_cdt, flic = TRUE
)

# B. Modelo 2 (Gravedad + Volumen Comercial Rezagado)
logit_comercio <- logistf(
  CDT_nuevo ~ comercio_t1 + pib_t1 + ln_dist + comlang_off + fta_t1,
  data = train_cdt, flic = TRUE
)

# C. Modelo 3 (Gravedad + Comercio + Actividad en Red de CDT previa)
logit_completo <- logistf(
  CDT_nuevo ~ comercio_t1 + lazo_t1 + deg_cdt_t1 + pib_t1 + ln_dist + comlang_off + fta_t1,
  data = train_cdt, flic = TRUE
)

cat("\n--- RESULTADO FIRTH: Modelo Gravedad  ---\n")
summary(logit_grav)
cat("\n--- RESULTADO FIRTH: Modelo Gravedad + Comercio ---\n")
summary(logit_comercio)
cat("\n--- RESULTADO FIRTH: Modelo Completo ---\n")
summary(logit_completo)

# Evaluación fuera de muestra (2019-2024) con direction = "<", IC 95% y AUC-PR
evaluar_firth <- function(mod, nombre) {
  preds <- predict(mod, newdata = test_cdt, type = "response")
  roc_obj <- roc(test_cdt$CDT_nuevo, preds, direction = "<", quiet = TRUE)
  ci_obj  <- ci.auc(roc_obj, quiet = TRUE)
  pr_obj  <- pr.curve(scores.class0 = preds[test_cdt$CDT_nuevo == 1],
                      scores.class1 = preds[test_cdt$CDT_nuevo == 0])
  tibble(
    Modelo     = nombre,
    AUC_ROC    = round(as.numeric(auc(roc_obj)), 4),
    IC95_Inf   = round(ci_obj[1], 4),
    IC95_Sup   = round(ci_obj[3], 4),
    AUC_PR     = round(pr_obj$auc.integral, 4)
  )
}

tabla_eval_cdt <- bind_rows(
  evaluar_firth(logit_grav,     "1. Solo Gravedad"),
  evaluar_firth(logit_comercio, "2. Gravedad + Comercio (t-1)"),
  evaluar_firth(logit_completo, "3. Gravedad + Comercio + Red CDT (t-1)")
)

roc.test(roc(test_cdt$CDT_nuevo, predict(logit_comercio, newdata = test_cdt, type = "response"), direction = "<", quiet = TRUE), 
         roc(test_cdt$CDT_nuevo, predict(logit_completo, newdata = test_cdt, type = "response"), direction = "<", quiet = TRUE))

cat("\n--- Desempeño Predictivo de Nuevos CDT (2019-2024) ---\n")
print(tabla_eval_cdt)
write.csv(tabla_eval_cdt, "Salidas/tablas/evaluacion_firth_cdt.csv", row.names = FALSE)

# Exportar candidatos para Colombia en 2024 (sin CDT vigente)
candidatos_col_2024 <- panel_cdt %>%
  filter(year == 2024 & (a == "COL" | b == "COL")) %>%
  mutate(
    Socio = ifelse(a == "COL", b, a),
    Prob_Firth = round(predict(logit_completo, newdata = ., type = "response"), 4)
  ) %>%
  arrange(desc(Prob_Firth)) %>%
  select(Socio, Prob_Firth, comercio_t1, lazo_t1, fta_t1, deg_cdt_t1)

cat("\n--- Top 15 Candidatos a CDT con Colombia según Modelo Longitudinal ---\n")
print(head(candidatos_col_2024, 15))
write.csv(head(candidatos_col_2024, 15), "Salidas/tablas/ranking_candidatos_colombia_2024.csv", row.names = FALSE)