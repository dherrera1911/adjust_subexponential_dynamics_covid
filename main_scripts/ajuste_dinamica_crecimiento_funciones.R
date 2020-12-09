####################################
#  Aplicacion del metodo descrito en Chowell (2017) Fitting dynamic
#  models to epidemic outbreaks with quantified uncertainty:
#  A primer for parameter uncertainty, identifiability and forecasts
#  
#  El metodo primero ajusta una funcion a los datos cumulativos de
#  la epidemia mediante regresion no lineal por minimos cuadrados.
#  Luego hace un bootstrap con data generada a partir de este
#  modelo para cuantificar la incertidumbre en los parametros.
#  El bootstrap es generado sampleando los casos nuevos
#  de cada día a partir de una Poisson con media en los
#  casos nuevos de ese día predicho por el modelo ajustado.
#
#
# Codigo original escrito por Daniel Herrera 2 de Abril 2020
# dherrera1911@gmail.com
#
###################################

####################
# devolver media y ci de las predicciones del bootstrap
###################

rango_predicciones <- function(predicciones_boot, rango_ci = 0.9) {
  predicciones_mean <- colMeans(predicciones_boot)
  predicciones_ordenadas <- apply(predicciones_boot, 2, sort)
  total_bootstrap <- dim(predicciones_ordenadas)[1]
  ci_min <- predicciones_ordenadas[round(total_bootstrap*(1-rango_ci)),]
  ci_max <- predicciones_ordenadas[round(total_bootstrap*rango_ci),]
  return(list(media = predicciones_mean, ci_min = ci_min, ci_max = ci_max))
}

##################
# ajustar exponencial a casos cumulativos
#################

modelo_exponencial <- function(casos_cumulativos) {
  # poner casos_cumulativos con dia en dataframe
  datos_reg <- data.frame(x = c(1:length(casos_cumulativos)),
                          y = casos_cumulativos)
  # estimar parametros iniciales
  fiteo_lineal <- lm(log(y+0.5) ~ 1 + x, data = datos_reg)$coefficients
  pendiente <- fiteo_lineal[2]
  # ajustar nls
  reg <- nls(y ~ a*exp(b*x), data=datos_reg, start = list(a=casos_cumulativos[1],
                                                          b=pendiente))
  return(reg)
}

####################
# ajustar polinomica a casos cumulativos
###################

modelo_poly <- function(casos_cumulativos, grado_poly=2) {
  # poner casos_cumulativos con dia en dataframe
  datos_reg <- data.frame(x = c(1:length(casos_cumulativos)),
                          y = casos_cumulativos)
  # ajustar modelo lineal
  reg <- lm(y ~ poly(x, grado_poly, raw=TRUE), data=datos_reg)
  return(reg)
}

####################
# Ajustar un modelo a los datos, y calcular la incertidumbre
# del ajuste y de las predicciones a futuro mediante bootstrap
###################

ajustar_modelo_bootstrap <- function(casos_cumulativos,
                                     tipo_modelo = "polinomico",
                                     grado_poly = 2,
                                     metodo_bootstrap = "parametric",
                                     bootstrap_n = 1000,
                                     dias_pred = 1,
                                     rango_ci = 0.9) {
  # dias a predecir
  dias_siguientes <- c(1:dias_pred) + length(casos_cumulativos)

  # ajustar modelo
  if (tipo_modelo == "polinomico") {
    modelo <- modelo_poly(casos_cumulativos, grado_poly = grado_poly)
  } else if (tipo_modelo == "exponencial") {
    modelo <- modelo_exponencial(casos_cumulativos)
  }
  # predecir proximos dias con modelo de mejor ajuste
  predicciones_mejor <- predict(modelo, newdata = data.frame(x=dias_siguientes)) 

  # generar secuencia de casos nuevos para muestrar con Poisson
  if (metodo_bootstrap %in% c("parametric", "parametrico")) {
    casos_modelo <- predict(modelo)
    casos_diarios <- c(casos_modelo[1], diff(casos_modelo))
  } else if (metodo_bootstrap %in% c("non parametric", "no parametrico")) {
    casos_diarios <- c(casos_cumulativos[1], diff(casos_cumulativos))
  }
  # a veces el polinomio da casos negativos el primer dia. Omitir
  casos_validos <- casos_diarios >= 0 # omitir casos nuevos negativos
  casos_diarios <- casos_diarios * as.integer(casos_validos) # poner 0 para que ande el sampleo

  # bootstrap sobre secuencias de casos confirmados
  n_params <- length(coefficients(modelo))
  parametros_boot <- matrix(nrow=0, ncol=n_params)
  ajuste_datos_boot <- matrix(nrow=0, ncol=length(casos_cumulativos)) 
  predicciones_boot <- matrix(nrow=0, ncol=dias_pred)
  n <- 1
  while (n <= bootstrap_n) {
    # samplear nuevos datos a partir del modelo principal
    sample_casos <- cumsum(rpois(length(casos_diarios), casos_diarios))
    sample_casos[!casos_validos] <- NA #cambiar dias no validos por NA
    # ajustar modelo a datos sampleados. Poner try por si falla nls
    try(
    expr = {
      if (tipo_modelo == "polinomico") {
        result <- modelo_boot <- modelo_poly(sample_casos, grado_poly = grado_poly)
      } else if (tipo_modelo == "exponencial") {
        modelo_boot <- modelo_exponencial(sample_casos)
      }
      # guardar parametros del modelo
      parametros_boot <- rbind(parametros_boot, coefficients(modelo_boot))
      # obtener y guardar el ajuste del modelo a los dias usados
      ajuste_datos_boot <- rbind(ajuste_datos_boot, predict(modelo_boot))
      # predecir dias siguientes y guardar
      df_dias <- data.frame(x = dias_siguientes, y = NA)
      predicciones_boot <- rbind(predicciones_boot, predict(modelo_boot, newdata = df_dias))
      n <- n + 1
    }, silent=TRUE)
  }

  # rango de ajuste a los datos
  ajuste_datos_rango <- rango_predicciones(ajuste_datos_boot, rango_ci = rango_ci)
  # rango de predicciones futuras
  predicciones_stats <- rango_predicciones(predicciones_boot, rango_ci = rango_ci)

  return(list(modelo = modelo, predicciones_modelo = predicciones_mejor,
              ajuste_datos_rango = ajuste_datos_rango,
              predicciones_boot = predicciones_stats,
              parametros_boot = parametros_boot))
}
    

