#########################################################
# Script para estudiar la tendencia de crecimiento de los casos
# cumulativos (basado en el primer script de Matias)
# 
# Se ajusta una exponencial o una polinomica a los casos cumulativos.
# Se hace un bootstrap muestreando una poisson para cada dia con
# media igual a los casos nuevos de ese dia. El bootstrap parametrico
# usa los casos nuevos predichos por el modelo, y el no parametrico
# usa los casos nuevos empiricos. Este metodo se basa en el
# descrito por G. Chowell (2017) Infectious Disease Modeling. 
# Las incertidumbres en prediccion y ajuste se basan en los
# rangos de valores ajustados/predichos por los modelos ajustados
# a estas muestras generadas.
# 
# Tambien se hacen los ajustes en una ventana movil para predecir
# el dia siguiente a la ventana, y se comparan las predicciones con
# los datos observados.
# hace lo mismo pero para dias pasados,
# 
# Descarga los datos de John Hopkins y permite elegir una lista
# de paises para hacer el ajuste.
#
# Codigo original escrito por Daniel Herrera 2 de Abril 2020
# dherrera1911@gmail.com
# 
########################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
source("ajuste_dinamica_crecimiento_funciones.R")

###################
# parametros modelo
####################

lista_paises <- c("Uruguay", "Argentina", "Slovakia", "Czechia", "Croatia")
dias_prediccion <- 12 # cantidad de días usados par predecir mañana. NaN = usar todos
ventana_temporal <- 6 # ventana temporal para testear ajuste con datos ya conocidos
grado_poly <- 2 # grado del polinomio a ajustar
n_predicciones <- 3 # cantidad de dias a predecir a partir de hoy
metodo_bootstrap = "parametrico" #parametrico samplea dias del ajuste.
#"no parametrico" samplea de los datos diarios observados
n_bootstrap = 500
rango_ci <- 0.9
tamano_texto <- 3

#######################
# descargar los datos
######################

baseURL <- "https://raw.githubusercontent.com/datasets/covid-19/master/data/time-series-19-covid-combined.csv"

covid_datos <- read.csv(baseURL) %>%
  dplyr::filter(., Country.Region %in% lista_paises) %>%
  droplevels(.) %>%
  dplyr::mutate(., Date = as.Date(as.character(Date))) %>%
  dplyr::filter(., Confirmed != 0) %>% #empezar conteo desde primer caso
  as_tibble(.)

# arreglar datos de algunos erroneos de Uruguay
cumulativo_uruguay <- c(4, 6, 8, 29, 50, 79, 94, 110, 135, 158, 162,
                        189, 217, 238, 274, 304, 310, 320, 338, 350,
                        369, 386, 400)
fecha_primer_dato <- as.Date("2020-03-13")
uruguay_temp <- data.frame(Date = c(1:length(cumulativo_uruguay)-1) + fecha_primer_dato,
                           Country.Region = "Uruguay",
                           Confirmed = cumulativo_uruguay)
covid_datos <- dplyr::filter(covid_datos, !(Country.Region == "Uruguay" &
                                            Date %in% uruguay_temp$Date)) %>%
  dplyr::bind_rows(., uruguay_temp)

# arreglar datos de Argentina.
cumulativo_arg <- c(1, 1, 2, 8, 9, 12, 17, 19, 21, 31, 34, 45,
                   56, 65, 79, 97, 128, 158, 225, 266, 301, 387,
                   502, 589, 690, 745, 820, 966, 1054, 1133, 1265,
                   1353, 1451)
fecha_primer_dato <- as.Date("2020-03-03")
arg_temp <- data.frame(Date = c(1:length(cumulativo_arg)-1) + fecha_primer_dato,
                           Country.Region = "Argentina",
                           Confirmed = cumulativo_arg)
covid_datos <- dplyr::filter(covid_datos, !(Country.Region == "Argentina" &
                                            Date %in% arg_temp$Date)) %>%
  dplyr::bind_rows(., arg_temp)

#httr::GET("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", httr::authenticate(":", ":", type="ntlm"), httr::write_disk(tf <- tempfile(fileext = ".csv")))
#allDat <- read.csv(tf) %>%
#  as_tibble(.) %>%
#  dplyr::mutate(., dateRep = as.Date(dateRep, tryFormats = c("%d/%m/%Y"))) %>%
#  dplyr::arrange(., dateRep, by_group = countriesAndTerritories)
#

#######################
# funcion para ajustar funcion a datos cumulativos
######################

ajuste_tendencia_cumulativos <- function(casos_cumulativo, dias_prediccion = NaN,
                                         ventana_temporal = 6,
                                         tipo_modelo = "polinomico",
                                         grado_poly = 2,
                                         n_predicciones = 3,
                                         metodo_bootstrap = "parametric",
                                         bootstrap_n = 1000,
                                         rango_ci = 0.9){
  # agregar conteo de dias
  casos_cumulativo$dia <- c(1:nrow(casos_cumulativo))
  dia_ultimo <- max(casos_cumulativo$dia)
  if (is.nan(dias_prediccion)) {
    dia_cutoff <- 1
    dias_prediccion <- nrow(casos_cumulativo)
  } else {
    dia_cutoff <- dia_ultimo - dias_prediccion + 1
  }

  # regresion del cumulativo de casos 
  casos_a_usar <- tail(casos_cumulativo$Confirmed, n = dias_prediccion)
  modelo <- ajustar_modelo_bootstrap(casos_a_usar,
                                     tipo_modelo = tipo_modelo,
                                     grado_poly = grado_poly,
                                     metodo_bootstrap = metodo_bootstrap,
                                     bootstrap_n = bootstrap_n,
                                     dias_pred = n_predicciones,
                                     rango_ci = rango_ci)

  # poner predicciones en data frame
  casos_futuros <- tibble(dia = c(1:n_predicciones) + dia_ultimo,
                          Confirmed = modelo$predicciones_modelo,
                          ci_min = modelo$predicciones_boot$ci_min,
                          ci_max = modelo$predicciones_boot$ci_max)
  casos_futuros$Date <- max(casos_cumulativo$Date) + c(1:n_predicciones)

  # extraer datos del ajuste
  aic_modelo <- round(AIC(modelo$modelo), digits=5)
  casos_manana <- round(modelo$predicciones_boot$media[1])
  min_casos_manana <- (max(modelo$predicciones_boot$ci_min[1],
                           casos_cumulativo$Confirmed[dia_ultimo]))
  max_casos_manana <- modelo$predicciones_boot$ci_max[1]
  nuevos_manana <- round(casos_futuros$Confirmed[1] - casos_cumulativo$Confirmed[dia_ultimo])
  min_nuevos_manana <- round(min_casos_manana -
                           tail(casos_cumulativo$Confirmed, n=1))
  max_nuevos_manana <- round(max_casos_manana -
                           tail(casos_cumulativo$Confirmed, n=1))

  # hacer labels del eje de x con los dias
  dias_labels <- seq(1, dia_ultimo+n_predicciones,
                     round((dia_ultimo + n_predicciones)/5))
  nombres_dias <- casos_cumulativo$Date[1] + dias_labels
  nombres_dias <- str_remove(as.character(nombres_dias), "2020-")
  # graficar datos y ajuste
  plot_fit <- ggplot(data = casos_cumulativo, aes(x = dia, y = Confirmed)) +
    geom_point() +
    geom_point(data = casos_futuros, aes(x = dia, y = Confirmed), color = "red") +
    geom_errorbar(data = casos_futuros, aes(ymin = ci_min,
                                            ymax = ci_max),
                  color = "red") +
    geom_text(aes(-Inf, +Inf, hjust=0, vjust=4,
                  label = paste("AIC = ", as.character(aic_modelo))),
              size = tamano_texto) +
    geom_text(aes(-Inf, +Inf, hjust=0, vjust=6,
                  label = sprintf("Casos manana = %d", casos_manana)),
              size = tamano_texto) +
    geom_text(aes(-Inf, +Inf, hjust=0, vjust=8,
                  label = sprintf("Nuevos manana = %d (CI %d: %d - %d)",
                                  nuevos_manana, round(rango_ci*100),
                                  min_nuevos_manana,
                                  max_nuevos_manana)),
              size = tamano_texto) +
    ylab("Casos confirmados") +
    xlab("Dia") +
    scale_x_discrete(limit = dias_labels, labels = nombres_dias)

  # agregar linea de tendencia estimada e intervalos de confianza
  curva_ajuste <- data.frame(x = seq(1, dias_prediccion, 0.1))
  curva_ajuste$y <- predict(modelo$modelo, newdata = curva_ajuste)
  curva_ajuste$x <- curva_ajuste$x + dia_cutoff - 1
  curva_rangos <- data.frame(dia = c(dia_cutoff:dia_ultimo),
                             ymin = modelo$ajuste_datos_rango$ci_min,
                             ymax = modelo$ajuste_datos_rango$ci_max)
  plot_fit <- plot_fit +
    geom_line(data = curva_ajuste, aes(x = x, y = y), color="blue") +
    geom_ribbon(data = curva_rangos, aes(x = dia, ymin=ymin,
                                         ymax=ymax), color="blue",
                inherit.aes = FALSE, alpha = 0.2)

  # testear el modelo con los datos ya conocidos, usando ventana movil
  dias_testear <- c((ventana_temporal+1):nrow(casos_cumulativo))
  predicciones_ventana <- data.frame()
  for (dt in dias_testear) {
    temp_df <- data.frame(dia = dt, Confirmed = NA, ymin = NA, ymax = NA)
    # seleccionar datos regresion
    dia_i <- dt - ventana_temporal
    dia_f <- dt - 1
    datos_ventana <- dplyr::filter(casos_cumulativo, dia %in% c(dia_i:dia_f))[["Confirmed"]]
    # ajustar modelo
    modelo_ventana <- ajustar_modelo_bootstrap(datos_ventana,
                                       tipo_modelo = tipo_modelo,
                                       grado_poly = grado_poly,
                                       metodo_bootstrap = metodo_bootstrap,
                                       bootstrap_n = bootstrap_n,
                                       dias_pred = 1,
                                       rango_ci = rango_ci)
    # guardar prediccion 
    temp_df$Confirmed <- modelo_ventana$predicciones_modelo
    temp_df$ymin <- modelo_ventana$predicciones_boot$ci_min
    temp_df$ymax <- modelo_ventana$predicciones_boot$ci_max
    #sacar predicciones del modelo y pegar en la tabla
    predicciones_ventana <- rbind(predicciones_ventana, temp_df)
  } 

  # graficar resultados de ventana temporal junto con casos
  plot_ventana <- ggplot(data = casos_cumulativo, aes(x = dia, y = Confirmed)) +
    geom_point() +
    geom_point(data = predicciones_ventana, aes(x = dia, y = Confirmed),
               color = "red") +
    geom_errorbar(data = predicciones_ventana, aes(ymin = ymin, ymax = ymax),
                  color = "red", width = 0) +
    xlab("Dia") +
    ylab("Casos confirmados") +
    scale_x_discrete(limit = dias_labels, labels = nombres_dias) +
    xlim(0, dia_ultimo+1)

  # poner las graficas juntas
  all_plots <- ggarrange(plot_fit, plot_ventana, ncol = 2)

  return(all_plots)
}

# hacer ajuste polinomico y exponencial para cada pais
ajustes_paises_poly <- list()
ajustes_paises_exp <- list()
for (pais in lista_paises) {
  datos_pais <- dplyr::filter(covid_datos, Country.Region == pais)
  # ajuste exponencial
  ajustes_paises_exp[[pais]] <- ajuste_tendencia_cumulativos(casos_cumulativo = datos_pais,
                                       dias_prediccion = dias_prediccion,
                                       ventana_temporal = ventana_temporal,
                                       tipo_modelo = "exponencial", 
                                       grado_poly = NaN,
                                       n_predicciones = n_predicciones,
                                       metodo_bootstrap = "parametric",
                                       bootstrap_n = n_bootstrap,
                                       rango_ci = 0.9)
  # ajuste polinomico
  ajustes_paises_poly[[pais]] <- ajuste_tendencia_cumulativos(casos_cumulativo = datos_pais,
                                       dias_prediccion = dias_prediccion,
                                       ventana_temporal = ventana_temporal,
                                       tipo_modelo = "polinomico", 
                                       grado_poly = 2,
                                       n_predicciones = n_predicciones,
                                       metodo_bootstrap = "parametric",
                                       bootstrap_n = n_bootstrap,
                                       rango_ci = 0.9)
}

s <- 12
ggsave("Uruguay_polinomio.png", ajustes_paises_poly$Uruguay, width=1.5*s, height=s, units="cm")
ggsave("Uruguay_exponencial.png", ajustes_paises_exp$Uruguay, width=1.5*s, height=s, units="cm")
ggsave("Argentina_polinomio.png", ajustes_paises_poly$Argentina, width=1.5*s, height=s, units="cm")
ggsave("Argentina_exponencial.png", ajustes_paises_exp$Argentina, width=1.5*s, height=s, units="cm")

