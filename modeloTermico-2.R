#' ---
#' title: "Difusión del calor en noches de helada: experimentos y modelo simplificado"
#' author: "Ana Laura Diedrichs"
#' 
#' ---
#' 
#' \pagenumbering{gobble}
#' Introducción al modelo y sus objetivos 
#' ---
#' Las heladas constituyen uno de los accidentes de tiempo que causan grandes pérdidas 
#' económicas, a la agricultura en la Argentina y gran parte del mundo, e impacto
#' social al verse afectado los cultivos de los productores, debido a que no son fenómenos 
#' locales sino extensivos. Existen varias definiciones de helada como considerar helada a las 
#' temperaturas mínimas menores a 0 °C. La más apropiada desde el punto de vista agronómico
#' es considerar a la helada como el evento meteorológico que ocurre cuando los cultivos 
#'         y otras plantas experimentan daño por congelación. El daño causado por heladas ocurre 
#' cuando las temperaturas están debajo de un límite tolerable para
#' los cultivos. El umbral de resistencia de las plantas al frío varía de acuerdo 
#' al estado fenológico en el que se encuentren (floración, frutos o yemas presentes, etc).
#' 
#'  El descenso térmico que caracteriza una noche de helada está determinado por 
#'  un balance calórico negativo
#' del entorno. Intervienen varias causas de pérdida del calor del sistema que denominaremos 
#' como flujos salientes
#' o entrantes cuya unidad es $W/m^2$ (watts por metro cuadrado). El presente trabajo pretende
#' presentar un modelo simplificado de propagación de calor.
#' 
#' **OBJETIVO: Evaluar la propagación de calor $Q$ en la superficie del suelo teniendo en cuenta los flujos de propagación
#' de calor en el suelo $Q_s$, pérdida de calor por radiación $Q_r$ y por convección $Q_c$, 
#' resultando como ecuación de balance la siguiente **
#' 
#' $$Q = Q_s - Q_r + Q_c$$
#' 
#' Para ello consideraremos
#' una disminución de la temperatura proporcional a su coeficiente de cambio negativo. Definimos tres temperaturas en grados Kelvin:
#'  
#' * $T_p$ a 1 m de profundidad y su coeficiente $c_p$, 
#' * $T_s$ en superficie del suelo y su coeficiente $c_s$
#' * $T_a$ temperatura a 1 m de altura sobre superficie y su coeficiente $c_a$.
#' 
#' Esto lo realizamos para simular temperaturas en disminución. 
#' En un sistema más realista estas serían variables de entrada.
#' 
#' Los flujos entrantes y salientes serían:
#' 
#' * Propagación de calor hacia la superficie del suelo calculada como 
#' $$Q_s=-k* Area_{subsuelo} * (T_s - T_p)$$,
#'  donde k es la constante de conductividad característica propia del material.
#'  
#' * Irradiación de onda larga durante la noche emitida por el suelo calculada como 
#' $$Q_r=-\epsilon * \sigma *Area_{superficie} * T_s^4 $$ 
#' , donde $\epsilon$ es el coeficiente de emisividad del suelo, 
#' $\sigma = 5.67 * 10^{-8} \frac{J}{m^2 s K}$ la constante de Stefan-Boltzmann.
#'  
#' * Pérdida de calor por convección determinada como 
#' $$Q_c = h * Area_{superficie} *  (T_a - T_p) $$, siendo  $h$ una 
#' constante de convección
#' 
#' ![Gráfico a modo ilustrativo del modelo de difusión](modelo_difusion.png)
#' 
#' Este modelo presenta muchas limitaciones y es sólo una aproximación simplista sobre la difusión de calor. Las limitaciones
#' del modelo son las siguientes:
#' 
#' * se trabaja 
#' considerando una dimensión unidimensional determinando ecuaciones ordinarias lineales y 
#' evitando considerar el problema
#' de difusión de calor en todas las dimensiones que requeriría resolución numérica 
#' (por ej. usando diferencias finitas),
#' * el coeficiente k como constante sin contemplar su variación con la temperatura
#' * no se considera la influencia de la humedad en el medio
#' * consideramos a la capa superficie de suelo muy delgada como para que su volumen no influya en el problema, y 
#' por lo tanto, aplanamos el problema simplificandolo a una dimensión considerando la superficie del suelo como 
#' una fina capa entre dos rectángulos, como describe la figura.
#' 
#' Para la realización de este laboratorio se ha utilizado R (http://www.r-project.org/), un software open source 
#' con funcionalidades similares a Matlab.
#' Utilizaremos la librería *deSolve* que permite resolver ecuaciones diferenciales. Como parámetro uno puede pasarle el método de resolución numérica por ejemplo Euler, Runge Kutta, etc.
#' 
#' 
#' Constantes y condiciones iniciales
#' ---
#' 
#' Determinamos los coeficientes de disminución de temperaturas las constantes del problema
c_p <- -0.010 
c_s <- -0.014
c_a <- -0.018

#' Declaramos la constante de difusión calórica dependiente del material para suelo arenoso húmedo
k <- 0.017
#' Declaramos algún valor de emisividad para el suelo
epsilon <- 0.002
#' área de la subsuelo: 1 m profundiad por 5 m de ancho o largo de la superficie
area_subsuelo <- 1 * 5
#' área de la superficie: 5 m de largo por uno de ancho
area_suelo <- 1 * 5
#' Definimos la constante de Stefan-Boltzmann
sigma <- 5.67 * 10^-8
#' Definimos la constante de convección del aire (unidad W/(m2*K))
h <- 2
parameters <- c(c_p=c_p, c_s=c_s, c_a=c_a, k=k, epsilon=epsilon, 
                area_subsuelo=area_subsuelo, area_suelo=area_suelo, sigma=sigma)
#' temperatura profunda a 25 °C pero en Kelvin es 290
T_p_i <-290 
#' temperatura de superficie a 29 °C pero en Kelvin es 280
T_s_i <-280
#' temperatura de aire a 35 °C pero en Kelvin es 278
T_a_i <- 278

Q_s <-  -k * area_subsuelo * (T_s_i - T_p_i)
Q_r <- epsilon * area_suelo * sigma * (T_s_i^4)
Q_c <- h * area_suelo * (T_a_i - T_s_i )
Q <- Q_s - Q_r + Q_c

state <- c(T_p=T_p_i, T_s=T_s_i, T_a=T_a_i, Q_s=Q_s, Q_r=Q_r,Q_c=Q_c, Q=Q) # estados iniciales

TOTAL_TIME <<- 10  # simulamos 10 horas para ir observando
times <- seq(1, TOTAL_TIME, by = 1)

#' 
#' Resolución del modelo 
#' ---
#' 
#' Cargamos la librería *deSolve* que permite resolver ecuaciones diferenciales. 
#' Como parámetro uno puede pasarle el método de resolución numérica por ejemplo Euler, Runge Kutta, etc.
#' 

library(deSolve)
# 
DifusionCalor<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dT_p <- c_p * T_p
    dT_s <- c_s * T_s
    dT_a <- c_a * T_a 
    Q_s <-  -k * area_subsuelo * (T_s - T_p)
    Q_r <- epsilon * area_suelo * sigma * (T_s^4)
    Q_c <- h * area_suelo * (T_a - T_s )
    Q <- Q_s - Q_r + Q_c
    
    list(c(dT_p,dT_s,dT_a,Q_s,Q_r,Q_c,Q))
  }) 
}
#' Llamamos a la función *ode* de la librería *deSolver*, que resuelve ecuaciones diferenciales ordinarias
resultado <- ode(y = state, times = times, 
                 func = DifusionCalor, parms = parameters, 
                 method="euler")

#' Mostramos los resultados del cálculo realizado, utilizando redondeo a dos decimales.
print(round(resultado,2))
#' Graficamos la variación de temperatura en el tiempo, convirtiendo los grados Kelvin a °C
library(graphics)
plot.new()
matplot(resultado[,2:4]-275, type = c("l"),pch=1,col=2:4,xlab = "horas",ylab="°C") 
legend("topright",legend=colnames(resultado[,2:4]),col=2:4, pch=1) 

#' El comando anterior construye la gráfica y le agrega una leyenda o referencia
#' Graficamos todos los flujos Q
matplot(resultado[,5:8],type=c("l"),pch=1,col=2:6,xlab ="horas",ylab = "Energía emitida W/m2")
legend("topleft",legend=colnames(resultado[,5:8]), col=2:6, pch=1)
#' Graficamos individualmente cada uno de los flujos 
matplot(resultado[,5], 
        type = c("l"),
        pch=1,col = 6,
        xlab = "horas", 
        ylab = paste(colnames(resultado)[5]," W/m2"))

matplot(resultado[,6], type = c("l"),pch=1,col = 5,xlab = "horas", 
        ylab = paste(colnames(resultado)[6]," W/m2"))

matplot(resultado[,7], type = c("l"),pch=1,col = 4,xlab = "horas", 
        ylab = paste(colnames(resultado)[7]," W/m2"))

#' 
#' Conclusiones
#' ---
#' 
#' En los resultados podemos notar la liberación de calor de los tres flujos y la principal incidencia
#' de la convección, que presenta una gran influencia en la pérdida de calor de la superficie del suelo.
#' Para realizar un modelo más realista habría que considerar un modelo tridimensional y por consiguiente 
#' utilizar otros métodos de resolución numéricos para calcular el gradiente de diferencias de temperaturas.
#' Otro aspecto es la turbulencia que no ha sido considerada en este problema.
#' 
#'  REFERENCIAS
#'  
#' [1] Soetaert, K., Meysman, F., & Petzoldt, T. (2010). 
#' Solving differential equations in R. AIP Conference Proceedings, 1281(December),
#'  31–34. doi:10.1063/1.3498463 
#'  
#'  
#' [2] Seinfeld J. H. and Pandis S. N. “Atmospheric Chemistry
#'   and Physics from Air Pollution to Climate Change” pag. 801.
#'   
#'   

