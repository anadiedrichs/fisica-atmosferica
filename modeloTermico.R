#' ---
#' title: "Difusión del calor en noches de helada"
#' author: "Ana Laura Diedrichs"
#' date: "junio, 2015"
#' ---
#' 
#' Introducción
#' ---
#' 
#' 
#' El descenso térmico que caracteriza una noche de helada está determinado por un balance calórico negativo
#' del entorno. Intervienen varias causas de pérdida del calor del sistema que denominaremos como flujos salientes
#' o entrantes cuya unidad es $W/m^2$ (watts por metro cuadrado).
#' 
#' Los elementos que componen el balance calórico están determinados por la ecuación fundamental que describe 
#' el proceso:
#' 
#'$$Q_G + SH + LH + Q^*+S_C=0$$
#'
#'En las noches de heladas el ambiente se enfría, ocasionando que el suelo pierda el calor adquirido durante
#'el día gracias al sol. El flujo por conducción en el suelo ($Q_G$) es el flujo de energía desde el subsuelo a la superficie
#'(positivo) y viceversa (negativo). Este flujo de calor por conducción del suelo se puede cuantificar
#'usando la ley de Fourier, que indica que el calor $Q_G$ es proporcional al gradiente de temperatura
#'(primera ley de conducción del calor): 
#'
#'$$Q_G=-k* Area * \nabla T $$. 
#'
#'La constante $k$ es la conductividad térmica del suelo.
#'El gradiente de temperatura que representa el flujo de calor entre suelo y ambiente se caracteriza como la siguiente fórmula 
#'$\nabla T =\frac{dT}{dt}=\frac{dT}{dx}+\frac{dT}{dy}$. Dada la presencia de derivadas parciales, la misma se resolverá 
#'mediante diferencias finitas y considerando dimensiones 2D para simplificar el problema.
#'
#'
#'El flujo de calor sensible $SH$ está asociado a la diferencia de temperatura que existe entre la temperatura
#'a nivel de superficie y la temperatura a nivel de aire. Su valor es positivo si la energía fluye 
#'desde el aire hacia la superficie y negativo en caso contrario. Se calcula como 
#'$SH=C_H*\overline{W}*(\overline{T}-T_G)$, siendo
#'$W$ la velocidad promedio del viento, $T$ la temperatura a nivel del aire, 
#'$T_G$ la temperatura a nivel del suelo, y $C_H$ es un coeficiente que en condiciones de estabilidad
#'puede variar entre $10^{-3}$ a $5*10^{-3}$. Se considera que $\overline{T}-T_G = \nabla T $, es decir,
#'reemplazar la expresión por el gradiente de temperatura quedaría determinada como
#'$$SH=C_H*\overline{W}*\nabla T$$
#'
#'El flujo de calor latente ($LH$) está asociado con los cambios de fase del agua, es decir, cuando hay 
#'evaporación o condensación sobre la superficie del suelo. Se calcula como
#'$$LH=C_E*\overline{W}*(\overline{q}-q_G)$$, siendo $q$ la humedad del aire, $q_G$ la humedad relativa a
#'nivel del suelo, y $C_E$ un coeficiente de estabilidad similar a $C_H$. Dado que esta ecuación representa 
#'el aporte/sustracción ante los cambios de estado del agua y como simulamos en un entorno de helada, si 
#'la temperatura es menor o igual a cero, se calculará $LH$; de lo contrario $LH=0$.
#'
#'El flujo por radiación neta ($Q^*$) considera la suma de la radiación emitida/recibida en los dos rangos
#'espectrales: radiación de onda corta (SW) y radiación de onda larga. Durante una noche de helada, dada
#'la ausencia de radiación solar, la ecuación se simplifica considerando sólo la radiación de onda larga.
#'$$ Q^*=LW \downarrow - LW \uparrow $$
#'
#'$LW \downarrow$ es el flujo de onda larga que el suelo recibe y es emitida casi en su totalidad por toda la 
#'atmósfera, dependiente de la temperatura, humedad y cobertura nubosa. Las nubes juegan
#'un papel muy importante durante la noche ya que modulan la temperatura de la superficie mediante la emisión
#'radiación infrarroja. 
#'
#'$LW \uparrow$ representa el flujo
#'de radiación de onda larga emitida desde la superficie terrestere hacia la atmósfera por lo que depende
#'de la temperatura del suelo y su emisividad. Interesa señalar que es responsable del calentamiento
#'y enfriamiento del aire, debido a que el aire no absorve la radiación solar por ser de onda corta, pero 
#'sí absorbe la del suelo que es onda larga. Está regida por la ley de Stefan-Boltzmann 
#'$LW \uparrow = \epsilon * \sigma * T_G^4 $, donde $\epsilon$ es el coeficiente de emisividad del suelo, 
#'$\sigma = 5.67 * 10^{-8} \frac{J}{m^2 s K}$  la constante de Stefan-Boltzmann y T_G la temperatura en grados
#'Kelvin del suelo.
#'
#'El flujo de calor por convección varía proporcionalmente según la constante de convección, definimos
#'$S_C = h * \nabla T$
#'
#' CONSIDERACIONES DEL MODELO
#'
#' * La capacidad térmica y emisividad del suelo no variará con el tiempo.
#' * Se considera una superficie delgada en la interface suelo - atmósfera de manera que no tenga
#' masa ni capacidad calorífica, para que los flujos que entran y salen de calor de esta superficie se 
#' compensen según la ley de conservación de energía.
#' 
#' Queda por considerar lo siguiente: 
#' 
#' * la temperatura del suelo variará linealmente con el tiempo
#' * se simularán 100 unidades de tiempo
#' * Viento y humedad serán constantes en el tiempo.
#' * Se brinda variación de la temperatura del suelo $T_G$ y ambiental $T$
#' 
#' 
#' OBJETIVOS DEL MODELO 
#' 
#' Evaluar la variación de la temperatura del suelo de acuerdo a la difusión del calor desde el subsuelo al suelo,
#' presentado como: $Q_G=-k* Area * \nabla T $. Para empezar a codificar y simular el comportamiento, utilizamos la 
#' librería ReacTran [1]
#'
#'  
library(ReacTran)

#' GRILLA DE DISCRETIZACION 2D
#' 
x.grid  <- setup.grid.1D(x.up = 0, L = 10, N = 100)
y.grid  <- setup.grid.1D(x.up = 0, L = 10, N = 100)
grid2D <- setup.grid.2D(x.grid, y.grid)

#' imprimimos los valores que guarda la variable grid2D seteados por defecto
#print(grid2D)

#' Vamos a simular la propagación del calor desde el subsuelo al suelo basados en la siguiente fórmula:
#' $$Q_G=-k* Area * \nabla T$$
#' Utilizo una grilla unidimensional de propagación de calor
#' 
#' Constante de difusión calórica dependiente del material para suelo arenoso húmedo
k <- 0.017
#' área de la hsuperficie
area <- 1* Grid$N
D <- k * area # k * Area
#' ¿Qué sería Q en mi problema?
Q <- 1
Cext <- 1
C.upp <- 0.001

Grid <- setup.grid.1D(x.down=2,x.up=-1,N=1000,L=10)

pde1D <-function(t, C, parms) {
  tran <- tran.1D(C = C, D = D,
                  C.up = Cext,
                  #C.up = Cext,
                  dx = Grid)$dC
  list(tran-Q) # return value: rate of change
}

library(rootSolve)

std<-steady.1D(y =rep(0,Grid$N),
                 func = pde1D, parms = NULL, nspec = 1)


plot(Grid$x.mid, std$y, type = "l", lwd = 2, 
     main = "steady-state PDE", xlab = "x", ylab = "T", col = "red")

#' Cargamos a memoria la librería que resuelve ecuaciones diferenciales.
require(deSolve)

#' Simulamos 12 horas (720 minutos) simulamos

times <- seq(0, 720, by = 1) 

out <- ode.1D(y = rep(1, Grid$N),
                times = times, func = pde1D,
                parms = NULL, nspec = 1)
#' En el siguiente gráfico 
#' mostramos la variación de la temperatura en el tiempo a distintas
#' profundidades 
image(out, xlab = "minutes",ylab = "Distance, m",main = "Suelo arenoso húmedo", add.contour=TRUE)
#' resolución utilizando diferencias finitas
#' sin usar la librería ReacTran

#'
#' Código ejemplo de la libreria ReacTran para 2D
#' 
#' Suelo arenoso seco
k <- 0.00017
D <- k * area # k * Area

out <- ode.1D(y = rep(1, Grid$N),
              times = times, func = pde1D,
              parms = NULL, nspec = 1)
#' En el siguiente gráfico mostramos la variación de la temperatura en el tiempo a distintas
#' profundidades 
image(out, xlab = "minutes",ylab = "Distance, m",main = "Suelo arenoso seco k=0.00017", add.contour=TRUE)

#' Alta conductividad en el suelo
k <- 0.4
D <- k * area # k * Area

out <- ode.1D(y = rep(1, Grid$N),
              times = times, func = pde1D,
              parms = NULL, nspec = 1)
#' En el siguiente gráfico mostramos la variación de la temperatura en el tiempo a distintas
#' profundidades 
image(out, xlab = "minutes",ylab = "Distance, m",main = "Suelo k=0.4", add.contour=TRUE)

#' MODELO DE DIFUSIÓN DEL CALOR EN EL AIRE
#' -------
#' A continuación declaramos una función para la 
#' la difusión del aire respecto a su temperatura
#' Regresa la temperatura en sistema mts
#' **REVISAR PAPER PROFESOR RAUL PEREZ** 
#' Función de difusión del aire respecto a la temperatura
#' pasamos t: tiempo o vector índice
#' h: altura o indice del vector de altura a considerar
#' temp: matriz de temperaturas vs t y h
#' pression: matriz de presion vs t y h
#' Regresa un vector con los valores de difusión **VERIFICAR**
difussion.air <- function(t,h,temp, pression)
{
 return (0.211/pression[t,]) * (temp[t,]/273)^194 
}
#' Valores iniciales del problema
#' 
#' 
Q <- 1
Cext <- 1
C.upp <- 0.001
#' Grilla vertical 
Grid <- setup.grid.1D(x.down=2,x.up=-1,N=1000,L=10)
#' Simulamos 12 horas (720 minutos) simulamos
times <- seq(0, 720, by = 1) 
l <- length(times)
temperature <- matrix(data=runif(Grid$N*l,min=-3,max=2),nrow=l,ncol=Grid$N)
pression <- matrix(data=runif(Grid$N*l,min=1000,max=1200),nrow=nrow(times),ncol=Grid$N)

my.parms <- NULL

pde1D <-function(t, C, parms) {
  #ARREGLAR EL TERMINO DE DIFUSSION
  
  D <- (0.211/pression[t,]) * (temp[t,]/273)^194 
  
  tran <- tran.1D(C = C, D = D,
                  C.up = Cext,
                  #C.up = Cext,
                  dx = Grid)$dC
  list(tran-Q) # return value: rate of change
}

std<-steady.1D(y =rep(0,Grid$N),
               func = pde1D, parms = my.parms, nspec = 1)

plot(Grid$x.mid, std$y, type = "l", lwd = 2, 
     main = "steady-state PDE", xlab = "x", ylab = "T", col = "red")


out <- ode.1D(y = rep(1, Grid$N),
              times = times, func = pde1D,
              parms = NULL, nspec = 1)
#' En el siguiente gráfico 
#' mostramos la variación de la temperatura en el tiempo a distintas
#' profundidades 
image(out, xlab = "minutes",ylab = "Distance, m",main = "Suelo arenoso húmedo", add.contour=TRUE)


#'  REFERENCIAS
#'  
#' [1] Código extraído del pdf titulado como
#' Soetaert, K., Meysman, F., & Petzoldt, T. (2010). 
#' Solving differential equations in R. AIP Conference Proceedings, 1281(December),
#'  31–34. doi:10.1063/1.3498463
#'  
