#' ---
#' title: "Modelo de balance calórico de helada radiativa"
#' author: "Ana Laura Diedrichs"
#' date: "febrero, 2015"
#' output: pdf_document
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
#'(primera ley de conducción del calor): $$Q_G=-k* Area * \nabla T $$. La constante $k$ es la conductividad térmica del suelo.
#'El gradiente de temperatura que representa el flujo de calor entre suelo y ambiente se caracteriza como la siguiente fórmula 
#'$\nabla T =\frac{dT}{dt}=\frac{dT}{dx}+\frac{dT}{dy}$. Dada la presencia de derivadas parciales, la misma se resolverá 
#'mediante diferencias finitas y considerando dimensiones 2D para simplificar el problema.
#'
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
#' masa ni capacidad calorífica, de manera que los flujos que entran y salen de calor de esta superficie se 
#' compensen según la ley de conservación de energía.
#' * la temperatura del suelo variará linealmente con el tiempo
#' * se simularán 100 unidades de tiempo
#' * Viento y humedad serán constantes en el tiempo.
#' * Se brinda variación de la temperatura del suelo $T_G$ y ambiental $T$
#' 
#' 
#' OBJETIVOS DEL MODELO 
#' 
#' Evaluar la evolución temporal de $Q_G, SH, LH, Q^*$ y $S_C$ ante la variación de temperatura del aire
#' y el suelo en diferentes casos (por ejemplo distintos tipos de suelo, distintas variaciones de temperatura)

#library(deSolve)

parameters <- c(e=0.7) # <--- e emisividad, k constante termica, h humedad
state <- c(P_t = 10) # <--- población inicial de 10 personas
times <- seq(0, 10, by = 1) # <--- simulamos tan sólo 10 tiempos

BalanceCalorico<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    P_t1 = k * P_t
    
    list(c(P_t1))
  }) # end with(as.list ...
}
#'
#' Código extraído del pdf titulado como
#' Soetaert, K., Meysman, F., & Petzoldt, T. (2010). 
#' Solving differential equations in R. AIP Conference Proceedings, 1281(December),
#'  31–34. doi:10.1063/1.3498463
library(ReacTran)
Grid <- setup.grid.1D(N=1000,L=10)

pde1D <-function(t, C, parms) {
  tran <- tran.1D(C = C, D = D,
                  C.down = Cext, dx = Grid)$dC
  list(tran - Q) # return value: rate of change
}

D<-1
Q<-1
Cext<-20

library(rootSolve)

print(system.time(
  std<-steady.1D(y = runif(Grid$N),
                 func = pde1D, parms = NULL, nspec = 1)
))

plot (Grid$x.mid, std$y, type = "l", lwd = 2, main = "steady-state PDE", xlab = "x", ylab = "C", col = "red")

require(deSolve)
times <- seq(0, 100, by = 1)
system.time(
  out <- ode.1D(y = rep(1, Grid$N),
                times = times, func = pde1D,
                parms = NULL, nspec = 1)
)

image(out, xlab = "time, days",ylab = "Distance, cm",main = "PDE", add.contour=TRUE)

