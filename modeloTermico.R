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
#'$$Q_G + SH + LH + Q^*=0$$
#'
#'En las noches de heladas el ambiente se enfría, ocasionando que el suelo pierda el calor adquirido durante
#'el día gracias al sol. El flujo por conducción en el suelo ($Q_G$) es el flujo de energía desde el subsuelo a la superficie
#'(positivo) y viceversa (negativo). Este flujo de calor por conducción del suelo se puede cuantificar
#'usando la ley de Fourier, que indica que el calor $Q_G$ es proporcional al gradiente de temperatura
#'(primera ley de conducción del calor): $Q_G=-k* Area * \Delta T $. 
#'La constante $k$ es la conductividad térmica 
#'del suelo.
#'
#'El flujo de calor sensible $SH$ está asociado a la diferencia de temperatura que existe entre la temperatura
#'a nivel de superficie y la temperatura a nivel de aire. Su valor es positivo si la energía fluye 
#'desde el aire hacia la superficie y negativo en caso contrario. Se calcula como 
#'$$SH=C_H*\overline{W}*(\overline{T}-T_G)$$, siendo
#'$W$ la velocidad promedio del viento, $T$ la temperatura a nivel del aire, 
#'$T_G$ la temperatura a nivel del suelo, y $C_H$ es un coeficiente que en condiciones de estabilidad
#'puede variar entre $10^{-3}$ a $5*10^{-3}$
#'
#'El flujo de calor latente ($LH$) está asociado con los cambios de fase del agua, es decir, cuando hay 
#'evaporación o condensación sobre la superficie del suelo. Se calcula como
#'$$LH=C_E*\overline{W}*(\overline{q}-q_G)$$, siendo $q$ la humedad del aire, $q_G$ la humedad relativa a
#'nivel del suelo, y $C_E$ un coeficiente de estabilidad similar a $C_H$.
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
#'$LW \uparrow = \epsilon * \sigma * T^4 $, donde $\epsilon$ es el coeficiente de emisividad del suelo, 
#'$\sigma = 5.67 * 10^{-8} \frac{J}{m^2 s K}$  la constante de Stefan-Boltzmann y T la temperatura en grados
#'Kelvin.
#'
#'
#' CONSIDERACIONES DEL MODELO
#'
#' * La capacidad térmica y emisividad del suelo no variará con el tiempo.
#' * Se considera una superficie delgada en la interface suelo - atmósfera de manera que no tenga
#' masa ni capacidad calorífica, de manera que los flujos que entran y salen de calor de esta superficie se 
#' compensen según la ley de conservación de energía.
#' * la temperatura del suelo variará linealmente con el tiempo
#' * se simularán 100 unidades de tiempo
#' * 
#' 

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

