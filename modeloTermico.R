#' ---
#' title: "Difusión del calor en noches de helada: experimentos"
#' author: "Ana Laura Diedrichs"
#' date: "junio, 2015"
#' ---
#' 
#' 
#' OBJETIVOS DEL MODELO 
#' ---
#' 
#' Evaluar la variación de la temperatura del suelo de acuerdo a la difusión del calor desde el subsuelo al suelo,
#' presentado como: $Q_G=-k* Area * \nabla T $. Para empezar a codificar y simular el comportamiento usando R , utilizamos la 
#' librería ReacTran [1]. Esta librería cuenta con rutinas para el desarrollo de modelos de reacción y transporte, con advección 
#' y difusión en una, dos o tres dimensiones. Principalmente provee:
#' 
#' * Funciones que subdividen el espacio en un número discreto de celdas en una grilla
#' * Aproximación por diferencias finitas o volúmenes finitos del término de transporte advectivo-difusivo
#'
#' $$ \frac{\partial T}{\partial t}
#' = -\frac{1}{A_x\xi_x} \left( \frac{\partial}{\partial x} A_x
#'           ( -D \frac{\partial \xi_x T}{\partial x} ) - 
#'            + \frac{\partial}{\partial x} (A_x\, v\,  \xi_x\, T) \right) $$
#'
#'En la fórmula anterior D es el coeficiente de difusión, $v$ la tasa de advección, $A_x$ el área de la superficie y
#'$\xi$ la fracción de volumen. Considerando que $A$,$\xi$,$D$ y $v$ son constantes en x, la fórmula podría ser reescrita como:
#'
#' $$ \frac{\partial T}{\partial t} =  -D \frac{\partial^2 T}{\partial x^2}  - u
#'            \frac{\partial T}{\partial x}  $$
#'
#' Podemos observar que en el primer término tenemos una derivada de segundo orden y el segundo de primer orden. El
#' primer término representa la tasa de difusión y el segundo término la tasa de consumo.
#' Vamos a simular la propagación del calor desde el subsuelo al suelo basados en la siguiente fórmula:
#' $$Q_G=-k* Area * \nabla T$$
#' 
#' Declaramos la constante de difusión calórica dependiente del material para suelo arenoso húmedo
k <- 0.017
#' área de la hsuperficie
area <- 1* Grid$N
D <- k * area # k * Area
#' Determinamos un consumo constante valor uno, al igual que la condición de contorno superior
Q <- 1
Cext <- 1
C.upp <- 0.001

# La siguiente línea carga a memoria la librería ReacTran.
library(ReacTran)

#' Utilizo una grilla unidimensional de propagación de calor
Grid <- setup.grid.1D(x.down=2,x.up=-1,N=1000,L=10)

pde1D <-function(t, C, parms) {
  tran <- tran.1D(C = C, D = D,
                  C.up = Cext,
                  dx = Grid)$dC
  list(tran-Q) # return value: rate of change
}

library(rootSolve)

std<-steady.1D(y =rep(0,Grid$N), func = pde1D, parms = NULL, nspec = 1)


plot(Grid$x.mid, std$y, type = "l", lwd = 2, main = "steady-state PDE", xlab = "x", ylab = "T", col = "red")

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

#' CONCLUSIONES MODELO DIFUSIÓN CALOR EN SUELO
#' -------
#' Para distintos valores de la constante de conductividad, propia de las características del material, observamos que mientras
#' más bajo sea el valor de k más resistencia ofrece para conducir calor de un extremo al otro.
#' 
#' MODELO DE DIFUSIÓN DEL CALOR EN EL AIRE
#' -------
#' La difusión del aire en función de la temperatura está dada por la siguiente fórmula [2,3], en unidades mts, siendo
#' $T$ temperatura y $P$ presión 
#' 
#' $$ D_v(T)= \frac{0.211}{P} (\frac{T}{273})^(1.94) [\frac{m^2}{s}]$$
#' 
#' $$ \frac{\partial T}{\partial t} =  -D \frac{\partial^2 T}{\partial x^2}  - u
#'            \frac{\partial T}{\partial x}  $$
#'            
#'            
#' A continuación declaramos una función la difusión del aire respecto a su temperatura. Para cal
#' Regresa la temperatura en sistema mts
#' **REVISAR PAPER PROFESOR RAUL PEREZ** 
#' 
#' 
#' Grilla vertical 
N_size <- 1000
Grid <- setup.grid.1D(x.down=2,x.up=-1,N=N_size,L=10)

#' Tomamos archivo de radiosondeo
mendoza.16.7.2015 <- read.csv("~/phd-repos/cursos/fisica/fisica-atmosferica/radiosondeos/mendoza-16-7-2015.csv")
alturas <- as.data.frame(cbind(mendoza.16.7.2015[1:91,2],mendoza.16.7.2015[2:92,2]))
#' index: indice en la grilla vertical discretizada
#' x1 vector de variables de alturas dado dataset
#' y1 vector de valores de presión o temperatura
#utilizamos las primeras 10 alturas, y por lo tanto, los primeros 10 valores de presión o temperatura
h.min <- mendoza.16.7.2015[1,"HGHT_m"]
h.max <- mendoza.16.7.2015[10,"HGHT_m"]
escalon <- (h.max - h.min)/1000
alturas <- mendoza.16.7.2015[1:10,"HGHT_m"]

ecuacion.recta <- function(index.N, y1)
{
  h <- escalon * index.N
#  which()
}

#' Función de difusión del aire respecto a la temperatura
#' pasamos t: tiempo o vector índice
#' h: altura o indice del vector de altura a considerar
#' temp: matriz de temperaturas vs t y h
#' pression: matriz de presion vs t y h
#' Regresa un vector con los valores de difusión **VERIFICAR**
difussion.air <- function(t,h,temp, pression)
{
 return (0.211/pression[t,h]) * (temp[t,h]/273)^194 
}
#' Where T is the temperature (Kelvin scale), and p the atmospheric pressure (Pa).
#' Valores iniciales del problema
#' 

#' 
Q <- 1
Cext <- 1
C.upp <- 0.001

#' Simulamos 12 horas (720 minutos) simulamos
times <- seq(0, 720, by = 1) 
l <- length(times)
#' Guardaremos los valores de temperatura y presión
#' para calcular la función de difusión
#' 
temperature <- seq(from=-2,to=-12,length.out=Grid$N)
pression <- 1000 #as.vector(runif(Grid$N,min=1000,max=1200))

#'Calculamos la difución en sistema métrico mts
#'considerando la difusión del aire
#'
pde1D <-function(t, C, parms) {
  
  D <- (0.211/pression) * (temperature[1]/273)^194 
  
  tran <- tran.1D(C = C, D = D,
                  C.up = Cext,
                  dx = Grid)$dC
  list(tran-Q) # return value: rate of change
}


std<-steady.1D(y =rep(0,Grid$N), func = pde1D, parms = NULL, nspec = 1)

#plot(Grid$x.mid, std$y, type = "l", lwd = 2, 
#     main = "steady-state PDE", xlab = "x", ylab = "T", col = "red")


out <- ode.1D(y = rep(1, Grid$N),
              times = times, func = pde1D,
              parms = NULL, nspec = 1)
#' En el siguiente gráfico 
#' mostramos la variación de la temperatura en el tiempo a distintas
#' alturas 
image(out, xlab = "minutes",ylab = "Distance, m",main = "Difusión Calor en Aire", add.contour=TRUE)


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

