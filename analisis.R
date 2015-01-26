
#' ---
#' title: "Principal Component Analysis of 6 vertical temperature sensors"
#' author: "Ana Laura Diedrichs"
#' date: "May 3rd, 2014"
#' output: pdf_document
#' ---
#' PRINCIPAL COMPONENT ANALISIS
#' 
#' 
#' librería a utilizar. Constantes para cargar el archivo de datos

library("FactoMineR")

SEP <- ";"
FILE_NAME <- "120726-minimal.csv"

#+COLNAMES <- c("s1","s2","s3","s4","s5","s6","s1_t","s2_t","s3_t","s4_t","s5_t","s6_t")
#+TOTAL_SENSORS<-6

#' 
#' Cargamos el dataset a memoria. Realizamos un vistazo de su contenido. Contiene 6 variables.
#' 


data <- as.data.frame(read.csv(FILE_NAME,sep=SEP))  
dataset <- data[,5:10]
head(dataset)
#' 
#' Renombramos las columnas para mayor claridad
#' 

names.col <- c("time","s_0","s_0_4","s_0_75","s_1_50","s_2","s_3")
data.table <- as.data.frame(cbind(paste(data[,2],data[,3],sep=" "),data[,5:10]))
colnames(data.table)<- names.col
#' 
#' Crearemos un objeto timeSeries (xts) para manipular mejor los datos
#' 

library(lubridate)
t <- ymd_hms(data.table[,1])  #<--convierte string a date
library(xts)
ts<- xts(x=data.table[,2:7],order.by=t)
my.color<- c("blue","red","yellow","green","orange","cyan")

plot(as.zoo(ts), plot.type="s", col=my.color, lty=1,lwd=2,ylab="temperatures")
legend(x="topleft", legend=names.col, col=my.color,lty=1:6) 
#TODO verificar que la leyenda del gráfico corresponda con las variables


#'
#' A continuación graficamos la serie temporal de los sensores por separado para mayor claridad
#' 

 require(graphics)
library(timeSeries)
par(mfrow=c(1, 1))
#línea que marca el CERO
lines2 <- function(X, Y, type, xlab, ylab, col, pch, lty, lwd, cex) {
  lines(x=X, y=Y, col=col)
  abline(h=0, col = "brown", lwd=2)
}
plot(as.zoo(ts), plot.type="m", col = .colorwheelPalette(3),panel=lines2)
#' 
#' Gráfico de cajas para tener un vistazo de los datos
#' 
boxplot(data.table[,2:7])
#' ---
#' ANEXO ANALISIS DE COMPONENTES PRINCIPALES 
#' 
#' principal component analysis
#' ---

#' metodo para analisis de componentes principales 
fit <- princomp(dataset, cor=TRUE)
# print variance accounted for
summary(fit) 
# pc loadings
loadings(fit) 
plot(fit,type="lines") # scree plot
#' fit$scores # the principal components
biplot(fit)

# PCA Variable Factor Map
result <- PCA(dataset) # graphs generated automatically 


