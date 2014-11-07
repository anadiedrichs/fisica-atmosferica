
#' --
#' title: Principal Component Analysis of 6 vertical temperature sensors
#' author: Ana Laura Diedrichs
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

#' ---
#' Dataset loaded in memory. Un vistazo de su contenido. Contiene 6 variables
#' ---


data <- as.data.frame(read.csv(FILE_NAME,sep=SEP))  
dataset <- data[,5:10]
head(dataset)

#' ---
#' normalizo el dataset. Leí que es necesario escalarlo ¿es correcto esto?
#' ---

dataset <- ((max(dataset)-min(dataset))/max(dataset)) * dataset

#' ---
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

