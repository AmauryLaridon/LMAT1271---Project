################################################################################
#                                                                              #
#                              LMAT 1271                                       #
#                           Projet 2022-2023                                   #
#         Auteurs : Marie Determe, Augustin Lambotte, Amaury Laridon           #
#         Script R utiliser pour le projet instructions et ressources          #
#       disponibles à : https://github.com/AmauryLaridon/LMAT1271-Project      #
#                                                                              #
################################################################################

rm(list=ls(all=TRUE)) # Nettoyage mémoire R Studio
set.seed(2023)

######################### - Partie 1.3 Simulations - ###########################

n <- 20 # size of the sampling
nbr_replication <- 100 # number of replication of the sampling process
a_0 <- 0.5 # alpha parameter of the density function of T
b_0 <- 1.2 # beta parameter of the density function of T

#Density function of random variable T
f_a_0b_0 <- function(t) {
  if(t > 0 & t < 1) {
    return(a_0*b_0*t^(a_0-1)*(1-t^a_0)^(b_0-1))
  } else {
    return(0)
  }
}

#Quantile of order p, a and b are unknown parameters
q_p_a_b <- function(p,a,b) {
  return ((1-(1-p)^(1/b))^(1/a))
}

#Quantile of order p
q_p <- function(p) {
  return (q_p_a_b(p,a_0,b_0))
}

#First estimator of q_0.95
q_s <- function(T_i) {
  size <- length(T_i)
  T_sort <- sort(T_i)
  return ((T_sort[ceiling(0.95*size)]+T_sort[floor(0.95*size)])/2)
}

#Estimator of a_m and b_m, using the method of moments
q_m <- function(T_i) {
  #On cherche une racine à la première équation de 1.2.b pour obtenir 
  #a_M, puis on calcule b_M et q_M
  f1 <- function(a) {
    return ((sum(T_i^(2*a)))*(n+sum(T_i^(a)))*(sum(T_i^(a)))^(-2)-2)
  }
  solution <- uniroot(f1,c(1e-9,10))
  a_m <- solution$root
  b_m <- n/(sum(T_i^(a_m)))-1
  return (q_p_a_b(0.95,a_m,b_m))
}

#Estimator of a_l and b_l, using the maximum likelihood estimator
q_l <- function(T_i) {
  #On cherche une racine à la première équation de 1.2.c pour obtenir a_L
  #puis on calcule b_L et q_L
  f1 <- function(a) {
    return (1/a + (1/(sum(log(1-T_i^a)))+(1/n))*sum((T_i^(a)
                                                     *log(T_i))/(1-T_i^(a)))+(1/n)*sum(log(T_i)))
  }
  solution_l <- uniroot(f1,c(1e-9,10))
  a_l <- solution_l$root
  b_l <- -n/(sum(log(1-T_i^(a_l))))
  
  return (q_p_a_b(0.95,a_l,b_l))
}

### - 1.3 (a) - ###

#Generation of a iid sample of size 20 from the density f_a_0b_0 and 
#estimation of q_0.95 in 3 different ways
generate_sample <- function() {
  U <- runif(n)
  T_i <- sapply(U,q_p)
  q_s_sample <- q_s(T_i)
  q_m_sample <- q_m(T_i)
  q_l_sample <- q_l(T_i)
  return (c(q_s_sample,q_m_sample,q_l_sample))
}

### - 1.3 (b) - ###

generate_estimation_q_p <- function(affichage) {
  estimation_q_p <- replicate(nbr_replication,generate_sample())
  estimation_q_s <- estimation_q_p[1,]
  estimation_q_m <- estimation_q_p[2,]
  estimation_q_l <- estimation_q_p[3,]
  
  
  ### - 1.3 (c) - ###
  
  #Approximation de l'espérance par la moyenne empirique
  esperance_q_s <- mean(estimation_q_s)
  esperance_q_m <- mean(estimation_q_m)
  esperance_q_l <- mean(estimation_q_l)
  
  #Calcul du biais
  bias_s <- esperance_q_s - q_p_true
  bias_m <- esperance_q_m - q_p_true
  bias_l <- esperance_q_l - q_p_true
  
  #Calcul de la variance
  variance_s <- mean(estimation_q_s^2)-esperance_q_s^2
  variance_m <- mean(estimation_q_m^2)-esperance_q_m^2
  variance_l <- mean(estimation_q_l^2)-esperance_q_l^2
  
  #Calcul du MSE
  mse_s <- bias_s^2 + variance_s
  mse_m <- bias_m^2 + variance_m
  mse_l <- bias_l^2 + variance_l
  
  if (affichage == TRUE) {
    hist(estimation_q_s,breaks=10, main="")
    title(main = paste("Histogram of estimated q_s\n", "n =", n, " N =", 
                    nbr_replication, " a_0 =", a_0, " b_0 =", b_0), 
                    sub= paste("Biais =", round(bias_s, digits=4), 
                    " Variance=", round(variance_s, digits=4), " MSE=",
                    round(mse_s, digits=4)))
    boxplot(estimation_q_s, xlab = "q_s", ylab="Value")
    title(main = paste("Boxplot of estimated q_s\n","n =", n, " N =", 
                    nbr_replication, " a_0 =", a_0, " b_0 =", b_0),
                    sub= paste("Biais =", round(bias_s, digits=4), 
                    " Variance=", round(variance_s, digits=4), " MSE=",
                    round(mse_s, digits=4)))
    
    
    hist(estimation_q_m,breaks=10, main = "")
    title(main = paste("Histogram of estimated q_m\n","n =", n, " N =", 
                    nbr_replication, " a_0 =", a_0, " b_0 =", b_0), 
                    sub= paste("Biais =", round(bias_m, digits=4), 
                    " Variance=", round(variance_m, digits=4), " MSE=",
                    round(mse_m, digits=4)))
    boxplot(estimation_q_m, xlab = "q_m", main = "", ylab="Value")
    title(main = paste("Boxplot of estimated q_m\n","n =", n, " N =", 
                    nbr_replication, " a_0 =", a_0, " b_0 =", b_0), 
                    sub= paste("Biais =", round(bias_m, digits=4), 
                    " Variance=", round(variance_m, digits=4), " MSE=",
                    round(mse_m, digits=4)))
    
    
    hist(estimation_q_l,breaks=10, main = "")
    title(main = paste("Histogram of estimated q_l\n","n =", n, " N =", 
                    nbr_replication, " a_0 =", a_0, " b_0 =", b_0), 
                    sub= paste("Biais =", round(bias_l, digits=4), 
                     " Variance=", round(variance_l, digits=4), " MSE=",
                     round(mse_l, digits=4)))
    boxplot(estimation_q_l, xlab = "q_l" , main = "", ylab="Value")   
    title(main = paste("Boxplot of estimated q_l\n","n =", n, " N =", 
                    nbr_replication, " a_0 =", a_0, " b_0 =", b_0), 
                    sub= paste("Biais =", round(bias_l, digits=4), 
                    " Variance=", round(variance_l, digits=4), " MSE=",
                    round(mse_l, digits=4)))
    
  }
  
  return (c(c(bias_s,bias_m,bias_l),
            c(variance_s,variance_m,variance_l),c(mse_s,mse_m,mse_l)))
}

#################### - Execution - ######################

#Vraie valeur de q_0.95
q_p_true <- q_p(0.95)

#1.3.a
generate_sample()

#1.3.b et 1.3.c
generate_estimation_q_p(TRUE)

#1.3.d
n_value <- c(20,50,100,250,400)
results <- array(0, dim = c(3, 3, length(n_value)))
for (i in 1:length(n_value)) {
  n<-n_value[i]
  results[,,i]=t(generate_estimation_q_p(FALSE))
}

results
#Graphiques

print(results[3,,])

#Biais
plot(n_value,results[1,1,],xlab="n",ylab="Bias",main="Bias of q_s") #q_s
plot(n_value,results[2,1,],xlab="n",ylab="Bias",main="Bias of q_m") #q_m
plot(n_value,results[3,1,],xlab="n",ylab="Bias",main="Bias of q_l") #q_l

#Variance
plot(n_value,results[1,2,],xlab="n",ylab="Variance",main="Variance of q_s") #q_s
plot(n_value,results[2,2,],xlab="n",ylab="Variance",main="Variance of q_m") #q_m
plot(n_value,results[3,2,],xlab="n",ylab="Variance",main="Variance of q_l") #q_l

#MSE
plot(n_value,results[1,3,],xlab="n",ylab="MSE",main="MSE of q_s") #q_s
plot(n_value,results[2,3,],xlab="n",ylab="MSE",main="MSE of q_m") #q_m
plot(n_value,results[3,3,],xlab="n",ylab="MSE",main="MSE of q_l") #q_l

n_values <- c(20,100,400)
for (i in 1:3) {
  n <- n_values[i]
  factor <- sqrt(n)
  estimation_q_p <- replicate(nbr_replication,generate_sample())
  estimation_q_l <- estimation_q_p[3,]
  hist(factor*(estimation_q_l-q_p_true),main=paste("Histogram of n^(1/2)(q_l-q_0.95) for n=",n),xlabel="n^(1/2)(q_l-q_0.95)")
}

################################################################################
#                                                                              #
#                              LMAT 1271                                       #
#                           Projet 2022-2023                                   #
#         Auteurs : Marie Determe, Augustin Lambotte, Amaury Laridon           #
#         Script R utiliser pour le projet instructions et ressources          #
#       disponibles à : https://github.com/AmauryLaridon/LMAT1271-Project      #
#                                                                              #
################################################################################

rm(list=ls(all=TRUE)) # Nettoyage mémoire R Studio
set.seed(2023) # génération clef aléatoire

##### Chargement librairies #####
#install.packages("lattice")
#install.packages("ggplot2")
#install.packages("esquisse")
library(lattice)
library(ggplot2)
#library(esquisse)

######################### - Partie 2 Regression - ###########################

### Importation des données ### 

dataproject <- read.csv("/home/amaury/Bureau/LMAT1271 - Calcul des probabilités et analyse statistique/Projet/LMAT1271-Project/Data/dataproject.txt", sep=";")
#View(dataproject)
fuel <- dataproject$Y
hrspw <- dataproject$X
n <- length(fuel)

### Question (a) ###

mean_hp <- mean(hrspw)
mean_fuel <- mean(fuel)


## Calcul de S_xx ##

S_xx <- sum((hrspw - mean_hp)^2)

## Calcul de S_xy ##

S_xy <- sum((hrspw-mean_hp)*(fuel-mean_fuel))

## Calcul de S_yy ##

S_yy <- sum((fuel-mean_fuel)^2)

## Calcul des estimateurs de beta_0 et beta_1 à l'aide de leur définition par la méthode MCO ##

beta_1 <- S_xy/S_xx
beta_0 <- mean_fuel - beta_1*mean_hp 
print(beta_0)
print(beta_1)

## Calcul de l'estimateur de la variance des erreurs ## 

var_err <- (S_yy-(beta_1*S_xy))/(n-2)
std_err <- sqrt(var_err)

## Création du modèle linéaire sur base des paramètres estimés ## 

lm <- beta_0+(beta_1*hrspw)

## Calcul des estimateurs de la variance de beta_0 et beta_1 à l'aide de leur définition par la méthode MCO ##

std_beta_1 <- var_err/S_xx

## Calcul de la statique de test T_obs ## 

T_obs <- ((beta_1 - 0)/(std_err))*sqrt(S_xx)

## Prévision ## 

X.star <- 1000
Y.star <- beta_0 + beta_1*X.star

## Intervalle de prévision à 95% ##

alpha <- 0.05
t <- qt(p=alpha/2, df=n-2, lower.tail = F)

IP_plus <- beta_0+(beta_1*X.star)+(t*std_err)*sqrt(1+(1/n)+((X.star-mean_hp)^2)/S_xx)
IP_min <- beta_0+(beta_1*X.star)-(t*std_err)*sqrt(1+(1/n)+((X.star-mean_hp)^2)/S_xx)

## Calcul p valeur ##

p_value <- 2*(pt(-T_obs,df=98))
print(p_value)

## Affichage ## 

xyplot(fuel~hrspw,xlab="Puissance [hP]", ylab="Ln(Consommation d'essence) [l/100km]")
xyplot((beta_0+(beta_1*hrspw))~hrspw)

## Verification à l'aide de la fonction déjà implémentée lm() ##

lm1 <- lm(formula = fuel~hrspw, data = dataproject) #Calcul et création d'un modèle de regression
                                   #linéaire
summary(lm1)

model <- lm(fuel ~ hrspw)
new <- data.frame(hrspw=1000)
predict(lm1,new, interval="confidence") 



