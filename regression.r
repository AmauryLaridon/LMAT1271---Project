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












