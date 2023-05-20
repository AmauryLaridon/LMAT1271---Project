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
install.packages("esquisse")
library(lattice)
library(ggplot2)
library(esquisse)

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

S_xx <- 0
for (i in range(length(hrspw))){
  S_xx = S_xx + (hrspw[i]-mean_hp)^2
}

## Calcul de S_xy ##

S_xy <- 0
for (i in range(length(fuel))){
  S_xy = S_xy + ((hrspw[i]-mean_hp)*(fuel[i]-mean_fuel))
}

## Calcul de S_yy ##

S_yy <- 0
for (i in range(length(fuel))){
  S_yy = S_yy + (fuel[i]-mean_fuel)^2
}

## Calcul des estimateurs de beta_0 et beta_1 à l'aide de leur définition par la méthode MCO ##

beta_1 <- S_xy/S_xx
beta_0 <- mean_fuel - beta_1*mean_hp 
print(beta_0)
print(beta_1)

## Calcul de l'estimateur de la variance des erreurs ## 

var_err <- (S_yy-(beta_1*S_xx))/(n-2)

## Création du modèle linéaire sur base des paramètres estimés ## 

lm <- beta_0+(beta_1*hrspw)

## Affichage ## 

xyplot(fuel~hrspw,xlab="Puissance [hP]", ylab="Ln(Consommation d'essence) [l/100km]")
xyplot((beta_0+(beta_1*hrspw))~hrspw)

## Verification à l'aide de la fonction déjà implémentée lm() ##

lm1 <- lm(Y~X, data = dataproject) #Calcul et création d'un modèle de regression
                                   #linéaire
                                  
model <- lm(fuel ~ hrspw)
#new <- data.frame(hrspw=1000)
#predict(model,new, interval="confidence") 












