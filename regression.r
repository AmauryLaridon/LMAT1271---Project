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

######################### - Partie 2 Regression - ###########################

### Importation des données ### 

dataproject <- read.csv("~/Bureau/LMAT1271 - Calcul des probabilités et analyse statistique/Projet/Data/dataproject.txt", sep=";")
#View(dataproject)

### Question (a) ###

lm1 <- lm(Y~X, data = dataproject) #Calcul et création d'un modèle de regression
                                   #linéaire
                                   
fuel <- dataproject$Y
hrspw <- dataproject$X
model <- lm(fuel ~ hrspw)
new <- data.frame(hrspw=1000)
predict(model,new, interval="confidence") 












