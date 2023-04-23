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

### - 1.3 (a) - ###
n <- 20 # taille n de l'échantillon
a_0 <- 2 # valeur initiale arbitraire de alpha_0  
b_0 <- 3 # valeur initiale arbitraire de beta_0 

# définition de la pdf fonction f_{alpha, beta}
pdf_f <- function(t, alpha, beta) {
    if (t>0 & t<1) {
       f = alpha*beta*(t^(alpha-1))*((1-((t^alpha))^(beta-1)))
    }
    else {
       f = 0
    }
    return(f)
}

sample <- rep(0, n) # création d'un vecteur de longueur n vide dans lequel on va mettre de manière aléatoire des données générées à partir de f_{alpha,beta}

for (i in sample) {
   t = sample(x = -1:1:0.1, size = 1)
   print(t) 
}