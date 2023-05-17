rm(list=ls(all=TRUE)) # Nettoyage mÃ©moire R Studio
set.seed(2023)

n <- 20
nbr_replication <- 100
a_0 <- 0.5
b_0 <- 1.2

#Fonction de densité
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
  return (T_sort[ceiling(0.95*size)]+(T_sort[floor(0.95*size)]-T_sort[ceiling(0.95*size)])/2)
}

#Estimator of a_m and b_m, using the method of moments
q_m <- function(T_i) {
  #On cherche une racine à la première équation de 1.2.b pour obtenir a_M, puis on calcule b_M et q_M
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
  #On cherche une racine à la première équation de 1.2.c pour obtenir a_L, puis on calcule b_L et q_L
  f1 <- function(a) {
    return (1/a + (1/(sum(log(1-T_i)))+1/n)*sum(T_i^(a)*log(T_i)/(1-T_i^(a)))+1/n*sum(log(T_i)))
  }
  solution_l <- uniroot(f1,c(1e-9,10))
  a_l <- solution_l$root
  b_l <- -n/(sum(log(1-T_i^(a_l))))
  
  return (q_p_a_b(0.95,a_l,b_l))
}

#Generation of a iid sample of size 20 from the density f_a_0b_0 and estimation of q_0.95 in 3 different ways
generate_sample <- function() {
  U <- runif(n)
  T_i <- sapply(U,q_p)
  q_s_sample <- q_s(T_i)
  q_m_sample <- q_m(T_i)
  q_l_sample <- q_l(T_i)
  return (c(q_s_sample,q_m_sample,q_l_sample))
}

generate_estimation_q_p <- function(affichage) {
  estimation_q_p <- replicate(nbr_replication,generate_sample())
  estimation_q_s <- estimation_q_p[1,]
  estimation_q_m <- estimation_q_p[2,]
  estimation_q_l <- estimation_q_p[3,]
  
  if (affichage == TRUE) {
    hist(estimation_q_s,breaks=10)
    boxplot(estimation_q_s)
    
    hist(estimation_q_m,breaks=10)
    boxplot(estimation_q_m)
    
    hist(estimation_q_l,breaks=10)
    boxplot(estimation_q_l)    
  }
  
  #Approximation de l'espérance par la moyenne empirique
  esperance_q_s <- mean(estimation_q_s)
  esperance_q_m <- mean(estimation_q_m)
  esperance_q_l <- mean(estimation_q_l)
  
  #Calcul du biais
  bias_s <- esperance_q_s - q
  bias_m <- esperance_q_m - q
  bias_l <- esperance_q_l - q
  
  #Calcul de la variance
  variance_s <- mean(estimation_q_s^2)-esperance_q_s^2
  variance_m <- mean(estimation_q_m^2)-esperance_q_m^2
  variance_l <- mean(estimation_q_l^2)-esperance_q_l^2
  
  #Calcul du MSE
  mse_s <- bias_s^2 + variance_s
  mse_m <- bias_m^2 + variance_m
  mse_l <- bias_l^2 + variance_l
  
  matrice <- matrix(c(bias_s,variance_s,mse_s,bias_m,variance_m,mse_m,bias_l,variance_l,mse_l), ncol = 3, byrow = TRUE, dimnames = list(c("q_s", "q_m", "q_l"), c("Biais", "Variance", "MSE")))
  
  return (matrice)
}

q <- q_p(0.95)

#1.3.a
n<-20
generate_sample()

#1.3.b et 1.3.c
result <- generate_estimation_q_p(TRUE)
result

#1.3.d
n_value <- c(20,50,100,250,400)
results <- array(0, dim = c(3, 3, length(n_value)))
for (i in 1:length(n_value)) {
  n<-n_value[i]
  results[,,i]=t(generate_estimation_q_p(FALSE))
}

#Graphiques

print(results[,,])

#Biais
plot(n_value,results[1,1,],xlab="n",ylab="Bias",main="Bias of q_s") #q_s
plot(n_value,results[1,2,],xlab="n",ylab="Bias",main="Bias of q_m") #q_m
plot(n_value,results[1,3,],xlab="n",ylab="Bias",main="Bias of q_l") #q_l

#Variance
plot(n_value,results[2,1,],xlab="n",ylab="Variance",main="Variance of q_s") #q_s
plot(n_value,results[2,2,],xlab="n",ylab="Variance",main="Variance of q_m") #q_m
plot(n_value,results[2,3,],xlab="n",ylab="Variance",main="Variance of q_l") #q_l

#MSE
plot(n_value,results[3,1,],xlab="n",ylab="MSE",main="MSE of q_s") #q_s
plot(n_value,results[3,2,],xlab="n",ylab="MSE",main="MSE of q_m") #q_m
plot(n_value,results[3,3,],xlab="n",ylab="MSE",main="MSE of q_l") #q_l
