# Funcoes importantes utilizadas no mcmc
#
##Calcula soma do quadrado de Di= Y-b0-b1X1-...-bkXk
D2j=function(Y, X, betas){
  preditorlinear = X%*%betas
  r=(Y-preditorlinear)^2
  return(sum(r))
}
#Calcula D menos o j-esimo componente
D_menos_j = function(Y, X, betas){
  Dmenosj= matrix(NA,nrow=nrow(X),ncol=ncol(X))
  for(j in 1:ncol(X))
  {
    predlinear_menosj= (X[,-j])%*%(betas[-j])
    Dmenosj[,j]= Y - predlinear_menosj
  }
  return(Dmenosj)
}

## Calcula o produto XjD_menosj
somaXjDmenosj=function(Y, X, betas){
  aux=NULL
  rj= D_menos_j(Y, X, betas)
  for(j in 1:ncol(X))
  {
    aux[j]=sum(as.numeric(X[,j])*as.numeric(rj[,j]))
  }
  return(aux)
}

#soma do quadrado de Xij em i

somaXj2=function(X){
  aux= NULL
  for(j in 1:ncol(X))
  {
    aux[j]=sum(X[,j]^2)
  }
  return(as.numeric(aux))
}
