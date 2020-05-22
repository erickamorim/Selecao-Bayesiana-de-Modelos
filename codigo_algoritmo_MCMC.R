rm(list=ls())
#instale os pacotes
install.packages("coda")
install.packages("mvtnorm")
install.packages("leaps")
#carregue os pacotes
library(coda)
library(mvtnorm)
library(leaps)
#carregue o arquivo funcoeimportantes.R
source("C: /.../funcoesimportantes.R")
#carregue os dados
dados=read.table("C:/.../dadosvinhos.txt",header = T)
n=nrow(dados)
k=ncol(dados)
names(dados)=c("Y",paste0("X", 1:(k-1)))
X=matrix( c(rep(1,n),as.matrix(dados[,2:k])) ,n,k)
Y=as.vector(dados[,1])

#estatisticas descritivas
summary(X)
summary(Y)
sd(Y);sd(X[,2]);sd(X[,3]);sd(X[,4]);sd(X[,5]);sd(X[,6])
cor(X[,-1])
win.graph()
plot(dados)
win.graph()
par(mfrow=c(3,2))
plot(X[,2],Y, ylab ="Qualidade" , xlab = "Claridade", pch=16,cex.axis=1.5,cex.lab=1.5)
plot(X[,3],Y, ylab ="Qualidade" , xlab = "Aroma", pch=16,cex.axis=1.5,cex.lab=1.5)
plot(X[,4],Y, ylab ="Qualidade" , xlab = "Corpo", pch=16,cex.axis=1.5,cex.lab=1.5)
plot(X[,5],Y, ylab ="Qualidade" , xlab = "Sabor", pch=16,cex.axis=1.5,cex.lab=1.5)
plot(X[,6],Y, ylab ="Qualidade" , xlab = "Afinação", pch=16,cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(3,2))
hist(Y, main="", xlab="Qualidade",cex.axis=1.5, cex.lab=1.5, freq = F)
hist(X[,2], main="", xlab="Claridade",cex.axis=1.5, cex.lab=1.5, freq = F)
hist(X[,3], main="", xlab="Aroma",cex.axis=1.5, cex.lab=1.5, freq = F)
hist(X[,4], main="", xlab="Corpo",cex.axis=1.5, cex.lab=1.5, freq = F)
hist(X[,5], main="", xlab="Sabor",cex.axis=1.5, cex.lab=1.5, freq = F)
hist(X[,6], main="", xlab="Afinação",cex.axis=1.5, cex.lab=1.5, freq = F)
coef(lm(Y~X1+X2+X3+X4+X5,data = dados))
# selecao de modelos
w=leaps(X[,2:k],Y,method = "adjr2")
modelo=cbind(w$which,w$adjr2)
modelo[order(-modelo[,6]),]
step(lm(Y~.,data = dados), scope = list(lower = ~1, upper = ~X1+X2+X3+X4+X5), direction = "both", test="F")
#----------------------------------------------------------------------------------------------------------------
#fazendo uma transformacao das variaveis explicativas mudando sua escala
X=apply(X,2,function(x){x/max(x)})
# Iniciio da aplicacao Bayesiana
#hiperparametros
a=2.1
b=1.1
vbetaj=50
gama1j=1
gama2j=1
#total de iteracoes, burn in, lag, valores iniciais e locais onde serao salvos os valores da cadeia
burnin=20000
napost=10000
lag=1
total=burnin+lag*napost
seqlag=seq(burnin+1,total,lag)
sigma2=array(1,c(total,1)) #local onde serao salvas todas as amostras
ssigma2=NULL #local onde serao salvas as amostras para analises
qj=array(0.5,c(total,k)) #local onde serao salvas todas as amostras
svaloresqj=NULL #local onde serao salvas as amostras para analises
qj.star=array(0,c(total,k)) #local onde serao salvas todas as amostras
svaloresqj.star=NULL #local onde serao salvas as amostras para analises
sbetas=array(0,c(total,k)) #local onde serao salvas todas as amostras
svalores=NULL #local onde serao salvas as amostras para analises
z=array(0,c(total,k))
svaloresz=NULL
xquad=somaXj2(X)
for (t in 2:total){
  D_0=D_menos_j(Y, X,t(t(sbetas[t-1,])))
  somaD_menos_0=sum(D_0[,1])
  vbeta0.star=1/( (n/sigma2[t-1]) + (1/vbetaj) )
  Mbeta0.star=vbeta0.star*(somaD_menos_0/sigma2[t-1])
  sbetas[t,1]=rnorm(1,Mbeta0.star,sqrt(vbeta0.star))
  for(j in 2:k){
    sbetas_aux=sbetas[t-1,];
    sbetas_aux[1:(j-1)]=sbetas[t,1:(j-1)]
    D_j=D_menos_j(Y,X,t(t(sbetas[t-1,])))
    somaXjD_menos_j=somaXjDmenosj(Y, X,t(t(sbetas_aux)))
    vbetaj.star=1/( (xquad[j]/sigma2[t-1]) + (1/vbetaj) )
    Mbetaj.star=vbetaj.star*(somaXjD_menos_j[j]/sigma2[t-1])
    razao = dnorm(0,Mbetaj.star, sqrt(vbetaj.star),log = T) - dnorm(0,0,sqrt(vbetaj), log = T)
    qj.star[t,j]=qj[t-1,j]/(qj[t-1,j]+(1-qj[t-1,j])*exp(razao) )
    u=runif(1,0,1)
    if(u<qj.star[t,j]){z[t,j]=1; sbetas[t,j]=rnorm(1, Mbetaj.star, sqrt(vbetaj.star))}
    qj[t,j]=rbeta(1, gama1j + z[t,j], gama2j - z[t,j] + 1)
  }
  D2=D2j(Y, X, t(t(sbetas[t,])))
  a.star= a + (n/2)
  b.star= b + 0.5*D2
  sigma2[t]=1/rgamma(1,shape = a.star, scale = 1/b.star)
  if (sum(seqlag==t)!=0) {
    svalores=rbind(svalores,sbetas[t,])
    ssigma2=rbind(ssigma2,sigma2[t])
    svaloresqj.star=rbind(svaloresqj,qj.star[t,])
    svaloresqj=rbind(svaloresqj,qj[t,])
    svaloresz=rbind(svaloresz,z[t,])
  }
  flush.console()
  cat(100*t/total, "% \r")
}
#----------------------------------------------------------------------------------------------------------------
win.graph()
par(mfrow=c(1,2))
plot(density(ssigma2),main="",ylab="Densidade",xlab=expression(sigma^2),bty="n",cex.axis=1.5,cex.lab=1.5)
plot(ssigma2, main="",ylab=expression(sigma^2),xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svaloresqj[,2]),main="",ylab="Densidade",xlab="q1",bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svaloresqj[,2], main="",ylab="q1",xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svaloresqj[,3]),main="",ylab="Densidade",xlab="q2",bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svaloresqj[,3], main="",ylab="q2",xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svaloresqj[,4]),main="",ylab="Densidade",xlab="q3",bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svaloresqj[,4], main="",ylab="q3",xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svaloresqj[,5]),main="",ylab="Densidade",xlab="q4",bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svaloresqj[,5], main="",ylab="q4",xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svaloresqj[,6]),main="",ylab="Densidade",xlab="q5",bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svaloresqj[,6], main="",ylab="q5",xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
#----------------------------------------------------------------------------------------------------------------
win.graph()
par(mfrow=c(1,2))
plot(density(svalores[,1]),main="",ylab="Densidade",xlab=expression(beta[0]),bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svalores[,1], main="",ylab=expression(beta[0]),xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svalores[,2]),main="",ylab="Densidade",xlab=expression(beta[1]),bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svalores[,2], main="",ylab=expression(beta[1]),xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svalores[,3]),main="",ylab="Densidade",xlab=expression(beta[2]),bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svalores[,3], main="",ylab=expression(beta[2]),xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svalores[,4]),main="",ylab="Densidade",xlab=expression(beta[3]),bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svalores[,4], main="",ylab=expression(beta[3]),xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svalores[,5]),main="",ylab="Densidade",xlab=expression(beta[4]),bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svalores[,5], main="",ylab=expression(beta[4]),xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svalores[,6]),main="",ylab="Densidade",xlab=expression(beta[5]),bty="n",cex.axis=1.5,cex.lab=1.5)
plot(svalores[,6], main="",ylab=expression(beta[5]),xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
#----------------------------------------------------------------------------------------------------------------
win.graph()
par(mfrow=c(1,2))
plot(density(svaloresqj.star[,2]),main="",ylab="Densidade",xlab="q1.star",bty="n",cex.axis=1.5,cex.lab=1.5)
abline(v=mean(svaloresqj.star[,2]),col="red",lwd=2)
abline(v=0.5,col="red",lwd=2,lty=3)
plot(svaloresqj.star[,2], main="",ylab="q1.star",xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svaloresqj.star[,3]),main="",ylab="Densidade",xlab="q2.star",bty="n",cex.axis=1.5,cex.lab=1.5)
abline(v=mean(svaloresqj.star[,3]),col="red",lwd=2)
abline(v=0.5,col="red",lwd=2,lty=3)
plot(svaloresqj.star[,3], main="",ylab="q2.star",xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svaloresqj.star[,4]),main="",ylab="Densidade",xlab="q3.star",bty="n",cex.axis=1.5,cex.lab=1.5)
abline(v=mean(svaloresqj.star[,4]),col="red",lwd=2)
abline(v=0.5,col="red",lwd=2,lty=3)
plot(svaloresqj.star[,4], main="",ylab="q3.star",xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svaloresqj.star[,5]),main="",ylab="Densidade",xlab="q4.star",bty="n",cex.axis=1.5,cex.lab=1.5)
abline(v=mean(svaloresqj.star[,5]),col="red",lwd=2)
abline(v=0.5,col="red",lwd=2,lty=3)
plot(svaloresqj.star[,5], main="",ylab="q4.star",xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
win.graph()
par(mfrow=c(1,2))
plot(density(svaloresqj.star[,6]),main="",ylab="Densidade",xlab="q5.star",bty="n",cex.axis=1.5,cex.lab=1.5)
abline(v=mean(svaloresqj.star[,6]),col="red",lwd=2)
abline(v=0.5,col="red",lwd=2,lty=3)
plot(svaloresqj.star[,6], main="",ylab="q5.star",xlab="Número de iterações",type="l",cex.axis=1.5,cex.lab=1.5)
#----------------------------------------------------------------------------------------------------------------
medias.qjstar=colMeans(svaloresqj.star)
aux=as.mcmc(svaloresqj.star)
HPD=NULL
for(i in 1:ncol(aux))
{
  HPDqjstar=HPDinterval(aux[,i],prob=0.95)
  HPD=rbind(HPD,HPDqjstar)
}
HPD
var(svaloresqj.star[,1])
var(svaloresqj.star[,2])
var(svaloresqj.star[,3])
var(svaloresqj.star[,4])
var(svaloresqj.star[,5])
#
mediasbeta0=sum(svalores[,1])/sum(svalores[,1]!=0)
mediasbeta1=sum(svalores[,2])/sum(svalores[,2]!=0)
mediasbeta2=sum(svalores[,3])/sum(svalores[,3]!=0)
mediasbeta3=sum(svalores[,4])/sum(svalores[,4]!=0)
mediasbeta4=sum(svalores[,5])/sum(svalores[,5]!=0)
mediasbeta5=sum(svalores[,6])/sum(svalores[,6]!=0)
#
varbeta0=(sum(svalores[,1]^2)- (sum(svalores[,1]!=0))*(mediasbeta0^2) )/( (sum(svalores[,1]!=0))-1)
varbeta1=(sum(svalores[,2]^2)- (sum(svalores[,2]!=0))*(mediasbeta1^2) )/( (sum(svalores[,2]!=0))-1)
varbeta2=(sum(svalores[,3]^2)- (sum(svalores[,3]!=0))*(mediasbeta2^2) )/( (sum(svalores[,3]!=0))-1)
varbeta3=(sum(svalores[,4]^2)- (sum(svalores[,4]!=0))*(mediasbeta3^2) )/( (sum(svalores[,4]!=0))-1)
varbeta4=(sum(svalores[,5]^2)- (sum(svalores[,5]!=0))*(mediasbeta4^2) )/( (sum(svalores[,5]!=0))-1)
varbeta5=(sum(svalores[,6]^2)- (sum(svalores[,6]!=0))*(mediasbeta5^2) )/( (sum(svalores[,6]!=0))-1)
#
varbeta0;varbeta1;varbeta2;varbeta3;varbeta4;varbeta5
HPDbeta0=c(mediasbeta0-1.96*sqrt(varbeta0), mediasbeta0+1.96*sqrt(varbeta0) )
HPDbeta1=c(mediasbeta1-1.96*sqrt(varbeta1), mediasbeta1+1.96*sqrt(varbeta1) )
HPDbeta2=c(mediasbeta2-1.96*sqrt(varbeta2), mediasbeta2+1.96*sqrt(varbeta2) )
HPDbeta3=c(mediasbeta3-1.96*sqrt(varbeta3), mediasbeta3+1.96*sqrt(varbeta3) )
HPDbeta4=c(mediasbeta4-1.96*sqrt(varbeta4), mediasbeta4+1.96*sqrt(varbeta4) )
HPDbeta5=c(mediasbeta5-1.96*sqrt(varbeta5), mediasbeta5+1.96*sqrt(varbeta5) )
#estimativa dos betas originais
mediasbeta0;mediasbeta1/max(dados[,2]);mediasbeta2/max(dados[,3]);mediasbeta3/max(dados[,4])
mediasbeta4/max(dados[,5]);mediasbeta5/max(dados[,6])
c((mediasbeta1/max(dados[,2]))-1.96*sqrt(varbeta1/((max(dados[,2]))^2)),(mediasbeta1/max(dados[,2]))+1.96*sqrt(varbeta1/((max(dados[,2]))^2)));
c((mediasbeta2/max(dados[,3]))-1.96*sqrt(varbeta2/((max(dados[,3]))^2)),(mediasbeta2/max(dados[,3]))+1.96*sqrt(varbeta2/((max(dados[,3]))^2)));
c((mediasbeta3/max(dados[,4]))-1.96*sqrt(varbeta3/((max(dados[,4]))^2)),(mediasbeta3/max(dados[,4]))+1.96*sqrt(varbeta3/((max(dados[,4]))^2)));
c((mediasbeta4/max(dados[,5]))-1.96*sqrt(varbeta4/((max(dados[,5]))^2)),(mediasbeta4/max(dados[,5]))+1.96*sqrt(varbeta4/((max(dados[,5]))^2)));
c((mediasbeta5/max(dados[,6]))-1.96*sqrt(varbeta5/((max(dados[,6]))^2)),(mediasbeta5/max(dados[,6]))+1.96*sqrt(varbeta5/((max(dados[,6]))^2)))
#estimativas para sigma2
medias.sigma2=colMeans(ssigma2)
aux=as.mcmc(ssigma2)
HPDs2=NULL
for(i in 1:ncol(aux))
{
  HPDsigma2=HPDinterval(aux[,i],prob=0.95)
  HPDs2=rbind(HPDs2,HPDsigma2)
}
HPDs2
var(ssigma2)
#----------------------------------------------------------------------------------------------------------------
p1=sum(svaloresz[,2])/length(svaloresz[,2])
p2=sum(svaloresz[,3])/length(svaloresz[,3])
p3=sum(svaloresz[,4])/length(svaloresz[,4])
p4=sum(svaloresz[,5])/length(svaloresz[,5])
p5=sum(svaloresz[,6])/length(svaloresz[,6])