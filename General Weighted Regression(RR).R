library(MCPMod)
library(DoseFinding)
library(pryr)
library(investr)
library(Hmisc)
library(ggplot2)
library(numDeriv)
library(gridExtra)
library(ggpubr)
library(MASS)
library(parallel)
library(doSNOW)
inverse = function (y,param,model, lower = 0, upper = 5) {
  e<-uniroot(function (x) do.call(model,c(dose=x,as.list(param))) - y, lower = lower, upper = upper)[1]$root
  e<-ifelse(e>1,1,e)
  return(e)
}
optFun<-function(h,SUM,B)
{
  return(norm(as.matrix(SUM%*%h+B)))
  
}

f_inverse<- function (x,coefficients,model) {
  if(model=="linlog")
  {
    coefficients[3]=0.2;
    names(coefficients)[3]<-"off"
  }
  y=switch(model,"emax"= coefficients[3]*(x-coefficients[1])/(coefficients[2]-(x-coefficients[1])),"sigEmax"=coefficients[3]*((x-coefficients[1])/(coefficients[2]-(x-coefficients[1])))^(1/coefficients[4]),"linlog"={exp((x-coefficients[1])/coefficients[2])-coefficients[3]},"linear"=(x-coefficients[1])/coefficients[2]) 
  if(is.na(y))
    return(1)
  else 
    return(ifelse(y>1,1,ifelse(y<0,0.001,y)))
}
linlog<-function (dose, e0, delta, off = 0.2) 
{
  linear(log(dose + off), e0, delta)
}

z1<-function(x,model,coefficients,Delta)
{
  (1-min(abs((f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)-x)/f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)),0.9999))^2
}

z2<-function(x,model,coefficients,Delta)
{
  (1-min(abs((f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)-x)/f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)),0.9999)^2)^2
}

z3<-function(x,model,coefficients,Delta,dose=c(0,0.05,0.2,0.6,1))
{
  
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-sort(abs(dose-d_MED))
  h<-dist[dist!=0][1]
  (1-min(abs((d_MED-x)/h),0.9999)^2)
}

z4<-function(x,model,coefficients,Delta,dose=c(0,0.05,0.2,0.6,1))
{
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-sort(abs(dose-d_MED))
  h<-dist[dist!=0][1]
  (1-min(abs((d_MED-x)/h),0.9999)^2)^2
}

z5<-function(x,model,coefficients,Delta,dose=c(0,0.05,0.2,0.6,1))
{
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-sort(abs(dose-d_MED))
  h<-dist[dist!=0][2]
  (1-min(abs((d_MED-x)/h),0.9999)^2)
}

z6<-function(x,model,coefficients,Delta,dose=c(0,0.05,0.2,0.6,1))
{
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-sort(abs(dose-d_MED))
  h<-dist[dist!=0][2]
  (1-min(abs((d_MED-x)/h),0.9999)^2)^2
}

z7<-function(x,model,coefficients,Delta,dose=c(0,0.05,0.2,0.6,1))
{
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-abs(dose-d_MED)
  if(abs(x-d_MED)==min(dist))
    return(2)
  else 
    return(1)
}
model="emax"
Delta=0.4

phi1<-function(dose,resp,coefficients,model,Delta,W=NULL,doseR)
{ 
  weights <- c("weights0","weights1","weights2","weights3","weights4","weights5","weights6","weights7","weights8")
  weightsfn <- match(W, weights)
  w=switch(weightsfn,"1"=z0(dose,model,coefficients,Delta),
           "2"=z1(dose,model,coefficients,Delta),
           "3"=z2(dose,model,coefficients,Delta),
           "4"=z3(dose,model,coefficients,Delta), 
           "5"=z4(dose,model,coefficients,Delta,doseR),
           "6"=z5(dose,model,coefficients,Delta,doseR),
           "7"=z6(dose,model,coefficients,Delta,doseR),
           "8"=z7(dose,model,coefficients,Delta,doseR),
           "9"=z8(dose,model,coefficients,Delta,doseR)) 
  return(w*(resp-do.call(model,c(dose=dose,as.list(coefficients))))*numDeriv::grad(function(coefficients)do.call(model,c(dose=dose,as.list(coefficients))),coefficients))
  
}

irls.delta <- function(old, new) sqrt(sum((old - new)^2, na.rm = TRUE)/
                                        max(1e-5, sum(old^2, na.rm = TRUE)))

emaxGrad2<-function(dose,eMax,ed50)
{
  rbind(c(0,0,0),c(0,0,-dose/(dose + ed50)^2),c(0,-dose/(dose + ed50)^2,2*eMax*(dose/(dose + ed50)^3)))
  
}
sigEmaxGrad2<-function(dose,eMax,ed50,h)
{
  lg2 <- function(x) ifelse(x == 0, 0, log(x))
  den <- (ed50^h + dose^h)
  f1prime2<-jacobian(function(ed50){-ed50^(h - 1) * dose^h * h /(ed50^h + dose^h)^2},ed50)
  f12prime<-numDeriv::grad(function(h)-ed50^(h - 1) * dose^h * h /(ed50^h + dose^h)^2,h)
  f21prime<-numDeriv::grad(function(ed50)dose^h * ed50^h * lg2(dose/ed50)/(ed50^h + dose^h)^2,ed50)
  f2prime2<-numDeriv::grad(function(h)dose^h * ed50^h * lg2(dose/ed50)/(ed50^h + dose^h)^2,h)
  
  rbind(c(0,0,0,0),c(0,0,-ed50^(h - 1) * dose^h * h /den^2,dose^h * ed50^h * lg2(dose/ed50)/den^2),
        c(0,-ed50^(h - 1) * dose^h * h /den^2,eMax*f1prime2,eMax* f12prime),
        c(0,dose^h * ed50^h * lg2(dose/ed50)/den^2,eMax*f21prime,eMax* f2prime2))
  
}
linearGrad2<-function(dose)
  matrix(0,2,2)
linlogGrad2<-function(dose)
  matrix(0,2,2)

Weighted_MED<-function(data,model,Delta,tol=1e-03,weights=c("weights0","weights1","weights2","weights3","weights4","weights5","weights6","weights7","weights8"))
{  
  dataFit <- data.frame(dose = sort(unique(data$dose)), resp = as.numeric(tapply(data$resp,data$dose,mean)))
  dose<-dataFit$dose
  resp<-dataFit$resp
  
  beta<-list()
  w=vector()
  converged=FALSE
  FLAG=FALSE
  coefficients<-fitMod(dose,resp,data=dataFit,model=model,bnds=defBnds(1)[[model]],addArgs = list(off=0.2,scal=1))$coefs
  MED_orig<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  beta[[1]]=coefficients
  r<-list()
  r[[1]]<-ifelse(is.na(MED_orig),dose[1],MED_orig)
  convi<-0
  DelG<-DelG2<-DelW<-delw<-A1<-K1<-K2<-B1<-b<-list()
  
  
  for(itr in 1:100)
  { B=0;
  SUM=matrix(0, length(coefficients), length(coefficients))

  for(j in 1:length(dose))
  {
    #Algorithm 1
    DelG[[j]]=switch(model,"emax"=DoseFinding::emaxGrad(dose[j],coefficients[2],coefficients[3]),
                     "linear"=DoseFinding::linearGrad(dose[j]),
                     "sigEmax"=DoseFinding::sigEmaxGrad(dose[j],coefficients[2],coefficients[3],coefficients[4]),
                     "linlog"=DoseFinding::linlogGrad(dose[j],off=0.2))
    DelG2[[j]]=switch(model,"emax"=emaxGrad2(dose[j],coefficients[2],coefficients[3]),
                      "linear"=linearGrad2(dose[j]),
                      "sigEmax"=sigEmaxGrad2(dose[j],coefficients[2],coefficients[3],coefficients[4]),
                      "linlog"=linlogGrad2(dose[j]))
    delw[[j]]=as.vector(jacobian(function(coefficients)switch(weights,"weights0"=z0(dose[j],model,coefficients,Delta),
                                                              "weights1"=z1(dose[j],model,coefficients,Delta),
                                                              "weights2"=z2(dose[j],model,coefficients,Delta),
                                                              "weights3"=z3(dose[j],model,coefficients,Delta), 
                                                              "weights4"=z4(dose[j],model,coefficients,Delta),
                                                              "weights5"=z5(dose[j],model,coefficients,Delta,dose),
                                                              "weights6"=z6(dose[j],model,coefficients,Delta,dose),
                                                              "weights7"=z7(dose[j],model,coefficients,Delta,dose),
                                                              "weights8"=z8(dose[j],model,coefficients,Delta,dose)),coefficients))
    w[j]=switch(weights,"weights0"=z0(dose[j],model,coefficients,Delta),
                "weights1"=z1(dose[j],model,coefficients,Delta),
                "weights2"=z2(dose[j],model,coefficients,Delta),
                "weights3"=z3(dose[j],model,coefficients,Delta), 
                "weights4"=z4(dose[j],model,coefficients,Delta),
                "weights5"=z5(dose[j],model,coefficients,Delta,dose),
                "weights6"=z6(dose[j],model,coefficients,Delta,dose),
                "weights7"=z7(dose[j],model,coefficients,Delta,dose),
                "weights8"=z8(dose[j],model,coefficients,Delta,dose))  
    if(any(delw[[j]][delw[[j]]!=0]>=1e-03))
    {
     # DelW[[j]]=diag(delw[[j]])
      A1[[j]]=-w[j]*t(DelG[[j]])%*% (DelG[[j]])
     # K1[[j]]=t(DelG[[j]])%*%rep(1,length(coefficients))%*%DelW[[j]]
      K1[[j]]=t(DelG[[j]])%*%delw[[j]]
      K2[[j]]=w[j]*DelG2[[j]]
      
      B1[[j]]=A1[[j]]+(resp[j]-do.call(model,c(dose=dose[j],as.list(coefficients))))*(K1[[j]]+K2[[j]])
    }
    else 
      B1[[j]]=-w[j]*(t(DelG[[j]])%*% (DelG[[j]]))+(resp[j]-do.call(model,c(dose=dose[j],as.list(coefficients))))*(w[j]*DelG2[[j]])

    SUM=SUM+B1[[j]]
    b[[j]]<-phi1(dose[j],resp[j],coefficients,model,Delta,W=weights,doseR=dose) 
    B=B+b[[j]]
    
  }
  
  h1<-qr.coef(qr(SUM,tol=1e-10),-B)
  h<-try(nlminb(start=h1, optFun, SUM= SUM, B=B, lower = rep(-1,length(coefficients)), upper = rep(1,length(coefficients))))
  
  if(inherits(h, "try-error"))
  {
    beta[[itr+1]]=beta[[itr]]
    coefficients=beta[[itr+1]]
    break;
  } 
  else {
    beta[[itr+1]]=beta[[itr]]+h$par
    coefficients=beta[[itr+1]]
  }
  r[[itr+1]]<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  convi <- irls.delta(r[[itr+1]], r[[itr]])
  converged<-convi<=tol
  if( any(abs(beta[[itr+1]])==Inf)||any(is.na(beta[[itr+1]]))||any(na.omit(coefficients[-1])<0))
  {
    coefficients=beta[[1]]
    FLAG=TRUE
    break;
  }
  else if(converged)
  { 
    coefficients=beta[[itr]];
    break;
  }
  
  }
  
  
  MED_new=ifelse(converged==TRUE||FLAG==TRUE,f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model),MED_orig)
  
  return(list(coefs=coefficients,MED=MED_new,converged=converged))
  
}

Weighted_MED2<-function(data,model,Delta,tol=1e-03,weights=c("weights0","weights1","weights2","weights3","weights4","weights5","weights6","weights7","weights8"))
{  
  #dataFit <- data.frame(dose = sort(unique(data$dose)), resp = as.numeric(tapply(data$resp,data$dose,mean)))
  #dose<-dataFit$dose
  #resp<-dataFit$resp
  data<-data[order(data$dose),]
  dose<-data$dose
  resp<-data$resp
  beta<-list()
  w=vector()
  converged=FALSE
  FLAG=FALSE
  coefficients<-fitMod(dose,resp,data=dataFit,model=model,bnds=defBnds(1)[[model]],addArgs = list(off=0.2,scal=1))$coefs
  MED_orig<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  beta[[1]]=coefficients
  r<-list()
  r[[1]]<-ifelse(is.na(MED_orig),dose[1],MED_orig)
  convi<-0
  DelG<-DelG2<-DelW<-delw<-A1<-K1<-K2<-B1<-b<-list()
  
  
  for(itr in 1:100)
  { B=0;
  SUM=matrix(0, length(coefficients), length(coefficients))
  
  for(j in 1:length(dose))
  {
    #Algorithm 1
    DelG[[j]]=switch(model,"emax"=DoseFinding::emaxGrad(dose[j],coefficients[2],coefficients[3]),
                     "linear"=DoseFinding::linearGrad(dose[j]),
                     "sigEmax"=DoseFinding::sigEmaxGrad(dose[j],coefficients[2],coefficients[3],coefficients[4]),
                     "linlog"=DoseFinding::linlogGrad(dose[j],off=0.2))
    DelG2[[j]]=switch(model,"emax"=emaxGrad2(dose[j],coefficients[2],coefficients[3]),
                      "linear"=linearGrad2(dose[j]),
                      "sigEmax"=sigEmaxGrad2(dose[j],coefficients[2],coefficients[3],coefficients[4]),
                      "linlog"=linlogGrad2(dose[j]))
    delw[[j]]=as.vector(jacobian(function(coefficients)switch(weights,"weights0"=z0(dose[j],model,coefficients,Delta),
                                                              "weights1"=z1(dose[j],model,coefficients,Delta),
                                                              "weights2"=z2(dose[j],model,coefficients,Delta),
                                                              "weights3"=z3(dose[j],model,coefficients,Delta), 
                                                              "weights4"=z4(dose[j],model,coefficients,Delta),
                                                              "weights5"=z5(dose[j],model,coefficients,Delta,dose),
                                                              "weights6"=z6(dose[j],model,coefficients,Delta,dose),
                                                              "weights7"=z7(dose[j],model,coefficients,Delta,dose),
                                                              "weights8"=z8(dose[j],model,coefficients,Delta,dose)),coefficients))
    w[j]=switch(weights,"weights0"=z0(dose[j],model,coefficients,Delta),
                "weights1"=z1(dose[j],model,coefficients,Delta),
                "weights2"=z2(dose[j],model,coefficients,Delta),
                "weights3"=z3(dose[j],model,coefficients,Delta), 
                "weights4"=z4(dose[j],model,coefficients,Delta),
                "weights5"=z5(dose[j],model,coefficients,Delta,dose),
                "weights6"=z6(dose[j],model,coefficients,Delta,dose),
                "weights7"=z7(dose[j],model,coefficients,Delta,dose),
                "weights8"=z8(dose[j],model,coefficients,Delta,dose))  
    if(any(delw[[j]][delw[[j]]!=0]>=1e-03))
    {
      # DelW[[j]]=diag(delw[[j]])
      A1[[j]]=-w[j]*t(DelG[[j]])%*% (DelG[[j]])
      # K1[[j]]=t(DelG[[j]])%*%rep(1,length(coefficients))%*%DelW[[j]]
      K1[[j]]=t(DelG[[j]])%*%delw[[j]]
      K2[[j]]=w[j]*DelG2[[j]]
      
      B1[[j]]=A1[[j]]+(resp[j]-do.call(model,c(dose=dose[j],as.list(coefficients))))*(K1[[j]]+K2[[j]])
    }
    else 
      B1[[j]]=-w[j]*(t(DelG[[j]])%*% (DelG[[j]]))
    SUM=SUM+B1[[j]]
    b[[j]]<-phi1(dose[j],resp[j],coefficients,model,Delta,W=weights,doseR=dose) 
    B=B+b[[j]]
    
  }
  
  h1<-qr.coef(qr(SUM),-B)
  h<-try(nlminb(start=h1, optFun, SUM= SUM, B=B,scale = 1, lower = rep(0,length(coefficients)), upper = rep(5,length(coefficients))))
  if(inherits(h, "try-error"))
  {
    beta[[itr+1]]=beta[[itr]]
    coefficients=beta[[itr+1]]
    break;
  } 
  else {
    beta[[itr+1]]=beta[[itr]]+h$par
    coefficients=beta[[itr+1]]
  }
  r[[itr+1]]<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  convi <- irls.delta(r[[itr+1]], r[[itr]])
  converged<-convi<=tol
  if( any(abs(beta[[itr+1]])==Inf)||any(is.na(beta[[itr+1]]))||any(na.omit(coefficients[-1])<0))
  {
    coefficients=beta[[1]]
    FLAG=TRUE
    break;
  }
  else if(converged)
  { 
    coefficients=beta[[itr]];
    break;
  }
  
  }
  
  
  MED_new=ifelse(converged==TRUE||FLAG==TRUE,f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model),MED_orig)
  
  return(list(coefs=coefficients,MED=MED_new,converged=converged))
  
}
var_weighted<-function(coefficients,dose,model,data,S,weights,Delta,ssize)
{
  DelG<-matrix(NA,length(data$dose),length(coefficients))
  for(j in 1:length(data$dose))
    DelG[j,]=switch(model,"emax"=DoseFinding::emaxGrad(data$dose[j],coefficients[2],coefficients[3]),
                    "linear"=DoseFinding::linearGrad(data$dose[j]),
                    "sigEmax"=DoseFinding::sigEmaxGrad(data$dose[j],coefficients[2],coefficients[3],coefficients[4]),
                    "linlog"=DoseFinding::linlogGrad(data$dose[j],off=0.2))
  
  w<-switch(weights,"weights0"="sapply(data$dose,z0,model,coefficients,Delta)"
            ,"weights1"="sapply(data$dose,z1,model,coefficients,Delta)"
            ,"weights2"="sapply(data$dose,z2,model,coefficients,Delta)"
            ,"weights3"="sapply(data$dose,z3,model,coefficients,Delta)"
            ,"weights4"="sapply(data$dose,z4,model,coefficients,Delta)"
            ,"weights5"="sapply(data$dose,z5,model,coefficients,Delta,dose)"
            ,"weights6"="sapply(data$dose,z6,model,coefficients,Delta,dose)"
            ,"weights7"="sapply(data$dose,z7,model,coefficients,Delta,dose)"
            ,"weights8"="sapply(data$dose,z8,model,coefficients,Delta,dose)")  
  
  #N<-diag(ssize)
  #N<-N/nrow(data)
  #N=diag(length(dose))
  #N<-N/tr(N)
 #N<-diag(length(data$resp))
 #N<-N/tr(N)
  N<-diag(length(data$dose))
  N<-N/tr(N)
  Wt<-eval(parse(text=w))
  
  Wt<-diag(Wt)
 
   Wt<-Wt/tr(Wt)
  V<-t(DelG)%*%N%*%Wt^2%*%S%*%DelG
  W<--t(DelG)%*%N%*%Wt%*%DelG
  Winv<-try(solve(W),silent=TRUE)
  if(class(Winv)=="try-error")
  Winv<- ginv(W)
  Var<- Winv%*%V%*%t(Winv)/sum(ssize)
  
  
  return(Var)
}
