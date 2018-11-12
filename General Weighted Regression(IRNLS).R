library(MCPMod)
library(DoseFinding)
library(pryr)
library(investr)
library(Hmisc)
library(ggpubr)
library(ggplot2)
library(gridExtra)
#weights<-"(1-abs((f_inverse(coefficients[1]+0.4,coefficients,model)-dose)/f_inverse(coefficients[1]+0.4,coefficients,model))^2)^2"
inverse = function (y,param,model, lower = 0, upper = 5) {
  e<-uniroot(function (x) do.call(model,c(dose=x,as.list(param))) - y, lower = lower, upper = upper)[1]$root
  e<-ifelse(e>1,1,e)
  return(e)
}
irls.delta <- function(old, new) sqrt(sum((old - new)^2, na.rm = TRUE)/
                                        max(1e-5, sum(old^2, na.rm = TRUE)))

f_inverse<- function (x,coefficients,model) {
  if(model=="linlog")
  {
    coefficients[3]=0.2;
    names(coefficients)[3]<-"off"
  }
  y=switch(model,"emax"= coefficients[3]*(x-coefficients[1])/(coefficients[2]-(x-coefficients[1])),"sigEmax"=coefficients[3]*((x-coefficients[1])/(coefficients[2]-(x-coefficients[1])))^(1/coefficients[4]),"linlog"=coefficients[3]*{exp((x-linlog(0,coefficients[1],coefficients[2],coefficients[3]))/coefficients[2])-1},"linear"=(x-coefficients[2])/coefficients[1]) 
  if(is.na(y))
    return(1)
  else 
    return(ifelse(y>1,1,ifelse(y<0,0.001,y)))
}



z0<-function(x,model,coefficients,Delta)
{
  (1-min(abs((f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)-x)/f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)),0.9999)^2)
}
z1<-function(x,model,coefficients,Delta)
{
  (1-min(abs((f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)-x)/f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)),0.9999))^2
}
z2<-function(x,model,coefficients,Delta)
{
  (1-min(abs((f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)-x)/f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)),0.9999))
}

z3<-function(x,model,coefficients,Delta)
{
  (1-min(abs((f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)-x)/f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)),0.9999)^2)^2
}

z4<-function(x,model,coefficients,Delta)
{
  (1-min(abs((f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)-x)/f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)),0.9999)^3)^3
}
z5<-function(x,model,coefficients,Delta,dose)
{
  dose<-sort(unique(dose))
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-sort(abs(dose-d_MED))
  h<-dist[dist!=0][1]
  (1-min(abs((d_MED-x)/h),0.9999)^2)^2
}
z6<-function(x,model,coefficients,Delta,dose)
{
  
  dose<-sort(unique(dose))
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-sort(abs(dose-d_MED))
  h<-dist[dist!=0][1]
  (1-min(abs((d_MED-x)/h),0.9999)^2)
}
z7<-function(x,model,coefficients,Delta,dose)
{
  
  dose<-sort(unique(dose))
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-sort(abs(dose-d_MED))
  h<-dist[dist!=0][2]
  (1-min(abs((d_MED-x)/h),0.9999)^2)^2
}
z8<-function(x,model,coefficients,Delta,dose)
{
  
  dose<-sort(unique(dose))
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-sort(abs(dose-d_MED))
  h<-dist[dist!=0][2]
  (1-min(abs((d_MED-x)/h),0.9999)^2)
}


z9<-function(x,model,coefficients,Delta,dose=c(0,0.05,0.2,0.6,1))
{
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-abs(dose-d_MED)
  if(abs(x-d_MED)==min(dist))
    return(2)
  else 
    return(1)
}


Weighted_fitMod<- function(data,model=model,bnds=NULL,wts,maxit = 100,tol = 0.001,off=0.2,convergence_criteria="MED",Delta=Delta)
{ 
  #define bounds
  if(is.null(bnds)){bnds=defBnds(1)[[model]]}
   t<-environment(fitMod)
  #data with doses and mean responses
  dataFit <- data.frame(dose=sort(unique(data$dose)),resp=as.numeric(tapply(data$resp,data$dose,mean)))
 
  A = matrix(NA, nrow=maxit, ncol=length(formalArgs(model))+1+length(dose),dimnames=list(c(1:maxit),c(formalArgs(model)[-1],"MED","resid",paste("weightes_dose",dose," "))), byrow = TRUE)
   w<-switch(wts,"weights0"="sapply(dataFit$dose,z0,model,coefficients,Delta)"
            ,"weights1"="sapply(dataFit$dose,z1,model,coefficients,Delta)"
            ,"weights2"="sapply(dataFit$dose,z2,model,coefficients,Delta)"
            ,"weights3"="sapply(dataFit$dose,z3,model,coefficients,Delta)"
            ,"weights4"="sapply(dataFit$dose,z4,model,coefficients,Delta)"
            ,"weights5"="sapply(dataFit$dose,z5,model,coefficients,Delta,dose)"
            ,"weights6"="sapply(dataFit$dose,z6,model,coefficients,Delta,dose)"
            ,"weights7"="sapply(dataFit$dose,z7,model,coefficients,Delta,dose)"
            ,"weights8"="sapply(dataFit$dose,z8,model,coefficients,Delta,dose)"
            ,"weights9"="sapply(dataFit$dose,z9,model,coefficients,Delta,dose)")  
  
  #bounds
  bnds <- switch(model,"emax"=matrix(bnds,nrow=1),"sigEmax"=matrix(bnds,nrow=2))

   if(model == "emax"|model == "exponential"){
      dim <- 1
      if(!is.matrix(bnds))
        bnds <- matrix(bnds, nrow = 1)
    } else {
      dim <- 2
    }
  #compute resXY
  dose <- dataFit$dose
  resp <- dataFit$resp
  converged=FALSE
 for(itr in seq_len(maxit))
 {
  #variances and group sizes
  if(itr==1)
  {
  n <- as.vector(table(data$dose))
  weights <- n/sum(n)
  }
  else
  {wt=eval(parse(text=w))
   weights=(sqrt(wt)/sum(sqrt(wt)))^2
  }

  form <- paste("resp~1")
  m <- model.matrix(as.formula(form),dataFit)
  clinS <- diag(sqrt(weights))
  qrX <- qr(clinS %*% m)
  resXY <- as.numeric(qr.resid(qrX,sqrt(weights)*resp))
  
  #make a grid with starting values for pow
  #taken from N<-ifelse(dim==1,ctrl$gridSize$dim1,ctrl$gridSize$dim2) in DoseFinding

  N<-ifelse(dim==1,30,144)
  nodes<-do.call(t[["getGrid"]],list(N,bnds,dim))
  
 
  #compute Zmat
  getPred <- switch(model,"emax"=function(vec,x)emax(x,0,1,vec),
                          "sigEmax"=function(vec,x)sigEmax(x,0,1,vec[1],vec[2]))
  Zmat <- apply(nodes, 1, getPred, x = dose)
  Zmat <- clinS%*%Zmat
  resZmat <- qr.resid(qrX,Zmat)
  
  #compute RSS
  colsms1 <- colSums(resZmat*resXY)
  colsms2 <- colSums(resZmat*resZmat)
  RSSvec <- sum(resXY*resXY)-(colsms1*colsms1)/colsms2
  indMin <- which.min(RSSvec)
  strt <- nodes[indMin,]
  resid <- RSSvec[indMin]
  
  #starting values
if(dim == 1){
  N<-30
  dif <- (bnds[2]-bnds[1])/N
  bnds[1] <- max(c(strt-1.1*dif),bnds[1])
  bnds[2] <- min(c(strt+1.1*dif),bnds[2])
  if(bnds[1]==bnds[2])bnds<-c(0,10)
}
  #local optimization
  optFunc <- function(par,x,qrX,resXY,clinS,model){
   Z <- switch(model,"emax"=emax(x,0,1,par),
                          "sigEmax"=sigEmax(x,0,1,par[1],par[2]))
  Z <- clinS%*%Z
    resXZ <- try(qr.resid(qrX,Z))
    if (inherits(resXZ,"try-error")) {return(NA)}
    sumrsXYrsXZ <- sum(resXY*resXZ)
    sum(resXY*resXY)-sumrsXYrsXZ*sumrsXYrsXZ/sum(resXZ*resXZ)
  }
   if(dim == 1){ # one-dimensional models
    optobj <-optimize(optFunc, c(bnds[1],bnds[2]), x=dose,qrX=qrX,resXY=resXY,model=model,clinS=clinS)
    coefs <- optobj$minimum
    RSS <- optobj$objective
  } else {
    optobj <- try(nlminb(strt, optFunc, x=dose, qrX=qrX, resXY=resXY,
                         model = model,
                         lower = bnds[,1], upper = bnds[,2], clinS=clinS))
    if(inherits(optobj, "try-error")){
      coefs <- RSS <- NA}
   else {
      coefs <- optobj$par
      RSS <- optobj$objective
    }
  }
    
  
  #calculation of linear coefficients
  nam1 <- switch(model, emax = c("eMax", "ed50"),linlog=c("delta"),
                     sigEmax = c("eMax", "ed50", "h"),
                     logistic = c("eMax", "ed50", "delta"),
                     exponential = c("e1", "delta"),
                     betaMod = c("eMax", "delta1", "delta2"))

  f0 <-switch(model,"emax"=emax(dose,0,1,coefs),
                          "sigEmax"=sigEmax(dose,0,1,coefs[1],coefs[2]))

  X <- cbind(1,f0)
  par0 <- qr.coef(qr(clinS%*%X),clinS%*%resp)
  par <- c(par0,coefs)
  names(par)<- c("e0", nam1)
  df <- sum(weights)-length(par)
  RSS <- sum((data$resp-sapply(data$dose,function(d)do.call(model,c(dose=d,as.list(par)))))^2)
  m<-list(coefs=par, resid = RSS)
  A[itr,1:length(par)]=m$coefs
  A[itr,length(par)+1]=as.numeric(ifelse(inherits(try(f_inverse(do.call(model,c(dose=0,as.list(m$coefs)))+Delta,m$coefs,model),silent=TRUE), "try-error"),1,f_inverse(do.call(model,c(dose=0,as.list(m$coefs)))+Delta,m$coefs,model)))
      A[itr,length(par)+2]=as.numeric(m$resid)
      A[itr,(length(par)+3):ncol(A)]=as.numeric(weights)
  
  if(itr!=1)
      {
        if(convergence_criteria=="MED")
        {
          placEff1<-do.call(model, c(min(dose), as.list(m$coefs)))
          placEff2<-do.call(model, c(min(dose), as.list(coefficients)))
         
          converged <-  irls.delta(A[itr-1,length(par)+1], A[itr,length(par)+1])<=tol

          if(do.call(f_inverse,c(Delta+placEff1,list(m$coefs),model=model))<0)
          {
            m$coefs<-A[itr-1,][1:length(par)]
            converged=TRUE
            break;
          }  
         }
      }
    coefficients=m$coefs
       if(converged)
        break

  }

m$converged<-converged
if(converged!=TRUE)
m$coefs<-fitMod(dose,resp,data,model,bnds=defBnds(1)[[model]])$coefs
m$itr<-itr
m$data<-data
dat<-list(A,m)
names(dat)<-list("A","m")
return(dat)
}

var_weighted_IRNLS<-function(coefficients,dose,model,data,S,weights,Delta,ssize,predtype)
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
            ,"weights8"="sapply(data$dose,z8,model,coefficients,Delta,dose)"  
            ,"weights9"="sapply(dataFit$dose,z9,model,coefficients,Delta,dose)")  
  #N<-diag(ssize)
  #N<-N/nrow(data)
  #N=diag(length(dose))
  #N<-N/tr(N)
  #N<-diag(length(data$resp))
  #N<-N/tr(N)
  N<-diag(length(data$dose))
  #N<-N/tr(N)
  Wt<-eval(parse(text=w))
  
  Wt<-diag(Wt)
  
  Wt<-Wt/tr(Wt)
  V<-t(DelG)%*%N%*%Wt%*%S%*%DelG
  Var<- V/sum(ssize)
  return(Var)
}
confintwIRLS_MED<-function(level,coefficients,model,dose,data,weights,S,Delta,Est,ssize)
{
  b<-jacobian(function(coefficients)f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
              ,coefficients)
  
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  Confint<-matrix(0,nrow=1,ncol=2)
  var3<-(b%*%var_weighted_IRNLS(coefficients,dose,model,data,S,weights,Delta,ssize)%*%t(b))
  
  Confint=Est+fac*sqrt(var3)
  return(Confint)
}


  