source("C:/Users/saha/Desktop/Weighted Regression MED/General Weighted Regression-IRNLS_SC3.R")
source("C:/Users/saha/Desktop/Weighted Regression MED/General Weighted Regression-Freiman_SC3.R")

library(numDeriv)
var_fitted<-function(coefficients,dose,model,data,S,ssize)
{
  DelG<-matrix(NA,length(dose),length(coefficients))
  for(j in 1:length(dose))
    DelG[j,]=switch(model,"emax"=DoseFinding::emaxGrad(dose[j],coefficients[2],coefficients[3]),
                    "linear"=DoseFinding::linearGrad(dose[j]),
                    "sigEmax"=DoseFinding::sigEmaxGrad(dose[j],coefficients[2],coefficients[3],coefficients[4]),
                    "linlog"=DoseFinding::linlogGrad(dose[j],off=0.2))
  
  N<-diag(rep(ssize,length(dose)))
  N<-N/tr(N)
  var2<-solve(t(DelG)%*%N%*%solve(S)%*%DelG)
  return(var2)
}



z6<-function(x,model,coefficients,Delta,dose=c(0,0.05,0.2,0.6,1))
{
  d_MED<-f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
  dist<-sort(abs(dose-d_MED))
  h<-dist[dist!=0][2]
  (1-min(abs((d_MED-x)/h),0.9999)^2)^2
}

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

confintwRR_MED<-function(level,coefficients,model,dose,data,weights,S,Delta,Est,ssize)
{
  
  b<-jacobian(function(coefficients)f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
              ,coefficients)
  
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  Confint<-matrix(0,nrow=1,ncol=2)
  var3<-(b%*%var_weighted(coefficients,dose,model,data,S,weights,Delta,ssize)%*%t(b))
  
  Confint=Est+fac*sqrt(var3)
  return(Confint)
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
confint_MED<-function(level,coefficients,model,fit,Delta)
{
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  b<-jacobian(function(coefficients)f_inverse(do.call(model,c(dose=0,as.list(coefficients)))+Delta,coefficients,model)
              ,coefficients)
  Confint<-matrix(0,nrow=1,ncol=2)
  var3<-b%*%vcov(fit)%*%t(b)
  Confint=TD(fit,Delta=Delta)+fac*sqrt(var3)
  return(Confint)
}


seed=readRDS("C:/Users/saha/Desktop/Weighted Regression MED/seedCovSigEmax.RDS")
Delta=0.4
dose=c(0,0.05,0.2,0.6,1)
model="emax"
nocl<-detectCores()
cl<-makeCluster(nocl) 
clusterEvalQ(cl,library(stats4))
clusterEvalQ(cl,library(DoseFinding))
clusterEvalQ(cl,library(MCPMod))
clusterEvalQ(cl,library(nloptr))
clusterEvalQ(cl,library(investr))
clusterEvalQ(cl,library(pryr))
clusterEvalQ(cl,library(Hmisc))
clusterEvalQ(cl,library(numDeriv))
clusterEvalQ(cl,library(psych))
clusterEvalQ(cl,library(parallel))
clusterEvalQ(cl,library(doSNOW))
clusterEvalQ(cl,library(MASS))
registerDoSNOW(cl)
SSIZE=c(25,50,75,100)
model="sigEmax"
nsim=2000
Delta=0.4
MED<- f_inverse(Delta+0.32,c( e0=0.32,eMax= 0.66,ed50=0.3,h=4),model=model)
system.time(CfitAll_sigEmax1<-foreach(s=1:4) %:%
              foreach(i=1:nsim)%dopar%
              { 
                set.seed(seed[i])
                data=MCPMod::genDFdata(model, c( e0=0.32,eMax= 0.66,ed50=0.3,h=4), c(0,0.05,0.2,0.6,1), SSIZE[s],0.65)
                # anovaMod <- lm(resp~factor(dose), data=data)
                #drFit <- coef(anovaMod)# placebo adjusted (=effect) estimates at doses
                # S <- vcov(anovaMod) # estimated covariance
                # gfit <- fitMod(dose, drFit, S=S, model = model,type = "general", bnds =defBnds(1)[[model]])
                fit <- fitMod(dose, resp,data,model = model, bnds = defBnds(1)[[model]])
                wfit<-Weighted_MED(data,model=model,Delta=Delta,weights="weights8")
                wRSS<-crossprod(data$resp-sapply(data$dose,function(d)do.call(model,c(dose=d,as.list(wfit$coefs)))))/(length(dose)*SSIZE[s]-length(wfit$coefs))
                wfit2<-Weighted_fitMod(data,model=model,wts="weights8",maxit = 100,tol = 0.001,off=0.2,convergence_criteria="MED",Delta=Delta)
                wRSS2<-crossprod(data$resp-sapply(data$dose,function(d)do.call(model,c(dose=d,as.list(wfit2$m$coefs)))))/(length(dose)*SSIZE[s]-length(+wfit2$m$coefs))
                MED_est<-f_inverse(Delta+wfit2$m$coefs[1],wfit2$m$coefs,model=model)
                SI=diag(as.numeric(wRSS2),SSIZE[s]*length(dose))
                CfitwIRNLS<-confintwIRLS_MED(level,wfit2$m$coefs,model,dose,data,weights="weights8",S=SI,Delta,MED_est,ssize=rep(SSIZE[s],length(dose)))
                Cfit1<-confint_MED(level=level,fit$coefs,model=model,fit,Delta=Delta)
                
                # Cfit2<-confint_MED(level=level,gfit$coefs,model=model,gfit,Delta=Delta)
                SI=diag(as.numeric(wRSS),SSIZE[s]*length(dose))
                CfitwRR<-confintwRR_MED(level=level,wfit$coefs,model=model,dose,data,weights="weights8",S=SI,Delta=Delta,wfit$MED,ssize=rep(SSIZE[s],length(dose)))
                #emxBobj <- bFitMod(dose, drFit, S=S,model=model, type = "bootstrap",nSim=1000, bnds =defBnds(1)[[model]])
                #tds <- TD(emxBobj, Delta=Delta)
                # Cfitb<-quantile(tds, c(1-level,level), na.rm = TRUE)
                cbind(Cfit1,CfitwRR,CfitwIRNLS)
              })
system.time(CfitAll_sigEmax2<-foreach(s=1:4) %:%
              foreach(i=1:nsim)%dopar%
              { 
                set.seed(seed[i])
                data=MCPMod::genDFdata(model, c( e0=0.32,eMax= 0.66,ed50=0.3,h=4), c(0,0.25,0.5,0.75,1), SSIZE[s],0.65)
               # anovaMod <- lm(resp~factor(dose), data=data)
                #drFit <- coef(anovaMod)# placebo adjusted (=effect) estimates at doses
               # S <- vcov(anovaMod) # estimated covariance
               # gfit <- fitMod(dose, drFit, S=S, model = model,type = "general", bnds =defBnds(1)[[model]])
                fit <- fitMod(dose, resp,data,model = model, bnds = defBnds(1)[[model]])
                wfit<-Weighted_MED(data,model=model,Delta=Delta,weights="weights8")
                wRSS<-crossprod(data$resp-sapply(data$dose,function(d)do.call(model,c(dose=d,as.list(wfit$coefs)))))/(length(dose)*SSIZE[s]-length(wfit$coefs))
                wfit2<-Weighted_fitMod(data,model=model,wts="weights8",maxit = 100,tol = 0.001,off=0.2,convergence_criteria="MED",Delta=Delta)
                wRSS2<-crossprod(data$resp-sapply(data$dose,function(d)do.call(model,c(dose=d,as.list(wfit2$m$coefs)))))/(length(dose)*SSIZE[s]-length(+wfit2$m$coefs))
                MED_est<-f_inverse(Delta+wfit2$m$coefs[1],wfit2$m$coefs,model=model)
                SI=diag(as.numeric(wRSS2),SSIZE[s]*length(dose))
                CfitwIRNLS<-confintwIRLS_MED(level,wfit2$m$coefs,model,dose,data,weights="weights8",S=SI,Delta,MED_est,ssize=rep(SSIZE[s],length(dose)))
                Cfit1<-confint_MED(level=level,fit$coefs,model=model,fit,Delta=Delta)
                
               # Cfit2<-confint_MED(level=level,gfit$coefs,model=model,gfit,Delta=Delta)
                SI=diag(as.numeric(wRSS),SSIZE[s]*length(dose))
                CfitwRR<-confintwRR_MED(level=level,wfit$coefs,model=model,dose,data,weights="weights8",S=SI,Delta=Delta,wfit$MED,ssize=rep(SSIZE[s],length(dose)))
                #emxBobj <- bFitMod(dose, drFit, S=S,model=model, type = "bootstrap",nSim=1000, bnds =defBnds(1)[[model]])
                #tds <- TD(emxBobj, Delta=Delta)
               # Cfitb<-quantile(tds, c(1-level,level), na.rm = TRUE)
                cbind(Cfit1,CfitwRR,CfitwIRNLS)
              })

CfitAll2<-CfitAll_sigEmax2
c=rep(0,4)

Cov<-matrix(0,nrow=4,ncol=3)
colnames(Cov)<-c("Wald1","Weighted","IRNLS")
rownames(Cov)<-SSIZE
IntervalCont<-function(x,vec)
{
  if(any(is.na(vec))==TRUE) return(NA)
  if(vec[1]==vec[2])return(0)
  else
    if(vec[1]<=x & x<=vec[2])return(1)
  else return(0)
}


for(k in 1:length(SSIZE))
{
  dat<-sapply(1:nsim, function(j)apply(CfitAll2[[k]][[j]],2,IntervalCont,x=MED))
  # dat1<-dat[,complete.cases(t(dat))] 
  Cov[k,]=sapply(1:3,function(j)sum(dat[j,], na.rm=TRUE))/nsim
}

#Bootstrap
source("C:/Users/saha/Desktop/Jyoti_Run/Prof_Lik_Baynv2.R") 
clusterExport(cl,c("curveInt","doseInt","globsearch","profileLR","ll","bootconfint"))
rnge <- c(0,1)
doseGrid <- seq(rnge[1],rnge[2],diff(rnge)/100)
dim(doseGrid) <- c(1,length(doseGrid))
nsim=2000

system.time(Cfit<-foreach(s=1:4)%:%
              foreach(i=1:nsim, .combine='cbind') %dopar%
              {set.seed(seed[i])
                data=MCPMod::genDFdata(model, c( e0=0.32,eMax= 0.66,ed50=0.3,h=4), c(0,0.05,0.2,0.6,1), SSIZE[s],0.65)
                CI_effect_boot<-bootconfint(data,model=model,nsamples=1000, nsteps=100, curvetype="effect", alpha=0.05,bounds=NULL)
                ind1_boot<-which(abs(CI_effect_boot$CI[1,]-Delta)==min(abs(CI_effect_boot$CI[1,]-Delta),na.rm=T))
                ind2_boot<-which(abs(CI_effect_boot$CI[2,]-Delta)==min(abs(CI_effect_boot$CI[2,]-Delta),na.rm=T))
                return(doseGrid[1,c(ind2_boot,ind1_boot)])})
clusterExport(cl,c("curveInt","doseInt","globsearch","profileLR","ll","bootconfint"))
rnge <- c(0,1)
doseGrid <- seq(rnge[1],rnge[2],diff(rnge)/100)
dim(doseGrid) <- c(1,length(doseGrid))
nsim=2000

system.time(Cfit2<-foreach(s=1:4)%:%
              foreach(i=1:nsim, .combine='cbind') %dopar%
              {set.seed(seed[i])
                data=MCPMod::genDFdata(model, c( e0=0.32,eMax= 0.66,ed50=0.3,h=4), c(0,0.25,0.5,0.75,1), SSIZE[s],0.65)
                CI_effect_boot<-bootconfint(data,model=model,nsamples=1000, nsteps=100, curvetype="effect", alpha=0.05,bounds=NULL)
                ind1_boot<-which(abs(CI_effect_boot$CI[1,]-Delta)==min(abs(CI_effect_boot$CI[1,]-Delta),na.rm=T))
                ind2_boot<-which(abs(CI_effect_boot$CI[2,]-Delta)==min(abs(CI_effect_boot$CI[2,]-Delta),na.rm=T))
                return(doseGrid[1,c(ind2_boot,ind1_boot)])})
MED<- ifelse(model=="emax",f_inverse(Delta+0.32,c( e0=0.32,eMax=0.74,ed50=0.14),"emax"),f_inverse(Delta+0.32,c( e0=0.32,eMax=0.66,ed50=0.3,h=4),"sigEmax"))

cov_prob=sapply(1:4,function(i)sum(apply(Cfit2[[i]],2,IntervalCont,x=MED))/nsim)

#PLApproach
source("C:/Users/saha/Desktop/Jyoti_Run/Prof_Lik_Baynv3.R") 
nsim=2000
system.time(Cfit3<-foreach(s=1:4)%:%
              foreach(i=1:nsim, .combine='cbind') %dopar%
              {set.seed(seed[i])
                data=MCPMod::genDFdata(model, c( e0=0.32,eMax=0.66,ed50=0.3,h=4), c(0,0.05,0.2,0.6,1), SSIZE[s],0.65)
                CI_effect<-tryCatch(curveInt(data,curvetype="effect",plot=F,model=model)$intervals,error=function(e) NA)
                if(is.na(any(c(CI_effect))))
                  return(c(999,999))
                else{
                  ind1<-which(abs(CI_effect[2,]-Delta)==min(abs(CI_effect[2,]-Delta),na.rm=T))
                  ind2<-which(abs(CI_effect[3,]-Delta)==min(abs(CI_effect[3,]-Delta),na.rm=T))
                  rm(intervals)
                  return(CI_effect[1,c(ind2,ind1)])}})

MED<- ifelse(model=="emax",f_inverse(Delta+0.32,c( e0=0.32,eMax=0.74,ed50=0.14),"emax"),f_inverse(Delta+0.32,c( e0=0.32,eMax=0.66,ed50=0.3,h=4),"sigEmax"))

cov_prob=sapply(1:4,function(i)sum(apply(Cfit4[[i]],2,IntervalCont,x=MED))/nsim)

nsim=2000
system.time(Cfit5<-foreach(s=1:4)%:%
              foreach(i=1:nsim, .combine='cbind') %dopar%
              {set.seed(seed[i+1000])
                data=MCPMod::genDFdata(model, c( e0=0.32,eMax=0.66,ed50=0.3,h=4), c(0,0.25,0.5,0.75,1), SSIZE[s],0.65)
                CI_effect<-tryCatch(curveInt(data,curvetype="effect",plot=F,model=model)$intervals,error=function(e) NA)
                if(is.na(any(c(CI_effect))))
                  return(c(999,999))
                else{
                  ind1<-which(abs(CI_effect[2,]-Delta)==min(abs(CI_effect[2,]-Delta),na.rm=T))
                  ind2<-which(abs(CI_effect[3,]-Delta)==min(abs(CI_effect[3,]-Delta),na.rm=T))
                  rm(intervals)
                  return(CI_effect[1,c(min(ind2),min(ind1))])}})


cov_prob=sapply(1:4,function(i)sum(apply(Cfit5[[i]],2,IntervalCont,x=MED))/nsim)

