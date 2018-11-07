library(DoseFinding)
library(psych)
s=20

seed<-readRDS("C:/Users/saha/Desktop/Weighted Regression MED/Code for Thesis/SeedforFig3.1.RDS")
set.seed(seed)
data<-MCPMod::genDFdata("emax", c(e0 = 0.3216, eMax = 0.7463, ed50 = 0.1422), c(0,0.05,0.2,0.6,1), s,0.65)

dose=unique(data$dose)
anovaMod <- lm(resp~factor(dose), data=data)
drFit <- coef(anovaMod)[2:length(dose)] # placebo adjusted (=effect) estimates at doses
S <- vcov(anovaMod)[2:length(dose),2:length(dose)] # estimated covariance
doseFit<-dose[-1]
doseSeq<-seq(0,max(dose),length=50)
ssize=rep(s,length(dose))
fit <- fitMod(dose,resp,data,model = "emax",bnds = defBnds(1)[["emax"]])

gfit <- fitMod(doseFit, drFit, S=S, model = "emax",placAdj=TRUE,type = "general", bnds = c(0.01, 2))
#plot(gfit, CI = TRUE, plotData = "meansCI")
level=0.95
x=gfit
if(class(x)=="DRMod")
  obj=x
if(class(x)=="MCPMod")
  obj<-x$mods[[1]]
doseNam<- attr(obj,"doseRespNam")[1]
respNam<- attr(obj,"doseRespNam")[2]
type<- attr(obj,"type")
placAdj<- attr(obj,"placAdj")
pList<-as.list(data)
pList$dos<-data[[doseNam]]
pList$mns<-data[[respNam]]
sdv<-sqrt(diag(data$S))
quant<-qnorm(1-(1-level)/2)
pList$lbndm<-pList$mns-quant*sdv
pList$ubndm<-pList$mns+quant*sdv

predtype=ifelse(placAdj,"effect-curve","ls-means")
CI=TRUE

pred<-predict(x,predType=predtype,doseSeq=doseSeq,se.fit=CI)
if(CI)
{
  quant<-qt(1-(1-level)/2,df=x$df)
  lbnd<-pred$fit-quant*pred$se.fit
  ubnd<-pred$fit+quant*pred$se.fit
  pred<-pred$fit
}
plotdf<-data.frame(rep(doseSeq,3),c(pred,lbnd,ubnd),rep(c("pred","LB","UB"),each=length(doseSeq)),attr(x,"model"))
names(plotdf)<-c(doseNam,respNam,"group","model")
source("C:/Users/saha/Desktop/Weighted Regression MED/Code for Thesis/General Weighted Regression(RR).R")


Delta=0.2
wfit<-Weighted_MED(data,model="emax",Delta=Delta,weights="weights5")
pred2<-list()
#Snew=diag(vcov(anovaMod))
wRSS<-crossprod(data$resp-sapply(data$dose,function(d)do.call(model,c(dose=d,as.list(wfit$coefs)))))/(sum(ssize)-length(wfit$coefs))

#cholmat<-var_weighted(wfit$coefs,dose,model,data,S=diag(Snew),weights="weights3",Delta,ssize,predtype)
cholmat<-var_weighted(wfit$coefs,dose,model,data,S=diag(as.numeric(wRSS),sum(ssize)),weights="weights5",Delta,ssize)

pred2$fit<-sapply(doseSeq,function(d)do.call(model,c(dose=d,as.list(wfit$coefs))))-do.call(model,c(dose=0,as.list(wfit$coefs)))
gr<-emaxGrad(doseSeq,wfit$coef[2],wfit$coef[3])
#gr0<-emaxGrad(0,wfit$coef[2],wfit$coef[3])
#gr<-gr[,-1,drop=FALSE]
#gr0<-gr0[,-1]
#gr<-sweep(gr,2,gr0,"-")

diffMat<-cbind(-1,diag(length(doseSeq)-1))
pred2$se.fit<-c(0,sqrt(diag(diffMat%*%gr%*%cholmat%*%t(gr)%*%t(diffMat))))

#pred2$se.fit<-sqrt(rowSums((gr%*%t(chol(cholmat)))^2))/10
if(CI)
{
  quant<-qnorm(1-(1-level)/2)
  lbndnew<-pred2$fit-quant*pred2$se.fit
  ubndnew<-pred2$fit+quant*pred2$se.fit
  prednew<-pred2$fit
}

plotdf1<-data.frame(rep(doseSeq,3),c(prednew,lbndnew,ubndnew),rep(c("pred_W","LB_W","UB_W"),each=length(doseSeq)),attr(x,"model"))
names(plotdf1)<-c(doseNam,respNam,"group","model")
plotdf$Method<-rep("N_CI",nrow(plotdf))
plotdf1$Method<-rep("W_CI",nrow(plotdf1))
coefAct<-c(e0 = 0.3216, eMax = 0.7463, ed50 = 0.1422)
predAct<-sapply(doseSeq,function(d)do.call(model,c(dose=d,as.list(coefAct))))-do.call(model,c(dose=0,as.list(coefAct)))
plotdfAct<-data.frame(doseSeq,predAct,rep("predAct",each=length(doseSeq)),attr(x,"model"))
names(plotdfAct)<-c(doseNam,respNam,"group","model")
plotdfAct$Method<-rep("Actual",nrow(plotdfAct))


plotdf2<-rbind(plotdf,plotdf1,plotdfAct)
rng <-  range(plotdf2[[respNam]])
dff <- diff(rng)

 
ylim <- c(rng[1] - 0.05 * dff, rng[2] + 0.05 * dff)

## produce plot
lattice.theme <- trellis.par.get()
lwd <- lattice.theme$superpose.symbol$lwd
col <- lattice.theme$superpose.symbol$col
form <- as.formula(paste(respNam, "~", doseNam, "|model", sep=""))
trellis.par.set(superpose.line = list(col=c("black","green","red","blue","black","black","black") ,lty=c(3,3,1,1,2,2,1),
                superpose.symbol = list(cex = 1.3, pch = 20), 
                reference.line = list(col = "gray", lty ="dotted"))) 
MEDAct<-f_inverse(0.3216+Delta,coefAct,model)
print(
  xyplot(form, groups = plotdf2$group, data = plotdf2,auto.key = list(space = "top",points=F, lines = T,cex=.9,text = c("Effect_Curve_CI (Wald)", "Effect_Curve_CI (Weighted Rgn)","Effect_Curve (Wald)","Effect_Curve (Weighted Rgn)","MED=0.052 (Dose Axis)","Delta=0.2 (Response Axis)","Effect Curve (Actual)")),pList=pList,
         ylim = ylim,xlab="Dose",ylab="Response", panel = function(x, y, ..., pList){
           panel.grid(h = -1, v = -1, col = "lightgrey", lty = 2)
            panel.abline(v = MEDAct,h=Delta, lty = 2,lwd=1)
           quant <- qnorm(1 - (1 - level)/2)
           for(i in 1:length(pList$dos)){
             llines(rep(pList$dos[i], 2),
                    c(pList$lbndm[i], pList$ubndm[i]),
                    lty=1, col = 1, ...)
           }
          panel.xyplot(x, y, col=c(1,2,1,3,4,3,1),lty=c(3,1,3,2,1,2,1),pch=3,cex=2, type="l", ...)
          }))
plotdfnew<-plotdf2[plotdf2$group%in% c("pred","predW","predAct"),]


print(
  xyplot(form, groups = plotdf2$group, data = plotdf2,auto.key = list(space = "top",points=F, lines = T,cex=.9,text = c("Effect_Curve_CI (Wald)", "Effect_Curve_CI (Weighted Rgn)","Effect_Curve (Wald)","Effect_Curve (Weighted Rgn)","MED=0.052 (Dose Axis)","Delta=0.2 (Response Axis)","Effect Curve (Actual)")),pList=pList,
         ylim = ylim,xlab="Dose",ylab="Response", panel = function(x, y, ..., pList){
           panel.grid(h = -1, v = -1, col = "lightgrey", lty = 2)
           panel.abline(v = MEDAct,h=Delta, lty = 2,lwd=1)
           quant <- qnorm(1 - (1 - level)/2)
           for(i in 1:length(pList$dos)){
             llines(rep(pList$dos[i], 2),
                    c(pList$lbndm[i], pList$ubndm[i]),
                    lty=1, col = 1, ...)
           }
           panel.xyplot(x, y, col=c(1,2,1,3,4,3,1),lty=c(3,1,3,2,1,2,1),pch=3,cex=2, type="l", ...)
         }))


plotdfnew<-plotdf2[plotdf2$group%in% c("pred","pred_W","predAct"),]
plotdfnew$group<-factor(plotdfnew$group)

print(
  xyplot(form, groups = plotdfnew$group, data = plotdfnew,auto.key = list(points=F, lines =T,cex=0.9,col=c(2,4,1),pch=c(2,4,1),text = c("Effect_Curve (Classical Approach)","Effect_Curve (New Approach)","Effect Curve (True Model)")),pList=pList,
         ylim = ylim,xlab="Dose",ylab="Response", panel = function(x, y, ..., pList){
           panel.grid(h = -1, v = -1, col = "lightgrey", lty = 1)
           panel.abline(v = MEDAct,h=Delta, lty = 2,lwd=1)
           quant <- qnorm(1 - (1 - level)/2)
           for(i in 1:length(pList$dos)){
             llines(rep(pList$dos[i], 2),
                    c(pList$lbndm[i], pList$ubndm[i]),
                    lty=1, col = 1, ...)
           }
           panel.xyplot(x, y, col=c(2,4,1),lty=c(2,4,1),pch=3,cex=2, type="l", ...)
         }))