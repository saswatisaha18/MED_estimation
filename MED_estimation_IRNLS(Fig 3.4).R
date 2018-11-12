library(parallel)
library(doSNOW)
source("C:/Users/saha/Desktop/Weighted Regression MED/Code for thesis/General Weighted Regression(IRNLS).R")
weights=c("weights0","weights3","weights5","weights6","weights7","weights8","weights9")
Delta=0.4
dose=c(0,0.05,0.2,0.6,1)
model="emax"
indCfit=0
indCfitp=0
MED_orig<- f_inverse(Delta+0.2,c( e0=0.2,eMax=0.7,ed50=0.2),model)
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
clusterExport(cl,c("z0","z3","z5","z6","z7","z8","z9","Weighted_fitMod","f_inverse"))
#Segment 1: MED Estimation
s=1
SSIZE=25
nsim=5000
Delta=0.4
conrate<-rep(0,length(weights))
names(conrate)=weights
data<-list()

model="emax"
system.time(EstAll0<-foreach(i=1:nsim,.combine="rbind")%dopar%
{  #set.seed(seed[i])
  MEDWt<-function(w,data,Delta)            
  {wfit2<-Weighted_fitMod(data=data,model=model,wts=w,maxit = 100,tol = 0.001,convergence_criteria="MED",Delta=Delta)
  #if(wfit2$m$converged==TRUE)
  #conrate[match(w,weights)]=conrate[match(w,weights)]+1
  MED=f_inverse(Delta+wfit2$m$coefs[1],wfit2$m$coefs,model=model)
  return(as.numeric(MED))}
  data[[i]]=MCPMod::genDFdata(model, c( e0=0.2,eMax=0.7,ed50=0.2), c(0,0.05,0.2,0.6,1), SSIZE[s],0.65)
  MED_est=sapply(weights,function(w)MEDWt(w,data[[i]],Delta))    
  fit<-DoseFinding::fitMod(dose,resp,data[[i]],model=model,bnds=defBnds(1)[[model]])
  MED_est0<-f_inverse(Delta+fit$coefs[1],fit$coefs,model=model)
  c(fit=MED_est0,MED_est)
})
WeightsInd<-character()
for(i in 1:length(weights))
  WeightsInd<-c(WeightsInd,paste("MED_Weight",i,sep=""))
WeightsInd=c("MED_fitted",WeightsInd)
MED_Actual<-f_inverse(0.2+Delta, c( e0=0.2,eMax=0.7,ed50=0.2),model="emax")

MED_Simulation=EstAll0
colnames(MED_Simulation)=WeightsInd
rownames(MED_Simulation)=1:nsim
MED_Summary<-data.frame(MED=c(MED_Simulation),MED_Class=rep(WeightsInd,each=nsim))
p10 <- ggplot(MED_Summary, aes(x = MED_Class, y = MED,fill=MED_Class)) +geom_boxplot(colour = "black",position = position_dodge(1), width = 0.4,alpha=0.7) +scale_y_continuous(name = "MED",limits=c(0, 1)) +ggtitle("Boxplot of MED under different methods")+geom_hline(yintercept = MED_Actual)+scale_fill_discrete(breaks=levels(MED_Summary$MED_Class),labels=levels(MED_Summary$MED_Class))+theme(legend.text=element_text(size=7),legend.title = element_text(colour="blue", size=8, face="bold"),axis.text.x = element_text(angle=90, hjust=1,size=8))
M1<-M2<-M3<-vector()

for(i in 1:(length(weights)+1))
{M1[i]<-mean((MED_Summary[which(MED_Summary[,2]==levels(MED_Summary$MED_Class)[i]),][,1]-MED_Actual)/MED_Actual)
M2[i]<-median((MED_Summary[which(MED_Summary[,2]==levels(MED_Summary$MED_Class)[i]),][,1]-MED_Actual)/MED_Actual)
M3[i]<-IQR((MED_Summary[which(MED_Summary[,2]==levels(MED_Summary$MED_Class)[i]),][,1]-MED_Actual)/MED_Actual)
}
R_i<-cbind(M1,M2,M3)

R_i<-round(R_i,2)
colnames(R_i)<-c("Mean Rel Dev","Median Rel Dev","IQR Rel Dev")
rownames(R_i)<-levels(MED_Summary$MED_Class)
R_i<-as.data.frame(R_i)
R_i<-cbind(Method=rownames(R_i),R_i)
stable.p <- ggtexttable(R_i,rows=NULL, theme = ttheme(colnames.style = colnames_style(color = "black", size = 10,fill = "grey80", linewidth = 0.5, linecolor = "white", parse = FALSE),base_size = 10.5))
text.p <- ggparagraph(text = "Relative Devation MED(R_i)=(MED_est-MED_actual)/MED_actual", face = "italic", size = 12, color = "black")
r1<-ggpar(p10,orientation="horizontal", main = "Simulated From: emax",font.main = c(10,"bold.italic", "red"),font.x = c(8, "bold", "#2E9FDF"),font.y = c(8, "bold", "#E7B800"),font.xtickslab = 8,font.ytickslab=8)
grid.arrange(r1,stable.p,nrow=2,as.table=TRUE,heights=c(6,6))


model="sigEmax"
system.time(EstAll<-foreach(i=1:nsim,.combine="rbind")%dopar%
{ MEDWt<-function(w,data,Delta)            
  {wfit2<-Weighted_fitMod(data=data,model=model,wts=w,maxit = 100,tol = 0.001,convergence_criteria="MED",Delta=Delta)
  #if(wfit2$m$converged==TRUE)
  #conrate[match(w,weights)]=conrate[match(w,weights)]+1
  MED=f_inverse(Delta+wfit2$m$coefs[1],wfit2$m$coefs,model=model)
  return(as.numeric(MED))}
   data[[i]]=MCPMod::genDFdata(model, c( e0=0.2,eMax=0.615,ed50=0.4,h=4), c(0,0.05,0.2,0.6,1), SSIZE[s],0.65)
   MED_est=sapply(weights,function(w)MEDWt(w,data[[i]],Delta))    
   fit<-DoseFinding::fitMod(dose,resp,data[[i]],model=model,bnds=defBnds(1)[[model]])
   MED_est0<-f_inverse(Delta+fit$coefs[1],fit$coefs,model=model)
   c(fit=MED_est0,MED_est)
              })
WeightsInd<-character()
for(i in 1:length(weights))
  WeightsInd<-c(WeightsInd,paste("MED_Weight",i,sep=""))
WeightsInd=c("MED_fitted",WeightsInd)
MED_Actual<-f_inverse(0.2+Delta, c( e0=0.2,eMax=0.615,ed50=0.4,h=4),model="sigEmax")

  MED_Simulation=EstAll
  colnames(MED_Simulation)=WeightsInd
  rownames(MED_Simulation)=1:nsim
MED_Summary<-data.frame(MED=c(MED_Simulation),MED_Class=rep(WeightsInd,each=nsim))
p10 <- ggplot(MED_Summary, aes(x = MED_Class, y = MED,fill=MED_Class)) +geom_boxplot(colour = "black",position = position_dodge(1), width = 0.4,alpha=0.7) +scale_y_continuous(name = "MED",limits=c(0, 1)) +ggtitle("Boxplot of MED under different methods")+geom_hline(yintercept = MED_Actual)+scale_fill_discrete(breaks=levels(MED_Summary$MED_Class),labels=levels(MED_Summary$MED_Class))+theme(legend.text=element_text(size=7),legend.title = element_text(colour="blue", size=8, face="bold"),axis.text.x = element_text(angle=90, hjust=1,size=8))
M1<-M2<-M3<-vector()

for(i in 1:(length(weights)+1))
{M1[i]<-mean((MED_Summary[which(MED_Summary[,2]==levels(MED_Summary$MED_Class)[i]),][,1]-MED_Actual)/MED_Actual)
M2[i]<-median((MED_Summary[which(MED_Summary[,2]==levels(MED_Summary$MED_Class)[i]),][,1]-MED_Actual)/MED_Actual)
M3[i]<-IQR((MED_Summary[which(MED_Summary[,2]==levels(MED_Summary$MED_Class)[i]),][,1]-MED_Actual)/MED_Actual)
}
R_i<-cbind(M1,M2,M3)

R_i<-round(R_i,2)
colnames(R_i)<-c("Mean Rel Dev","Median Rel Dev","IQR Rel Dev")
rownames(R_i)<-levels(MED_Summary$MED_Class)
R_i<-as.data.frame(R_i)
R_i<-cbind(Method=rownames(R_i),R_i)
stable.p <- ggtexttable(R_i,rows=NULL, theme = ttheme(colnames.style = colnames_style(color = "black", size = 10,fill = "grey80", linewidth = 0.5, linecolor = "white", parse = FALSE),base_size = 10.5))
text.p <- ggparagraph(text = "Relative Devation MED(R_i)=(MED_est-MED_actual)/MED_actual", face = "italic", size = 12, color = "black")
r1<-ggpar(p10,orientation="horizontal", main = "Simulated From: sigEmax",font.main = c(10,"bold.italic", "red"),font.x = c(8, "bold", "#2E9FDF"),font.y = c(8, "bold", "#E7B800"),font.xtickslab = 8,font.ytickslab=8)
grid.arrange(r1,stable.p,nrow=2,as.table=TRUE,heights=c(6,6))



