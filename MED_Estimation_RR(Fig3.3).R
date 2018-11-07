source("C:/Users/saha/Desktop/Weighted Regression MED/Code for Thesis/General Weighted Regression(RR).R")

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

data<-list()
SSIZE=c(25,50,75,100)

seed=readRDS("C:/Users/saha/Desktop/Weighted Regression MED/Seed_Misp.RDS")
Delta=0.4
dose=c(0,0.05,0.2,0.6,1)
model="emax"
indCfit=0
indCfitp=0
MED_orig<- f_inverse(Delta+0.32,c( e0=0.32,eMax=0.74,ed50=0.14),model)
nocl<-detectCores()
cl<-makeCluster(nocl-1) 
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
clusterExport(cl,c("z1","z2","z3","z4","z5","z6","Weighted_MED","f_inverse","genDFdata","linlog","emax"))

#Segment 1: MED Estimation
nsim=5000
Delta=0.4
seed=sample(1:100000000,nsim)
doses <- c(0,0.05,0.2,0.6,1)
fmodels <- Mods(linear = NULL,linlog=NULL, emax = 0.2,
                exponential =0.28,sigEmax=c(0.4,4),
                doses=doses, placEff = 0.2, maxEff = 0.6,
                addArgs=list(off=0.2))
#saveRDS(seed,"Seed_Misp.RDS")
## calculate doses giving an improvement of 0.3 over placebo
MED<-TD(fmodels, Delta=0.4)
model="emax"

system.time(EstAll<-foreach(s=1:4)%:%
              foreach(i=1:nsim,.combine="rbind")%dopar% { 
                set.seed(seed[i])
                data=MCPMod::genDFdata(model, fmodels[[model]], doses, SSIZE[s],0.65)
                fit<-DoseFinding::fitMod(dose,resp,data,model=model,bnds=defBnds(1)[[model]],addArgs = list(off=0.2))
                MED_est0<-f_inverse(fit$coefs[1]+Delta,fit$coefs,model=model)
                wfit1<-Weighted_MED(data,model=model,Delta=Delta,weights="weights1")
                MED_est1<-f_inverse(wfit1$coefs[1]+Delta,wfit1$coefs,model=model)
                wfit2<-Weighted_MED(data,model=model,Delta=Delta,weights="weights2")
                MED_est2<-f_inverse(wfit2$coefs[1]+Delta,wfit2$coefs,model=model)
                wfit3<-Weighted_MED(data,model=model,Delta=Delta,weights="weights3")
                MED_est3<-f_inverse(wfit3$coefs[1]+Delta,wfit3$coefs,model=model)
                wfit4<-Weighted_MED(data,model=model,Delta=Delta,weights="weights4")
                MED_est4<-f_inverse(wfit4$coefs[1]+Delta,wfit4$coefs,model=model)
                wfit5<-Weighted_MED(data,model=model,Delta=Delta,weights="weights5")
                MED_est5<-f_inverse(wfit5$coefs[1]+Delta,wfit5$coefs,model=model)
                wfit6<-Weighted_MED(data,model=model,Delta=Delta,weights="weights6")
                MED_est6<-f_inverse(wfit6$coefs[1]+Delta,wfit6$coefs,model=model)
                MEDest<-c(MED_est0,MED_est1,MED_est2,MED_est3,MED_est4,MED_est5,MED_est6)
                return(MEDest)
              })

for(s in 1:4)
{
WeightsInd<-"MED_Fitted"
weights=c("weights1","weights2","weights3","weights4","weights5","weights6")
for(i in 1:length(weights))
  WeightsInd<-c(WeightsInd,paste("MED_Weight",i,sep=""))
MED_Simulation=EstAll[[s]]
MED_Actual<-f_inverse(0.2+Delta, c( e0=0.2,eMax=0.7,ed50=0.2),model="emax")
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
}

#SigEmax Model
model="sigEmax"

system.time(EstAll_SigEmax<-foreach(s=1:4)%:%
              foreach(i=1:nsim,.combine="rbind")%dopar% { 
                set.seed(seed[i])
                data=MCPMod::genDFdata(model, fmodels[[model]], doses, SSIZE[s],0.65)
                fit<-DoseFinding::fitMod(dose,resp,data,model=model,bnds=defBnds(1)[[model]],addArgs = list(off=0.2))
                MED_est0<-f_inverse(fit$coefs[1]+Delta,fit$coefs,model=model)
                wfit1<-Weighted_MED(data,model=model,Delta=Delta,weights="weights1")
                MED_est1<-f_inverse(wfit1$coefs[1]+Delta,wfit1$coefs,model=model)
                wfit2<-Weighted_MED(data,model=model,Delta=Delta,weights="weights2")
                MED_est2<-f_inverse(wfit2$coefs[1]+Delta,wfit2$coefs,model=model)
                wfit3<-Weighted_MED(data,model=model,Delta=Delta,weights="weights3")
                MED_est3<-f_inverse(wfit3$coefs[1]+Delta,wfit3$coefs,model=model)
                wfit4<-Weighted_MED(data,model=model,Delta=Delta,weights="weights4")
                MED_est4<-f_inverse(wfit4$coefs[1]+Delta,wfit4$coefs,model=model)
                wfit5<-Weighted_MED(data,model=model,Delta=Delta,weights="weights5")
                MED_est5<-f_inverse(wfit5$coefs[1]+Delta,wfit5$coefs,model=model)
                wfit6<-Weighted_MED(data,model=model,Delta=Delta,weights="weights6")
                MED_est6<-f_inverse(wfit6$coefs[1]+Delta,wfit6$coefs,model=model)
                MEDest<-c(MED_est0,MED_est1,MED_est2,MED_est3,MED_est4,MED_est5,MED_est6)
                return(MEDest)
              })
for(s in 1:4)
{
  WeightsInd<-"MED_Fitted"
  weights=c("weights1","weights2","weights3","weights4","weights5","weights6")
  for(i in 1:length(weights))
    WeightsInd<-c(WeightsInd,paste("MED_Weight",i,sep=""))
  MED_Simulation=EstAll[[s]]
  MED_Actual<-f_inverse(0.2+Delta, c( e0=0.2,eMax=0.66,ed50=0.4,h=4),model=model)
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
}
