source("C:/Users/saha/Desktop/Weighted Regression MED/Code for Thesis/General Weighted Regression(RR).R")
library(DoseFinding)
library(parallel)
library(doSNOW)

data<-list()
SSIZE=c(25,50)
seed=readRDS("C:/Users/saha/Desktop/Weighted Regression MED/Code for Thesis/Seed_Misp.RDS")
doses=c(0,0.05,0.2,0.6,1)
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
clusterEvalQ(cl,library(pbapply))
clusterEvalQ(cl,library(doSNOW))
clusterEvalQ(cl,library(MASS))
registerDoSNOW(cl)
clusterExport(cl,c("z1","z2","z3","z4","z5","z6","Weighted_MED","f_inverse","genDFdata","linlog"))

#Segment 1: MED Estimation
nsim=5000
Delta=0.4
fmodels <- Mods(linear = NULL,linlog=NULL, emax = 0.2,
                exponential =0.28,sigEmax=c(0.4,4),
                doses=doses, placEff = 0.2, maxEff = 0.6,
                addArgs=list(off=0.2))
## calculate doses giving an improvement of 0.3 over placebo
MED<-TD(fmodels, Delta)

model_orig="sigEmax"
model_fit1="emax"
#Replace by "weights4" and "weights5" to generate the MED with the corresponding weights (Fig A.1 and A.2)
system.time(EstAll<-foreach(s=1:2)%:%
              foreach(i=1:nsim)%dopar% { 
                set.seed(seed[i])
                data=MCPMod::genDFdata(model_orig, fmodels[[model_orig]], doses, SSIZE[s],0.65)
                fit<-DoseFinding::fitMod(dose,resp,data,model=model_orig,bnds=defBnds(1)[[model_orig]],addArgs = list(off=0.2))
                x0=do.call(model_orig,c(dose=0,as.list(fit$coefs)))
                MEDest0<-f_inverse(x0+Delta,fit$coefs,model=model_orig)
                wfit1<-Weighted_MED(data,model=model_orig,Delta=Delta,weights="weights5")
                x0=do.call(model_orig,c(dose=0,as.list(wfit1$coefs)))
                MEDest1<-f_inverse(x0+Delta,wfit1$coefs,model=model_orig)
                
                fit<-DoseFinding::fitMod(dose,resp,data,model=model_fit1,bnds=defBnds(1)[[model_fit1]],addArgs = list(off=0.2))
                x0=do.call(model_fit1,c(dose=0,as.list(fit$coefs)))
                MEDest2<-f_inverse(x0+Delta,fit$coefs,model=model_fit1)
                wfit1<-Weighted_MED(data,model=model_fit1,Delta=Delta,weights="weights5")
                x0=do.call(model_fit1,c(dose=0,as.list(wfit1$coefs)))
                MEDest3<-f_inverse(x0+Delta,wfit1$coefs,model=model_fit1)
                return(rbind(MEDest0,MEDest1,MEDest2,MEDest3))
              })

for(s in 1:2)
{
  MED_Simulation<-list()
  MED_Simulation = matrix(NA, nrow=nsim, ncol=4,dimnames=list(c(1:nsim),c("MED_Ft_Orig","MED_Wt_Orig","MED_Ft_Mis","MED_Wt_Mis")), byrow = TRUE)
  for(j in 1:4)
  {   
    MED_Simulation[,j]=unlist(sapply(1:nsim,function(i)c(EstAll[[s]][[i]][j])))
  }
  MED_Actual=MED[[model_orig]]
  MED_Summary<-data.frame(MED=c(MED_Simulation),MED_Class=rep(colnames(MED_Simulation),each=nsim))
  p10 <- ggplot(MED_Summary, aes(x = MED_Class, y = MED,fill=MED_Class)) +geom_boxplot(colour = "black",position = position_dodge(1), width = 0.4,alpha=0.7) +scale_y_continuous(name = "MED",limits=c(0, 1)) +ggtitle("Boxplot of MED under different methods")+geom_hline(yintercept = MED_Actual)+scale_fill_discrete(breaks=levels(MED_Summary$MED_Class),labels=levels(MED_Summary$MED_Class))+theme(legend.text=element_text(size=7),legend.title = element_text(colour="blue", size=8, face="bold"),axis.text.x = element_text(angle=90, hjust=1,size=8))
  M1<-M2<-M3<-vector()
  
  for(i in 1:4)
  {M1[i]<-mean((MED_Summary[which(MED_Summary[,2]==levels(MED_Summary$MED_Class)[i]),][,1]-MED_Actual)/MED_Actual,na.rm=TRUE)
  M2[i]<-median((MED_Summary[which(MED_Summary[,2]==levels(MED_Summary$MED_Class)[i]),][,1]-MED_Actual)/MED_Actual,na.rm=TRUE)
  M3[i]<-IQR((MED_Summary[which(MED_Summary[,2]==levels(MED_Summary$MED_Class)[i]),][,1]-MED_Actual)/MED_Actual,na.rm=TRUE)
  }
  R_i<-cbind(M1,M2,M3)
  
  R_i<-round(R_i,2)
  colnames(R_i)<-c("mean(R_i)","median(R_i)","IQR(R_i)")
  rownames(R_i)<-levels(MED_Summary$MED_Class)
  R_i<-as.data.frame(R_i)
  R_i<-cbind(Method=rownames(R_i),R_i)
  stable.p <- ggtexttable(R_i,rows=NULL, theme = ttheme(colnames.style = colnames_style(color = "black", size = 9,fill = "grey80", linewidth = 0.5, linecolor = "white", parse = FALSE),base_size = 6.5))
  text.p <- ggparagraph(text = "Relative Devation MED(R_i)=(MED_est-MED_actual)/MED_actual", face = "italic", size = 12, color = "black")
  r1<-ggpar(p10,orientation="horizontal", main = "Simulated From: sigEmax, Fitted on: sigEmax(ft,wt), linear(ft,wt)",font.main = c(10,"bold.italic", "red"),font.x = c(8, "bold", "#2E9FDF"),font.y = c(8, "bold", "#E7B800"),font.xtickslab = 8,font.ytickslab=8)
  grid.arrange(r1,stable.p,nrow=2,as.table=TRUE,heights=c(7,3))
}