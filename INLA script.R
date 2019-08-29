
########### LIBRARIES ###########

library(rgdal)
library(maptools)
library(rgeos)
library(spdep)
library(tmap)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(Rgraphviz)
library(INLA)
library(lme4)
library(psych)
library(arules)


########### IMPORT SHAPFEILE  ###########

setwd("~/Documents/THESIS/DATA/FIJI/Final folder/tempdir")

myshp=readOGR('.', 'dataframe4')

head(myshp@data)

xtabs(~VIL_NAM+TIK_NAM,myshp@data)
xtabs(~TIK_NAM+LMMA, myshp@data)

########### DATA  ###########

mydata=myshp@data

########### DESCRIPTIVE RESULTS BY LMMA  ###########

desc=describeBy(mydata,group = 'LMMA')
desc
write.csv(desc[[1]], file='descNOT_LMMA3.csv')
write.csv(desc[[2]], file='descLMMA3.csv')


########### CHECK DATA  ###########

par(mfrow=c(2,4))
for(i in 6:14){
  hist(mydata[,i], main=names(mydata)[i], breaks=50)
}

round(apply(mydata[,6:14]*1000,2, summary),3)

apply(mydata[,6:14],2, function(x) round(x, 3))  

########### AGGREGATE DATA INTO CATEGORIES  ###########

names(mydata)

mydata$CmmPOI_cat = discretize(mydata$CmmPOI_*1000, breaks=3, method='freq')
mydata$neigh_cat = discretize(mydata$neghbrs*1000, breaks=3, method='freq')
mydata$wtr_cat  =discretize(mydata$wtr_dst*1000, breaks=3, method='freq')
mydata$TrsPOI_cat=discretize(mydata$TrsPOI_*1000, breaks=3, method='freq')
mydata$CORAL_cat=discretize(mydata$CORAL_A,   breaks=3, method='freq')
mydata$MANGROV_cat  =discretize(mydata$MANGROV, breaks=3, method='freq')
mydata$cocnt_cat=discretize(mydata$cocnt_r,breaks=3, method='fixed')
mydata$coralyes=ifelse(mydata$CORAL_A>113,1,0)
mydata$trsdic=ifelse(mydata$TrsPOI_1>43,0,1)

for ( i in 15:21) {
  print(table(mydata[,i]))
  sum(table(mydata[,i]))
}

apply(mydata,2,function(x) sum(is.na(x)))

########### TEST MODELS IN BASIC FREQUENTIST GLM ###########
########### TEST WITH AND WITHOUT RANDOM EFFECT  ###########

###frequentist no random effect
fit1=glm(LMMA~   MrkPOI_cat+ neigh_cat + wtr_cat +TrsPOI_cat+CORAL_cat+ MANGROV_cat+ cocnt_cat,
         family='binomial', data=mydata)

summary(fit1)

backwards = step(fit1,trace = 0) 

fit1b=glmer(LMMA~   MrkPOI_cat+ neigh_cat+ wtr_cat +TrsPOI_cat +CORAL_cat+ MANGROV_cat+ cocnt_cat+(1|VIL_ID),
            family='binomial', data=mydata)

summary(fit1b)
lattice::dotplot(ranef(fit1b, which = "VIL_ID", condVar = TRUE), scales = list(y = list(alternating = 0)))


### frequentist random  effect ( varying intercept) at district
fit2=glmer(LMMA~   MrkPOI_cat+ neigh_cat+ wtr_cat +TrsPOI_cat +CORAL_cat+ MANGROV_cat+ cocnt_cat+(1| TIK_ID),
           family='binomial', data=mydata)
summary(fit2)
ranef(fit2)

fit2b=glmer(LMMA~ CORAL_cat+(1| TIK_ID),
            family='binomial', data=mydata)
summary(fit2b)
ranef(fit2b)

fit2b4=glm(LMMA~ CORAL_cat,
           family='binomial', data=mydata)
summary(fit2b4)

fit2b=glmer(LMMA~ CORAL_cat+ MANGROV_cat+ cocnt_cat+(1| TIK_ID),
            family='binomial', data=mydata)
summary(fit2b)
ranef(fit2b)

fit2b1=glmer(LMMA~ CORAL_cat+ wtr_cat +MANGROV_cat+ cocnt_cat+(1| TIK_ID),
             family='binomial', data=mydata)
summary(fit2b1)
ranef(fit2b1)


fit2b2=glmer(LMMA~ CORAL_cat+ MrkPOI_cat +MANGROV_cat+ cocnt_cat+(1| TIK_ID),
             family='binomial', data=mydata)
summary(fit2b2)
ranef(fit2b2)


fit2c=glmer(LMMA~coralyes+ MANGROV_cat+ cocnt_cat+(1| DSTRCT_I),
           family='binomial', data=mydata)
summary(fit2c)
ranef(fit2c)


lattice::dotplot(ranef(fit2, which = "TIK_ID", condVar = TRUE), scales = list(y = list(alternating = 0)))

### frequentist nested  random  effect  
fit3=glmer(LMMA~   MrkPOI_cat+ neigh_cat+ wtr_cat +TrsPOI_cat +CORAL_cat+ MANGROV_cat+ cocnt_cat+(1|TIK_NAM),
           family='binomial', data=mydata)
summary(fit3)
ranef(fit3)

########### AIC MODEL FIT EVALUATION  ###########

AIC(fit1)
AIC(fit1b)
AIC(fit2)
AIC(fit2b)
AIC(fit2b1)
AIC(fit2b2)
AIC(fit3)


########### INLA BAYESIAN MODELS ###########

mydata$DS_VIL=paste(mydata$TIK_ID, mydata$VIL_ID,sep='')

form1=LMMA~ CmmPOI_cat+ neigh_cat+ wtr_cat +TrsPOI_cat +CORAL_cat+ MANGROV_cat
form1b=LMMA~ CmmPOI_cat+ neigh_cat+ wtr_cat +TrsPOI_cat +CORAL_cat+ MANGROV_cat+f(VIL_ID ,model='iid',constr=T)

form2=LMMA~ CmmPOI_cat+ neigh_cat+ wtr_cat +TrsPOI_cat +coralyes + MANGROV_cat+ f(TIK_NAM,model='iid',constr=T)

formll=list(form1, form1b, form2)


res=list()
for ( i in 1:3){
  res[[i]]=inla(formll[[i]], family="binomial",data=mydata, control.compute = list(waic = TRUE,cpo=TRUE))
}
lapply(res, function(x)  x$waic$waic)


###### according to waic the best fit is for model form2
######## ODDS RATIO COVARIATE ESTIMATES #########

modelres=round(exp(res[[1]]$summary.fixed[,c(1,3,5)]),3)
modelres
write.csv(modelres, file='Model_Estimates_OR3.csv')



