#################R Practical Script###################
###################Danel Tompsett#####################

#First define your working directory
#This is where we should store the datasets and any R packages
#So that R will know where to find them

WorkDirec<-getwd()

#Lets load in the example datasets. We first need the R package "haven"
#Which we install and load here
#This can be tempramental, so please ask if you need assistance

install.packages("haven",lib=WorkDirec)
library(haven)

#Now lets read in the wide form of the data and look at it

datawide<-as.data.frame(read_dta("Practical2_contY_wide.dta"))
attach(datawide)

head(datawide)
summary(datawide)

##Lets look specifically at the outcome

summary(datawide$Y)

hist(datawide$Y)

##Quantile=Quantile plot. If normal should be a straight line.
qqnorm(datawide$Y)


##Lets check the treatment

table(A0=datawide$A0,A1=datawide$A1)

#Those that sustain or adhere to treatment are those on the diagonal
#Adherence is fairly low, about 56 percent

##A standard analysis of Y on treatment can be done with glm

summary(glm(Y~A0+L0+A1+L1,data=datawide,family="gaussian"))



###################Inverse Probability of Treatment Weighting

#We need to create the stabilised IPTW weights at each time
#To do this we fit glm models and take predicted valued

#Time 0

num1<-glm(A0~1,family=binomial,data=datawide)$fitted
num1[A0==0]<-1-num1[A0==0]

den1<-glm(A0~L0,family=binomial,data=datawide)$fitted
den1[A0==0]<-1-den1[A0==0]


#Standard IPTW weight at time zero is
iptw.0<-1/den1

#Stabilised weight is
stab.0<-num1/den1

#Time 1

num2<-glm(A1~1+A0,family=binomial,data=datawide)$fitted
num2[A1==0]<-1-num2[A1==0]

den2<-glm(A1~L1+A0+L0,family=binomial,data=datawide)$fitted
den2[A1==0]<-1-den2[A1==0]

iptw.1<-1/den2
stab.1<-num2/den2

##The IPTW weights are now the product of these at each time

iptw<-iptw.0*iptw.1
stab<-stab.0*stab.1

summary(iptw)
summary(stab)

#These weights are quite stable, in practice if extreme weights exist
#You could truncate them or set weights with extreme values to zero
# for example iptw[abs(iptw)>5.5]<-0 but not needed here.


##The MSM for the ITT effect is then 
glm(Y~A0,data=datawide,weights=stab.0)
#Note we only need to weights for time 0

## The Sustained Treatment Effect for both time periods is
glm(Y~A0+A1,data=datawide,weights=stab)


##########################Per Protocol Analysis#######################

#Adherence Weights
#First establish those who adhered to treatment


datawide$Adhere<-as.integer(ifelse(datawide$A0==datawide$A1,TRUE,FALSE))


glm(Adhere~1,family=binomial,data=datawide)
#NOTE: By our definition, everyone adheres at time 0
#Thus the adherence weights at time 0 are all 1
#Now predict the probability of adhering at time 1


den2<-glm(Adhere~L0+A0+L1,family=binomial,data=datawide)$fitted

#Adherence weights are now

adh<-1*(datawide$Adhere/den2)



#We now redo the MSM fit, but multiplying the iptw and adherence weights


glm(Y~A0,data=datawide,weights=stab.0*adh)








#####################Extra Analysis: Confidence Interval######
##You can also fit the MSM in long format, we need packages
install.packages("geepack",lib=WorkDirec)
library(geepack)
install.packages("broom",lib=WorkDirec)
library(broom)
#Add our weights to the data and then set it to long format
datawide$iptw<-iptw
datawide$stab<-stab


#rename the variables 
datawidenew<-datawide
names(datawidenew)<-c("id","A.0","L.0","A.1","L.1","Y","iptw","stab")
datalongipw<-reshape(datawidenew,direction="long",
varying=c("A.0","A.1","L.0","L.1"))

#We can now fit the data using generalised estimating equations
#Take the data at the first time period
datat0<-datalongipw[datalongipw$time==0,]

geeITT<-geeglm(Y~A,data=datat0,id=id,weights=stab.0)

#Why do it this way? We can use cluster robust standard errors
#to get a confidence interval on the causal estimands.

summary(geeITT)

tidy(geeITT, conf.int=TRUE)


geePP<-geeglm(Y~A,data=datat0,id=id,weights=stab.0*adh)


summary(geePP)

tidy(geePP, conf.int=TRUE)


####################ITT: G Computation###################


data<-datawide
#We first generate the Qmodel for time 0 
#This describes the relationship between Y, A0 and L0.
Qmodel<-glm(Y~1+A0+L0,family=gaussian,data=data)

#Now we take predictions from this model when we set A0=1
data1<-data
data1$A0<-1

#Our predicted values for Y[A0=1]
Y1pred<-predict(Qmodel,newdata=data1,type="response")

#Repeat this for when A0=0
data0<-data
data0$A0<-0


Y0pred<-predict(Qmodel,newdata=data0,type="response")

#Now take the average difference
summary(Y1pred)
summary(Y0pred)

mean(Y1pred)-mean(Y0pred)


#Standard errors and Confidence Intervals can be ontained by bootstrap
#We do not discuss this here for simplicity

#######################Extra Analysis: Sustained effect###########
#For the sustained efefct we need to pedict values for L1 for when
#A0=0 and A0=1
#Estimate L1~A0 relationship


L1model<-glm(L1~1+A0+L0,family=gaussian,data=data)
dataA0<-data
dataA0$A0<-0

data$L1A0<-predict(L1model,newdata=dataA0,type="response")


dataA1<-data
dataA1$A0=1

data$L1A1<-predict(L1model,newdata=dataA1,type="response")

#Define the Full Q model for Y
#Predict# 
Qmodel<-glm(Y~1+A0+A1+L0+L1,family=gaussian,data=data)

#Now predict Y(1,1) with A0=A1 and L1 at the predicted value
#when A0=1

data11<-data
data11$A0<-1
data11$A1<-1
data11$L1<-data$L1A1

Y11pred<-predict(Qmodel,newdata=data11,type="response")

#Repeat for Y(0,0)


data00<-data
data00$A0<-0
data00$A1<-0
data00$L1<-data$L1A0

Y00pred<-predict(Qmodel,newdata=data00,type="response")

#Now take the expected difference
summary(Y11pred)
summary(Y00pred)

mean(Y11pred)-mean(Y00pred)













