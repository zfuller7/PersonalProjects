###The goal of this document is to demonstrate different methods for survival analysis. 
####A concise presentation of results will be in a separate document.

#set the working directory
setwd()


#####Load data and libraries

library(dplyr)
library(readr)
library(survival)
library(survminer)
library(ggplot2)

#note that I set the working directory with setwd(). 
data <- read_csv("cirrhosis.csv")


###Status: C=Censored, CL=Censored due to liver transplant, D=Death

#####In order to make the code easily reproducible, use this block to set parameters. You will need to change labels on graphs for different data.

#time till event
time <- data$N_Days
#name of status column in dataset
status <- data$Status
#event of interest (in this example, death)
s <- "D"

#gender 
female <- "F"
male <- "M"


#####Kaplan-Meier Approach, without covariates

#Parameters to set
#time till event
time <- data$N_Days
#name of status column in dataset
status <- data$Status
#event of interest (in this example, death)
s <- "D"

#Calculations
data$SurvObj <- with(data, Surv(time, status == s))
km_surv <- survfit(SurvObj ~ 1, data = data, conf.type = "log-log")

#Standard Plot
plot(km_surv, col=c("red", "blue", "blue"), lwd=2,
     main="Cirrhosis patient Survival (Kaplan-Meier)",
     xlab="Survival Time in Days",
     ylab = "Survival Probability")

####Kaplan-Meier using ggsurvplot

ggsurvplot(
  km_surv,
  data=data,
  size=1,
  palette="navy",
  conf.int=TRUE,
  ggtheme= theme_minimal(),
  xlab="Time in Days")+
  labs(title = "Cirrhosis Patient Survival (Kaplan-Meier)")


####Survival Curves by Treatment
#####Using ggsurvplot only

bytreat <- survfit(SurvObj~Drug, data=data)
ggsurvplot(bytreat, data=data, conf.int=TRUE, pval=TRUE)+labs(title="Cirrhosis Survival by Treatment")

###From the plot, we see that gender is not a significant predictor of survival probability. At least not without including covariates.
###The p-value shown is calculated from the log-rank test and is used to compare the two curves.The issue here is that the curves cross, so log-rank may not perform well in this case.

####Use non-parametric test for the difference between the curves.

survdiff(SurvObj~Drug, data=data, rho=0)


###The p-value from the non-parametric test, still shows that there is not a significant difference between the curves.

####Find a model with significant covariates

#reduce data to complete cases.
data <- data[complete.cases(data),]
mdl <- survreg(SurvObj~.-N_Days-Status-ID, data=data)
set.seed(1234)
mdl_reduced <- capture.output(step(mdl)) #capture.output is just blocking the iterations of the step function from printing.
mdl_reduced[187:215]


####Now, compare a model with all significant covariates to one with all plus drug.

#
drug_model <- survreg(SurvObj~Age+Edema+Bilirubin+Albumin+Copper+SGOT+Prothrombin+Stage+Drug, data=data)

nodrug_model <- survreg(SurvObj~Age+Edema+Bilirubin+Albumin+Copper+SGOT+Prothrombin+Stage, data=data)

anova(nodrug_model, drug_model)

###What we see here is that even when significant covariates and included, the drug does not have a significant effect on survival outcome.

#Data Downloaded from 
#fedesoriano. (August 2021). Cirrhosis Prediction Dataset. Retrieved [Date Retrieved] from https://www.kaggle.com/fedesoriano/cirrhosis-prediction-dataset.