#Our tasks will be to:
  
#1) Visualize how hind limb length varies with size (i.e., allometry!).
#2)Visualize and asses how hind limb length vs. size relationship covary with ecological niche.
#3)Learn more complex operations in ggplot than we’ve undertaken so far.
#4)Evaluate this hindlimb-size relationship using standard frequentist models within 
#  and without a phylogenetic context.
#5)Using an information theory approach, assess the fit of phylogenetically corrected 
#  models of hind-limb variation under different modes of character evolution

library(tidyverse) # Remember to load your libraries!
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)


#Remember to set your working directory, too!
setwd("~/Schoolwork Fall 20/Organismal Biology/Module 3 Project")

#load data
anole <- read_csv("anole.dat.csv")
anole.eco <- read_csv("anole.eco.csv")

#must merge anole data tibble w anole.eco tibble
anole2 <- anole%>%
  left_join(anole.eco)%>%
  filter(!Ecomorph%in%c("U","CH"))%>%
  na.omit()%>%
  print()
anole.log <- anole2%>%
  mutate_at(c("SVL", "HTotal","PH","ArbPD"),log)

#does hind limb length vary with size? Let’s assess this initially by 
#plotting hind limb length vs size
anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_smooth(method="lm")

anole.lm <- lm(HTotal~SVL,anole2)
#here where asking R to build a linear model with HTotal predicted by SVL
#The second parameter specifies what data the variables are coming from, in this case, anole2.
coef(anole.lm)

anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_abline(slope=coef(anole.lm)[2],intercept=coef(anole.lm)[1],col="blue")

SVL2 <- seq(min(anole2$SVL),max(anole2$SVL),0.1)

pred.lm <-tibble(
  SVL=SVL2,
  H.pred=predict(anole.lm,newdata = data.frame(SVL=SVL2))
)

anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_point(data=pred.lm,aes(SVL,H.pred),col="blue")

summary(anole.lm)

anole.allo <- nls(HTotal~a*SVL^b, start=list(b=1, a=1),data = anole2)

summary(anole.allo)
#so far, we’ve assessed the null hypothesis that SVL has no effect on HTotal 
#and rejected it based on two statistical tests under different models, 
#one linear, the other allometric.

