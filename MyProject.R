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

#AICc from the MuMIn package
anole.aic <- AICc(anole.lm,anole.allo)

#aicw from the geiger package
anole.aicw <- aicw(anole.aic$AICc)

print(anole.aicw)

logLik(anole.lm)

## 'log Lik.' -148.5536 (df=3)

logLik(anole.allo)

## 'log Lik.' -147.4431 (df=3)

anole.log%>%
  ggplot(aes(HTotal,SVL,col=Ecomorph2))+geom_point()+geom_smooth(method="lm")

#model incuding ecomorph as variable:
anole.log.eco.lm <- lm(HTotal~SVL*Ecomorph2,anole.log)
summary(anole.log.eco.lm)
anova(anole.log.eco.lm)

#is a model that includes ecomorph just a good one, or is it better than without ecomorph?
anole.log.lm  <- lm(HTotal~SVL,anole.log)
anova(anole.log.lm)
anole.log.aic <- AICc(anole.log.lm,anole.log.eco.lm)
aicw(anole.log.aic$AICc)

anole.log <- anole.log %>%
  mutate(res=residuals(anole.log.lm))
anole.log%>%
  ggplot(aes(Ecomorph2,res))+geom_point()
#boxplot fo today
p.eco <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=res)) +geom_boxplot()
print(p.eco)

p.eco+ geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)#incoudes mean residual

anole.tree <- read.tree("anole.tre")
plot(anole.tree,cex=0.4) #our phylo tree
#PGLS under BM, no ecomorph
pgls.BM1 <- gls(HTotal ~SVL, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#PGLS under BM, w ecomorph
pgls.BM2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")


#PGLS under OU, no ecomorph
pgls.OU1 <- gls(HTotal ~SVL, correlation = corMartins(0,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#PGLS under OU, w, ecomorph
pgls.OU2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corMartins(0,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

anole.phylo.aic <- AICc(pgls.BM1,pgls.BM2,pgls.OU1,pgls.OU2)
aicw(anole.phylo.aic$AICc)
#consider whether Ecomorph is a significant factor in predicting the hinglimb-SVL 
#relationship in a phylogenetic context by preforming an ANOVA on our best fitting model, pgsl.BM2

anova(pgls.BM2)

anole.log <- anole.log%>%
  mutate(phylo.res=residuals(pgls.BM2))

p.eco.phylo <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=phylo.res)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)

print(p.eco.phylo)
#faceting
anole.log%>%
  dplyr::select(Ecomorph2,res,phylo.res)%>%
  pivot_longer(cols=c("res","phylo.res"))%>%
  print%>%
  ggplot(aes(x=Ecomorph2,y=value)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)+facet_grid(name~.,scales = "free_y")+ylab("residual")

#For this project report, write a script named “groupname_module3.R” and undertake the following operations.

#1 Combine the code above so that you can establish the anole.log data tibble.
#2 Using the log-transformed data, construct two simple linear models that assess the effect 
  #of perch diameter and height by including these as covariates in your models. Be sure to 
  #use + notation rather than *, assuming there is no interaction (there isn’t, trust me!).
#3 Explore how both perch diameter and height effect the hindlimb-SVL relationship by plotting the residuals of your simple linear models against these covariates. This will require mutating a data tibble to include residuals from both models. Please produce two separate plots.
#4 Under a BM model of trait evolution and using the tree provided, construct phylogenetic 
  #least squares models of the hindlimb-SVL relationships that include the unique 
  #combinations of these two covariates, i.e,

    #* A PGLS model with the hindlimb-SVL relationship + perch height
    #* A PGLS model with the hindlimb-SVL relationship + perch diameter
    #* A PGSL model with the hindlimb-SVL relationship + perch height + perch diameter
    #* 
#5 Assess the fit of each of these three models using AICc and AICw and comment on 
  #(with comments in the script) whether one or both of the covariates is a significant 
  #predictor of hinglimb length in a phylogenetic context.
#6 Produce a plot of your own design that concisely visualizes the effect of your covariate(s) 
  #on the hindlimb residuals of the best fitting PGLS model.
#7 Submit your script to the link. Submissions are due by 11:59 PM on Sunday, September 27th.
