#=========================== LDA2 - Non-Linear PK Modeling ============
#====== CAMPBELL MCDULING
#====== 06 JUNE 2023
rm(list = ls())
setwd("~/OneDrive - University of Cape Town/2023/MSc Biostatistics/Coursework Semester 1/LDA/A2")

# LOAD PACKAGES
library(dplyr)
library(nlme)
library(psych)
library(ggplot2)
library(lattice)
# library(ggeffects)
# library(emmeans)
library(ggpubr)

#==================== DATA MANAGEMENT =============
dat0 = read.csv("malariadata.csv")
head(dat0)
summary(dat0)

describe(dat0)
str(dat0)

#define variable types accordingly
dat0$site = as.factor(dat0$site)
dat0$arm = as.factor(dat0$arm)
dat0$pid = as.factor(dat0$pid)
dat0$gender = as.factor(dat0$gender)
dat0$country = as.factor(dat0$country)
dat0$PIoutcome = as.factor(dat0$PIoutcome)
summary(dat0)

names(dat0)[8] = "Pconc" #shorten names for convenience
names(dat0)[9] = "Sconc"

#create logged drug concentrations
dat0$l.Sconc = log(dat0$Sconc+0.000001)
dat0$l.Pconc = log(dat0$Pconc+0.000001)

#create weight and age quartiles
quart <- quantile(dat0$weight, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = T)
dat0$weight_quart <- cut(dat0$weight, breaks = quart, labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = T)
quart <- quantile(dat0$age, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = T)
dat0$age_quart <- cut(dat0$age, breaks = quart, labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = T)


#==================== DATA VALIDATION =============

#--------- create wide df
dat0.wide = reshape(dat0, timevar = "pday",  idvar="pid",direction="wide", v.names = c("Pconc", "Sconc"), 
                    drop=c("Hb", "pardens", "gamedens", "Studyyear"))

#REMOVE ALL PARTICIPANTS WITH NON-ZERO VISIT 1 DRUG CONC
#BOTH DRUGS
index.zer = which(dat0.wide$Pconc.0 == 0 & dat0.wide$Sconc.0 == 0) #retain only those with zero for both conc at visit 0
pid.zer = dat0.wide[index.zer,]$pid
index.zer.long = which(dat0$pid %in% pid.zer)
dat0 = dat0[index.zer.long,] #remove from long DF
dat0.wide = dat0.wide[index.zer,] #remove non-zeros from wide DF
n_distinct(dat0$pid) # SAMPLE SIZE REDUCED to 296 participants
n_distinct(dat0.wide$pid)


#--------- MISSING VALUES
#missing drug concentration values
sum(is.na(dat0$Pconc)) #308 missing Pconc observations
n_distinct(dat0[is.na(dat0$Pconc),]$pid) #over 180 subjects 
sum(is.na(dat0$Sconc)) #474 missing Sconc observations
n_distinct(dat0[is.na(dat0$Sconc),]$pid) #over 216 subjects 

#--------- MISSING BY TIME POINT

#missings per time point
sum(is.na(dat0.wide$Sconc.0)) 
sum(is.na(dat0.wide$Sconc.2)) 
sum(is.na(dat0.wide$Sconc.3)) 
sum(is.na(dat0.wide$Sconc.7)) 
sum(is.na(dat0.wide$Sconc.14))
sum(is.na(dat0.wide$Sconc.21)) 
sum(is.na(dat0.wide$Sconc.28)) 
sum(is.na(dat0.wide$Sconc.42))
sum(is.na(dat0.wide$Pconc.0)) 
sum(is.na(dat0.wide$Pconc.2)) 
sum(is.na(dat0.wide$Pconc.3)) 
sum(is.na(dat0.wide$Pconc.7)) 
sum(is.na(dat0.wide$Pconc.14))
sum(is.na(dat0.wide$Pconc.21)) 
sum(is.na(dat0.wide$Pconc.28)) 
sum(is.na(dat0.wide$Pconc.42))


#--------- which participants have very few drug concentration measurements
#REMOVE ALL WHO HAVE LESS THAN 4 measurements for each drug
index.missmoreP = which(rowSums(!is.na(dat0.wide[,seq(14, 31, 2)]))<4) #those with less than four S after base
dat0.wide[index.missmoreP,] 
index.missmoreS = which(rowSums(!is.na(dat0.wide[,seq(15, 31, 2)]))<4) #those less than four P after base
dat0.wide[index.missmoreS,] 

index.missmore = c(index.missmoreP, index.missmoreS)

pid.missmore = dat0.wide[index.missmore,]$pid
index.missmore.long = which(dat0$pid %in% pid.missmore)
dat0 = dat0[-index.missmore.long,]
dat0.wide = dat0.wide[-index.missmore,]
n_distinct(dat0$pid) #sample size reduced to 271
n_distinct(dat0.wide$pid)


dat0 %>% filter(is.na(weight))
dat0 = dat0 %>% filter(!is.na(weight)) #remove this participant - no weight measurement
dat0.wide = dat0.wide %>% filter(!is.na(weight))
n_distinct(dat0$pid)
n_distinct(dat0.wide$pid) # sample size reduced to 270

# write.csv(dat0, file="dat0.csv", row.names = F)
# write.csv(dat0.wide, file="dat0.wide.csv", row.names = F)

#==================== DATA EXPLORATION =============
#exploratory data:for descriptive stats - as original 

dat0 = read.csv("dat0.csv")
dat0.wide = read.csv("dat0.wide.csv")

dat0$site = as.factor(dat0$site)
dat0$arm = as.factor(dat0$arm)
dat0$pid = as.factor(dat0$pid)
dat0$gender = as.factor(dat0$gender)
dat0$weight_quart = as.factor(dat0$weight_quart)
dat0$age_quart = as.factor(dat0$age_quart)
summary(dat0)

#------- DESCRIPTIVE STATS (gender, age, weight, site, arm, study year)

#baseline statistics by study arm
dat.baseline = dat0 %>% filter(pday==0) 

#continuous variables
library(purrr)
dat.baseline[,c(2,8, 9, 11, 12)] %>% split(.$arm) %>% map(summary)

#discrete variables
dat.baseline[,c(1,2, 10)] %>% split(.$arm) %>% map(summary)

#descriptive plots - continuous covars by categorical
ggplot(dat.baseline, aes(x=age, col=arm)) + geom_density()
ggplot(dat.baseline, aes(x=weight, col=arm)) + geom_density()
ggplot(dat.baseline, aes(x=age, col=gender)) + geom_density()
ggplot(dat.baseline, aes(x=weight, col=gender)) + geom_density()

#correlations between continuous variables
#cor.plot(dat.baseline[,c(5:9, 11:12)])

#------- OUTCOMES (all zero at baseline)
### BIVARIATE
dat.baseline[,c(2,8,9, 17,18)] %>% split(.$arm) %>% map(summary) #baseline drug vs arm

#dist of drug conc disregarding timepoint
hist(dat0$Pconc)
hist(dat0$Sconc)
ggplot(dat0, aes(x=Pconc, fill=arm)) + geom_histogram() + theme_minimal() 
ggplot(dat0, aes(x=Sconc, fill=arm)) + geom_histogram() + theme_minimal()
table(dat0$arm)

#by visit - densities 
p1=ggplot(dat0.wide, aes(x=Pconc.1, col=arm)) + geom_density() + theme_minimal() 
p2=ggplot(dat0.wide, aes(x=Pconc.2, col=arm)) + geom_density() + theme_minimal() 
p3=ggplot(dat0.wide, aes(x=Pconc.3, col=arm)) + geom_density() + theme_minimal() 
p4=ggplot(dat0.wide, aes(x=Pconc.7, col=arm)) + geom_density() + theme_minimal() 
p5=ggplot(dat0.wide, aes(x=Pconc.14, col=arm)) + geom_density() + theme_minimal() 
p6=ggplot(dat0.wide, aes(x=Pconc.21, col=arm)) + geom_density() + theme_minimal() 
p7=ggplot(dat0.wide, aes(x=Pconc.28, col=arm)) + geom_density() + theme_minimal() 
p8=ggplot(dat0.wide, aes(x=Pconc.42, col=arm)) + geom_density() + theme_minimal() 
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=4,nrow=2)

s1=ggplot(dat0.wide, aes(x=Sconc.1, col=arm)) + geom_density() + theme_minimal() 
s2=ggplot(dat0.wide, aes(x=Sconc.2, col=arm)) + geom_density() + theme_minimal() 
s3=ggplot(dat0.wide, aes(x=Sconc.3, col=arm)) + geom_density() + theme_minimal() 
s4=ggplot(dat0.wide, aes(x=Sconc.7, col=arm)) + geom_density() + theme_minimal() 
s5=ggplot(dat0.wide, aes(x=Sconc.14, col=arm)) + geom_density() + theme_minimal() 
s6=ggplot(dat0.wide, aes(x=Sconc.21, col=arm)) + geom_density() + theme_minimal() 
s7=ggplot(dat0.wide, aes(x=Sconc.28, col=arm)) + geom_density() + theme_minimal() 
s8=ggplot(dat0.wide, aes(x=Sconc.42, col=arm)) + geom_density() + theme_minimal() 
ggarrange(s1,s2,s3,s4,s5,s6,s7,s8,ncol=4,nrow=2)



#---------------------- LONGITUDINAL EXPLORATORY ----
#===SPAGHETTI PLOTS

xyplot(Sconc ~ pday, dat0, type="l", groups = pid) 
#by study arm
xyplot(Sconc ~ pday|arm, dat0, type="l", groups = pid) 
#by study site
xyplot(Sconc ~ pday|site, dat0, type="l", groups = pid) 

xyplot(Pconc ~ pday, dat0, type="l", groups = pid) 
#by study arm
xyplot(Pconc ~ pday|arm, dat0, type="l", groups = pid) 
#by study site
xyplot(Pconc ~ pday|site, dat0, type="l", groups = pid) 



#mean log S concentrations
ggplot(dat0, aes(x = pday, y = l.Sconc, color = arm, group = arm)) +
  stat_summary(fun = "mean", geom = "line") +
  labs(x = "Day", y = "Mean log(S + 1e-06)", color = "Study Arm") +
  theme_minimal()

#Mean drug conc by study arm (both drugs)
ggplot(dat0, aes(x = pday, color = arm, group=arm)) +
  geom_line(aes(y = Sconc, linetype = "Sconc"), stat = "summary", fun = "mean") +
  geom_line(aes(y = Pconc, linetype = "Pconc"), stat = "summary", fun = "mean") +
  labs(x = "Day", y = "Mean drug concentrations", color = "Study Arm", linetype = "Drug") +
  scale_linetype_manual(values = c("Sconc" = "solid", "Pconc" = "dashed")) +
  theme_minimal()
#above but on log scale
ggplot(dat0, aes(x = pday, color = arm, group=arm)) +
  geom_line(aes(y = l.Sconc, linetype = "log(Sconc)"), stat = "summary", fun = "mean") +
  geom_line(aes(y = l.Pconc, linetype = "log(Pconc)"), stat = "summary", fun = "mean") +
  labs(x = "Day", y = "Mean log drug concentrations", color = "Study Arm", linetype = "Drug") +
  scale_linetype_manual(values = c("log(Sconc)" = "solid", "log(Pconc)" = "dashed")) +
  theme_minimal()
# BY SEX
ggplot(dat0, aes(x = pday, color = gender, group=gender)) +
  geom_line(aes(y = Sconc, linetype = "Sconc"), stat = "summary", fun.y = "mean") +
  geom_line(aes(y = Pconc, linetype = "Pconc"), stat = "summary", fun.y = "mean") +
  labs(x = "Day", y = "Mean drug concentrations", color = "Sex", linetype = "Drug") +
  scale_linetype_manual(values = c("Sconc" = "solid", "Pconc" = "dashed")) +
  theme_minimal()
ggplot(dat0, aes(x = pday, color = gender, group=gender)) +
  geom_line(aes(y = l.Sconc, linetype = "log(Sconc)"), stat = "summary", fun.y = "mean") +
  geom_line(aes(y = l.Pconc, linetype = "log(Pconc)"), stat = "summary", fun.y = "mean") +
  labs(x = "Day", y = "Mean log drug concentrations", color = "Sex", linetype = "Drug") +
  scale_linetype_manual(values = c("log(Sconc)" = "solid", "log(Pconc)" = "dashed")) +
  theme_minimal()
# BY SITE
ggplot(dat0, aes(x = pday, color = site, group=site)) +
  geom_line(aes(y = Sconc, linetype = "Sconc"), stat = "summary", fun.y = "mean") +
  geom_line(aes(y = Pconc, linetype = "Pconc"), stat = "summary", fun.y = "mean") +
  labs(x = "Day", y = "Mean drug concentrations", color = "Site", linetype = "Drug") +
  scale_linetype_manual(values = c("Sconc" = "solid", "Pconc" = "dashed")) +
  theme_minimal()
ggplot(dat0, aes(x = pday, color = site, group=site)) +
  geom_line(aes(y = l.Sconc, linetype = "log(Sconc)"), stat = "summary", fun.y = "mean") +
  geom_line(aes(y = l.Pconc, linetype = "log(Pconc)"), stat = "summary", fun.y = "mean") +
  labs(x = "Day", y = "Mean log drug concentrations", color = "Site", linetype = "Drug") +
  scale_linetype_manual(values = c("log(Sconc)" = "solid", "log(Pconc)" = "dashed")) +
  theme_minimal()

#by weight quartile
ggplot(dat0, aes(x = pday, color = weight_quart, group=weight_quart)) +
  geom_line(aes(y = Sconc, linetype = "Sconc"), stat = "summary", fun = "mean") +
  geom_line(aes(y = Pconc, linetype = "Pconc"), stat = "summary", fun = "mean") +
  labs(x = "Day", y = "Mean drug concentrations", color = "Weight quartile", linetype = "Drug") +
  scale_linetype_manual(values = c("Sconc" = "solid", "Pconc" = "dashed")) +
  theme_minimal()


# VARIANCE PROFILES
ggplot(dat0, aes(x = pday, color = arm, group=arm)) +
  geom_line(aes(y = Sconc, linetype = "Sconc"), stat = "summary", fun = "var") +
  geom_line(aes(y = Pconc, linetype = "Pconc"), stat = "summary", fun = "var") +
  labs(x = "Day", y = "Variance of drug concentrations", color = "Study Arm", linetype = "Drug") +
  scale_linetype_manual(values = c("Sconc" = "solid", "Pconc" = "dashed")) +
  theme_minimal()
ggplot(dat0, aes(x = pday, color = arm, group=arm)) +
  geom_line(aes(y = l.Sconc, linetype = "Sconc"), stat = "summary", fun = "var") +
  geom_line(aes(y = l.Pconc, linetype = "Pconc"), stat = "summary", fun = "var") +
  labs(x = "Day", y = "Variance of logged drug concentrations", color = "Study Arm", linetype = "Drug") +
  scale_linetype_manual(values = c("Sconc" = "solid", "Pconc" = "dashed")) +
  theme_minimal()

#============================= MODELING ========================
dat0 = read.csv("dat0.csv")
dat0.wide = read.csv("dat0.wide.csv")

dat0$site = as.factor(dat0$site)
dat0$arm = as.factor(dat0$arm)
dat0$pid = as.factor(dat0$pid)
dat0$gender = as.factor(dat0$gender)
dat0$weight_quart = as.factor(dat0$weight_quart)
dat0$age_quart = as.factor(dat0$age_quart)


dat0$Dose = 1 #create dose variable

#what is minimum conc, besides zero
temp = dat0 %>% filter(.$Pconc != 0)
min(temp$Pconc) #0.146
temp = dat0 %>% filter(.$Sconc != 0)
min(temp$Sconc) #0.125

dat0$Sconc[which(dat0$Sconc==0)] = 1e-4 #change all zero concentrations to very small +ive number
dat0$Pconc[which(dat0$Pconc==0)] = 1e-4

dat0.grp.s = groupedData(Sconc ~ pday|site/pid, dat0)
dat0.grp.p = groupedData(Pconc ~ pday|site/pid, dat0)

#------------- FIRST ORDER COMPARTMENT MODEL
#====================== SULFCONC=======
#----- LEVEL 1 MODEL
#fit for one subject
tst = dat0.grp.s[dat0.grp.s$pid=="MOB2004_015",]
nls(Sconc~SSfol(Dose, pday, lKe, lKa, lCl), data=tst)

#fit for all subjects
#contr = nls.control(minFactor = 1/5000, tol=1e-2)
fol1.lis = nlsList(Sconc ~ SSfol(Dose, pday, lKe, lKa, lCl)|pid, data=dat0.grp.s,
                    na.action = na.exclude)
fol1.lis
summary(fol1.lis)

#investigate random effects structure
plot(intervals(fol1.lis)) #much between-pid variability in clearance and elimination, absorption is similar
#investigate covariance structure
pairs(fol1.lis, id=0.05) #maybe negative relationship between lKe and lKa, positive between lKe and lCl
#residuals
plot(fol1.lis, pid~resid(.), abline=0, xlab="residuals")
densityplot(resid(fol1.lis), xlab="residuals")
plot(resid(fol1.lis))

#----- LEVEL 2 MODEL
##### RANDOM EFFECTS
nlm.contr = nlmeControl(minScale = 1e-20, maxIter = 50, tolerance = 1e-1)

#Model 1a - general +ive def
ran.S1a = list(lKe + lKa + lCl ~ 1)
nlm.S1a = nlme(fol1.lis,  random = ran.S1a, control = nlm.contr, method="REML")
summary(nlm.S1a) #magnitude of lKa and lCl RE are much smaller than lKe

#Model 1b - diagonal
ran.S1b = pdDiag(lKe + lKa + lCl ~ 1)
nlm.S1b = update(nlm.S1a, random=ran.S1b) 
summary(nlm.S1b) #magnitude of lKe is much smaller now

anova(nlm.S1a, nlm.S1b) #loglik is signif reduced, take diagonal model

# #Model 1c - blocked
# ran.S1c = pdBlocked(list(lKe~1, lKa+lCl~1))
# nlm.S1c = update(nlm.S1a, random=ran.S1c) #does not converge
# summary(nlm.S1c)
#
# anova(nlm.S1b, nlm.S1c)

#Model 1d - no elimination RE
ran.S1d = pdDiag(lKa + lCl ~ 1)
nlm.S1d = update(nlm.S1b, random=ran.S1d)
summary(nlm.S1d)

anova(nlm.S1b, nlm.S1d) # no sig diff between models, take simpler model (1d)

#Model 1e - only lKa RE
nlm.S1e = update(nlm.S1d, random=list(lKa~1))
# #Model 1f - only lCl RE
# nlm.S1f = update(nlm.S1d, random=list(lCl~1)) #singularity error

anova(nlm.S1d, nlm.S1e) #sig diff between models, stick with model 1d with two REs


plot(nlm.S1d) #non-constant variance
qqnorm(nlm.S1d) #big deviation from normality for residuals
densityplot(resid(nlm.S1d))
qqnorm(nlm.S1d, ~ ranef(.)) #normality for random effects

##### FIXED EFFECTS
#investigate fixed eff structure
re.nlm.S1d = ranef(nlm.S1d)
re.nlm.S1d$arm = dat0.wide$arm
re.nlm.S1d$gender = dat0.wide$gender
re.nlm.S1d$weight = dat0.wide$weight
re.nlm.S1d$age = dat0.wide$age
plot(re.nlm.S1d, form=lKa~arm+gender) #no significant deviation between arms or sexes
plot(re.nlm.S1d, form=lKa~age+weight) #appears to be a weight and age effect, weight is stronger
plot(re.nlm.S1d, form=lCl~arm+gender) #no significant deviation between arms or sexes
plot(re.nlm.S1d, form=lCl~age+weight) #appears to be a weight and age effect, strong weight effect
cor(dat0.wide$age, dat0.wide$weight, method="pearson")


# BUILD MODEL
nlm.contr2 = nlmeControl(minScale = 1e-20, tolerance = 0.2)

nlm.S1d.ML = update(nlm.S1d, method="ML")#create ML version for model comparisons 
summary(nlm.S1d.ML)
#extract fixed effects for starting values
nlm.S1d.fix = fixef(nlm.S1d.ML)
fix1 = list(lKe + lKa + lCl ~ arm)
start1 = c(nlm.S1d.fix[1], -0.02, nlm.S1d.fix[2], 0.04, nlm.S1d.fix[3], 0.08)

#Model 2a - fixed effect on study arm
nlm.S2a = update(nlm.S1d.ML, fixed=fix1, start=start1, control=nlm.contr2)
summary(nlm.S2a) #arm effects are sig


anova(nlm.S1d.ML, nlm.S2a) #logLik significantly reduced

#Model 2b - adjust for other covariates - forward selection
nlm.S2a.fix = fixef(nlm.S2a)
fix2 = list(lKe + lKa + lCl ~ arm+weight)
start2 = c(nlm.S2a.fix[1], nlm.S2a.fix[2], 0,  nlm.S2a.fix[3], nlm.S2a.fix[4], 0,
           nlm.S2a.fix[5], nlm.S2a.fix[6], 0)
nlm.S2b = update(nlm.S2a, fixed=fix2, start=start2)
summary(nlm.S2b) #weight significantly modifies elimination and maybe clearance, not absorption

anova(nlm.S2a, nlm.S2b) #significantly different, 2b is chosen
anova(nlm.S1d.ML, nlm.S2a, nlm.S2b)

#look at RE associations now that weight and arm are adjusted for
re.nlm.S2b = ranef(nlm.S2b)
re.nlm.S2b$arm = dat0.wide$arm
re.nlm.S2b$gender = dat0.wide$gender
re.nlm.S2b$weight = dat0.wide$weight
re.nlm.S2b$age = dat0.wide$age
plot(re.nlm.S2b, form=lKa.(Intercept)~arm+gender) #much the same as before
plot(re.nlm.S2b, form=lKa.(Intercept)~age+weight) #much of the association across both age and weight has been removed
plot(re.nlm.S2b, form=lCl.(Intercept)~arm+gender) #much the same as before
plot(re.nlm.S2b, form=lCl.(Intercept)~age+weight) #still association remaining

#Model 2c - try weight^2 also
fix3 = list(lKe + lKa + lCl ~ arm+weight+I(weight^2))
start3 = c(nlm.S2a.fix[1], nlm.S2a.fix[2], 0, 0, nlm.S2a.fix[3], nlm.S2a.fix[4], 0, 0,
           nlm.S2a.fix[5], nlm.S2a.fix[6], 0,0)
nlm.S2c = update(nlm.S2b, fixed=fix3, start=start3)
#nlm.S2c = nlme(Sconc ~ SSfol(Dose, pday, lKe, lKa, lCl), data=dat0.grp.s,
#               fixed=fix3, start=start3, random = ran.S1d, control=nlm.contr2, na.action = na.omit)
summary(nlm.S2c)

anova(nlm.S2b, nlm.S2c) #model significantly improved

anova(nlm.S2a, nlm.S2b, nlm.S2c)

#Model 2d - has weight effects only on elim and clear
fix4 = list(lKa~arm, lKe+lCl~arm+weight+I(weight^2))
nlm.S2c.fix = fixef(nlm.S2c)
start4 = c(nlm.S2c.fix[1], nlm.S2c.fix[2], nlm.S2c.fix[3], nlm.S2c.fix[4], 
           nlm.S2c.fix[5], nlm.S2c.fix[6], 
           nlm.S2c.fix[9], nlm.S2c.fix[10], nlm.S2c.fix[11], nlm.S2c.fix[12])
#nlm.S2d = update(nlm.S2c, fixed=fix4, start=start4) #singularity error

#=========== MODEL VALIDATION
#fit model with REML
nlm.S2c = update(nlm.S2c, method="REML")
summary(nlm.S2c)
re.nlm.S2c = ranef(nlm.S2c)
re.nlm.S2c$arm = dat0.wide$arm
re.nlm.S2c$gender = dat0.wide$gender
re.nlm.S2c$weight = dat0.wide$weight
re.nlm.S2c$age = dat0.wide$age

plot(re.nlm.S2c, form=lKa.(Intercept)~age+weight) 
plot(re.nlm.S2c, form=lCl.(Intercept)~age+weight) #association removed

#residuals
plot(nlm.S2c) #definitely heteroskedasticity
qqnorm(nlm.S2c) #definitely non-normal
densityplot(resid(nlm.S2c), xlab="Standardized residuals")

#random effects
qqnorm(nlm.S2c, ~ ranef(.)) #normality for random effects
d1=densityplot(unlist(ranef(nlm.S2c)[1]), xlab="lKa random effects")
d2=densityplot(unlist(ranef(nlm.S2c)[2]), xlab="lCl random effects")
ggarrange(d1,d2)
plot(ranef(nlm.S2c), xlab="")


#### MODEL VARIANCE
nlm.S3a = update(nlm.S2c, weights=varExp())
summary(nlm.S3a)
plot(nlm.S3a)
qqnorm(nlm.S3a, ~resid(.))
densityplot(resid(nlm.S3a))
plot(nlm.S2c)

#nlm.S3b = update(nlm.S2c, weights=varPower()) #error thrown

#nlm.S3c = update(nlm.S2c, weights=varConstPower()) #error thrown

nlm.S3d = update(nlm.S2c, weights=varIdent(~pday))
summary(nlm.S3d)
plot(nlm.S3d) #not fixed

#====================== PYRCONC =======
#----- LEVEL 1 MODEL
#fit for one subject
tst = dat0.grp.p[dat0.grp.p$pid=="MOB2004_015",]
nls(Pconc~SSfol(Dose, pday, lKe, lKa, lCl), data=tst)

#fit for all subjects
#lis.contr = nls.control(maxiter = 200, minFactor = 1/1024, tol=1e-2)
fol2.lis = nlsList(Pconc ~ SSfol(Dose, pday, lKe, lKa, lCl)|pid, data=dat0.grp.p,
                   na.action = na.exclude)
fol2.lis

#investigate random effects structure
plot(intervals(fol2.lis)) #much between-pid variability in all three parameters
plot(fol2.lis, pid~resid(.), abline=0)

#investigate covariance structure
pairs(fol2.lis, id=0.05) #maybe negative relationship between lKe and lKa

#----- LEVEL 2 MODEL
##### RANDOM EFFECTS
nlm.contr = nlmeControl(minScale = 1e-6, maxIter = 50, tol=1e-1)

#Model 1 - diagonal with all three parameters as RE
ran.P1 = pdDiag(lKa + lKe + lCl ~ 1)
#nlm.P1 = nlme(fol2.lis, random=ran.P1, control=nlm.contr, method = "REML") #cannot fit all three - singularity error (unidentifiable model)


#Model 1a - diagonal with two as RE
ran.P1a = pdDiag(lKa + lKe  ~ 1) 
nlm.P1a = nlme(fol2.lis, random = ran.P1a, control = nlm.contr, method = "REML")
summary(nlm.P1a) #magnitude of lKa RE is small

#Model 1b 
ran.P1b = pdDiag(lKe + lCl  ~ 1) 
nlm.P1b = update(nlm.P1a, random=ran.P1b)  
summary(nlm.P1b)

anova(nlm.P1a, nlm.P1b) #take model 1b on information criteria and mag of REs

#Model 1c
ran.P1c = pdDiag(lKa +  lCl ~ 1)
#nlm.P1c = update(nlm.P1a, random=ran.P1c)  #does not converge
#summary(nlm.P1c)


#Model 1c - +ive definite
ran.P1d = list(lKe + lCl ~ 1)
nlm.P1d = update(nlm.P1a, random=ran.P1d)
summary(nlm.P1d)

anova(nlm.P1a,nlm.P1b, nlm.P1d) # sig diff, take simpler model

#Model 1e - RE only on elimination
ran.P1e = list(lKe~1)
nlm.P1e = update(nlm.P1a, random=ran.P1e)
summary(nlm.P1e)



# #Model 1g - RE on clearance only
# ran.P1g = list(lCl~1)
# nlm.P1g = update(nlm.P1a, random=ran.P1g) #singularity error
# summary(nlm.P1g)


anova(nlm.P1a, nlm.P1b, nlm.P1d, nlm.P1e) #best model is 1b



plot(nlm.P1b) #fairly random scatter
qqnorm(nlm.P1b) #big deviation from normality for residuals
densityplot(resid(nlm.P1b))
qqnorm(nlm.P1b, ~ ranef(.)) #reasonable normality for random effects
densityplot(unlist(ranef(nlm.P1b)[1]), xlab="lKe")
densityplot(unlist(ranef(nlm.P1b)[2]), xlab="lKe")


##### FIXED EFFECTS
#investigate fixed eff structure
re.nlm.P1b = ranef(nlm.P1b)
re.nlm.P1b$arm = dat0.wide$arm
re.nlm.P1b$gender = dat0.wide$gender
re.nlm.P1b$weight = dat0.wide$weight
re.nlm.P1b$age = dat0.wide$age
plot(re.nlm.P1b, form=lKe~arm+gender) #no significant deviation between arms or sexes
plot(re.nlm.P1b, form=lKe~age+weight) #appears to be a weight and age effect, strong weight effect
plot(re.nlm.P1b, form=lCl~age+weight) #appears to be a weight and age effect, strong weight effect
plot(re.nlm.P1b, form=lCl~arm+gender) #appears to be a weight and age effect, strong weight effect

# BUILD MODEL
nlm.P1b.ML = update(nlm.P1b, method="ML")#create ML version for model comparisons 
summary(nlm.P1b.ML)
#extract fixed effects for starting values
nlm.P1b.fix = fixef(nlm.P1b.ML)
fix1 = list(lKe + lKa + lCl ~ arm)
start1 = c(nlm.P1b.fix[1], 0, 2, 0, nlm.P1b.fix[3], 0)
nlm.contr2 = nlmeControl(minScale = 1e-20, tol=0.1) #reduce minimum halving factor

#Model 2a - fixed effect on study arm
nlm.P2a = update(nlm.P1b.ML, fixed=fix1, start=start1, control=nlm.contr2)
summary(nlm.P2a) #none of additional effects are sig


anova(nlm.P1b.ML, nlm.P2a) 

#adjust for other covariates - forward selection
#weight
nlm.P2a.fix = fixef(nlm.P2a)
fix2 = list(lKe + lKa + lCl ~ arm+weight)
start2 = c(nlm.P2a.fix[1], nlm.P2a.fix[2], 0, nlm.P2a.fix[3], nlm.P2a.fix[4], 0, 
           nlm.P2a.fix[5], nlm.P2a.fix[6], 0)
#start2 = c(-1, -0.02,0.001, 2, -0.3, 0.0001, -7, 0.15, 0.0001)
#nlm.P2b = update(nlm.P1b.ML,  fixed=fix2, start=start2, control=nlm.contr2) #AAAAAHAHHHHHHHHHHHHH SINGULARITY AHHHHHHH
#summary(nlm.P2b) #ANGRY, fury, mania

#anova(nlm.P1b.ML, nlm.P2a, nlm.P2b) #frustration



#look at RE associations now that weight and arm are adjusted for
# re.nlm.P2b = ranef(nlm.P2b)
# re.nlm.P2b$arm = dat0.wide$arm
# re.nlm.P2b$gender = dat0.wide$gender
# re.nlm.P2b$weight = dat0.wide$weight
# re.nlm.P2b$age = dat0.wide$age
# plot(re.nlm.P2b, form=lKe.(Intercept)~arm+gender) #much the same as before
# plot(re.nlm.P2b, form=lKe.(Intercept)~age+weight) #the association across both age and weight has been removed
# plot(re.nlm.P1b, form=lKe~age+weight)  #before correction

#=========== MODEL VALIDATION
nlm.P2a.reml = update(nlm.P2a, method="REML")
summary(nlm.P2a.reml)
#residuals
plot(nlm.P2a.reml) #non-constant variance
qqnorm(nlm.P2a.reml) #large deviation from normality
densityplot(resid(nlm.P2a.reml), xlab="Standardized residuals")

#random effects
qqnorm(nlm.P2a.reml, ~ ranef(.)) #normality for random effects
d1=densityplot(unlist(ranef(nlm.P2a.reml)[1]), xlab="lKe random effects")
d2=densityplot(unlist(ranef(nlm.P2a.reml)[2]), xlab="lCl random effects")
ggarrange(d1, d2)
plot(ranef(nlm.P2a.reml))


#---------------------- MODEL VARIANCE ---------
varmod1 = varPower()
# nlm.P3a = update(nlm.P2a.reml, weights=varPower()) #singular
# nlm.P3b = update(nlm.P2a.reml, weights=varExp()) #singular
# nlm.P3c = update(nlm.P2a.reml, weights=varConstPower()) #singular

