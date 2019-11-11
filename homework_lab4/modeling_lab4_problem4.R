library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(broom)
library(tidyverse)

# Bottomly Eset Data
load(file="bottomly.Rdata")

pdata=pData(bottomly.eset)
dim(pdata)
head(pdata)

# Edata conversion to matrix, log transformation and row means
# I've used Aleksander Mućk's hints
edata <- as.matrix(exprs(bottomly.eset))
edata <- log2(as.matrix(edata) + 1)
edata <- edata[rowMeans(edata) > 10, ]

mod = model.matrix(~as.factor(strain), data=pdata)


# Linear model with technical variables
mod = lm(edata[1,] ~ as.factor(pdata$strain) + as.factor(pdata$lane.number))
print(mod)

mod = lm(t(edata) ~ as.factor(pdata$strain) + as.factor(pdata$lane.number))
names(mod)

dim(mod$coefficients)
rownames(mod$coefficients)

library(broom)
mod_tidy <- tidy(mod)

# Using ComBat to clean a dataser
batch = pdata$lane.number
combat_edata = ComBat(dat=edata, batch=pdata$lane.number, mod=model.matrix(~1, data=pdata), par.prior=TRUE, prior.plots=TRUE)
class(combat_edata)

# Estimating surrogate variables (SVs)
mod = model.matrix(~as.factor(strain),data=pdata)
mod0 = model.matrix(~1, data=pdata)
sva_output = sva(edata,mod,mod0,n.sv=2)
head(sva_output$sv)
summary(lm(sva_output$sv ~ pdata$lane.number))

# Fitting LM with surrogate variables
modsva = lm(t(edata) ~ as.factor(pdata$strain) + sva_output$sv)
modsva_tidy <- tidy(modsva)

modcombat = lm(t(combat_edata) ~ as.factor(pdata$strain))
modcombat_tidy <- tidy(modcombat)

est_compare <- tibble(
  LinearModel = mod_tidy %>% filter(term == "as.factor(pdata$strain)DBA/2J") %>% select("estimate") %>% unlist,
  ComBat = modcombat_tidy %>% filter(term == "as.factor(pdata$strain)DBA/2J") %>% select("estimate") %>% unlist,
  SVA = modsva_tidy %>% filter(term == "as.factor(pdata$strain)DBA/2J") %>% select("estimate") %>% unlist)

# problem 4.2 - Linear Model
ggplot(est_compare, aes(x=LinearModel, y=SVA)) +
  geom_point(col="darkgrey", alpha=.5, size=.5) + geom_abline(intercept=0, slope=1, col="darkred") + geom_smooth(method = "lm", se = TRUE)  + theme_bw()

# problem4.2 - ComBat 
ggplot(est_compare, aes(x=ComBat, y=SVA)) +
  geom_point(col="darkgrey", alpha=.5, size=.5) + geom_abline(intercept=0, slope=1, col="darkred") + geom_smooth(method = "lm", se = TRUE)  + theme_bw()

# Check pvalue
ggplot(modsva_tidy %>% filter(term == "as.factor(pdata$strain)DBA/2J")) + geom_histogram(aes(x=p.value), bins = 100, fill="darkorange")
