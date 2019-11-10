library(devtools)
library(Biobase)
library(BatchQC)
library(sva)
library(bladderbatch)
library(broom)
library(tidyverse)
library(edgeR)
library(Rtsne)
library(data.table)

# Bottomly Eset Data
load(file="bottomly.Rdata")

pdata=pData(bottomly.eset)
dim(pdata)
head(pdata)

edata=as.matrix(exprs(bottomly.eset))
dim(edata)
edata[1:5,1:5]


# Linear model with technical variables
# problem4.1
mod = lm(edata[1,] ~ as.factor(pdata$experiment.number) + as.factor(pdata$strain))
print(mod)

mod = lm(t(edata) ~ as.factor(pdata$experiment.number) + as.factor(pdata$strain))
names(mod)

dim(mod$coefficients)
rownames(mod$coefficients)

library(broom)
mod_tidy <- tidy(mod)
ggplot(mod_tidy) + geom_histogram(aes(x=estimate), bins = 100, fill="darkorange")

# problem4.2
batch = pdata$strain
combat_edata = ComBat(dat=edata, batch=pdata$strain, mod=model.matrix(~1, data=pdata), par.prior=TRUE, prior.plots=TRUE)
