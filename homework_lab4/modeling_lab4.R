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
library(gplots)

# Bladder Cancer Gene Expression
data(bladderdata)

pheno = pData(bladderEset)
edata = exprs(bladderEset)

dim(pheno)
dim(edata)

edata[1:5,1:10]

# Dive into data about samples
names(pheno)
head(pheno)


# problem1
# I've used http://bioconductor.org/packages/release/bioc/vignettes/BatchQC/inst/doc/BatchQC_examples.html
# The table is in Goździewska_problem1.1.pdf but I added all plots, because this function generates very intersting data
batch <- pheno$batch  
condition <- pheno$cancer
batchQC(edata, batch=batch, condition=condition, view_report=TRUE, interactive=TRUE)

# Linear model with technical variables
mod = lm(edata[1,] ~ as.factor(pheno$cancer) + as.factor(pheno$batch))
print(mod)

mod = lm(t(edata) ~ as.factor(pheno$cancer) + as.factor(pheno$batch))
names(mod)

dim(mod$coefficients)
rownames(mod$coefficients)

library(broom)
mod_tidy <- tidy(mod)
ggplot(mod_tidy) + geom_histogram(aes(x=estimate), bins = 100, fill="darkorange")

mod_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer")

ggplot(mod_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer")) + geom_histogram(aes(x=estimate), bins = 100, fill="darkorange")

ggplot(mod_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer")) + geom_histogram(aes(x=p.value), bins = 100, fill="darkorange")


#problem2.1
cluster <- kmeans(edata, centers = 5)

my_palette <- colorRampPalette(c("blue", "white", "green", "darkorange", "purple"))(n = 299)

# afterwards I convert png to pdf 
png("Goździewska_problem2.2.png",height=700,width=700)
heatmap.2(as.matrix(cluster$centers),
          main = "Bladder Cancer Data Clustered", 
          trace="none",        
          margins =c(12,9),    
          col=my_palette,       
          dendrogram="none",     
          scale = "row")
dev.off()

# Empirical Bayes
batch = pheno$batch
combat_edata = ComBat(dat=edata, batch=pheno$batch, mod=model.matrix(~1, data=pheno), par.prior=TRUE, prior.plots=TRUE)

class(combat_edata)
dim(combat_edata)
combat_edata[1:10,1:10]

library(gplots)
library(RColorBrewer)
my_palette <- colorRampPalette(c("blue", "white", "darkred"))(n = 299)

png("bladder_combat.png",height=700,width=700)
heatmap.2(combat_edata,
          main = "Bladder Cancer Data Cleaned by ComBat", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="none",     # only draw a row dendrogram
          scale = "row")
dev.off()

# problem2.2
# Clustering the rows of the data AFTER combat
cluster_combat <- kmeans(combat_edata, centers = 5)

my_palette <- colorRampPalette(c("blue", "white", "green", "darkorange", "purple"))(n = 299)

# afterwards I convert png to pdf 
png("Goździewska_problem2.2.png",height=700,width=700)
heatmap.2(as.matrix(cluster_combat$centers),
          main = "Bladder Cancer Data Clustered", 
          trace="none",        
          margins =c(12,9),    
          col=my_palette,       
          dendrogram="none",     
          scale = "row")
dev.off()


# problem3
pcor <- cor(edata, method = "pearson")

png("./Goździewska_problem3.png",height=700,width=700)
heatmap.2(as.matrix(pcor),
          main = "Bladder Cancer Data Clustered", 
          trace="none",        
          margins =c(12,9),    
          col=my_palette,       
          dendrogram="none",     
          scale = "row")
dev.off()

#########################################################

modcombat = lm(t(combat_edata) ~ as.factor(pheno$cancer))
library(broom)
modcombat_tidy <- tidy(modcombat)
ggplot(modcombat_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer")) + geom_histogram(aes(x=estimate), bins = 100, fill="darkorange")

ggplot(modcombat_tidy, aes(estimate, term)) +
  geom_point() +
  geom_vline(xintercept = 0)

est_compare <- tibble(
  LinearModel = mod_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer") %>% select("estimate") %>% unlist,
  ComBat = modcombat_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer") %>% select("estimate") %>% unlist)

ggplot(est_compare, aes(x=LinearModel, y=ComBat)) +
  geom_point(col="darkgrey", alpha=.5, size=.5) + geom_abline(intercept=0, slope=1, col="darkred") + geom_smooth(method = "lm", se = TRUE)  + theme_bw()

ggplot(modcombat_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer")) + geom_histogram(aes(x=p.value), bins = 100, fill="darkorange")


# Surrogate Variable Analysis SVA
mod = model.matrix(~as.factor(cancer), data=pheno)
num.sv(edata,mod,method="be")
num.sv(edata,mod,method="leek")


# Estimating surrogate variables SVs
mod = model.matrix(~as.factor(cancer),data=pheno)
mod0 = model.matrix(~1, data=pheno)
sva_output = sva(edata,mod,mod0,n.sv=2)

head(sva_output$sv)
summary(lm(sva_output$sv ~ pheno$batch))


# Visualizing and exploring SVs
sva_batch <- tibble(SV1=sva_output$sv[,1],SV2=sva_output$sv[,2],batch=as.factor(pheno$batch), cancer=as.factor(pheno$cancer),outcome=as.factor(pheno$outcome))
ggplot(sva_batch) + geom_point(aes(x=SV1,y=SV2, col=batch))

ggplot(sva_batch) + geom_point(aes(x=SV1,y=SV2, col=cancer))

ggplot(sva_batch) + geom_point(aes(x=SV1,y=SV2, col=outcome))

sva_batch <- tibble(sv1 = sva_output$sv[,1], sv2 = sva_output$sv[,2], batch = as.factor(pheno$batch))
sva_batch <- gather(sva_batch,"sv","value",-batch)

ggplot(sva_batch) + geom_violin(aes(x=batch,y=value)) + geom_jitter(aes(x=batch,y=value,col=batch)) + facet_wrap(~ sv, ncol = 1)


# Fitting a LM with surrogate variables
modsva = lm(t(edata) ~ as.factor(pheno$cancer) + sva_output$sv)
modsva_tidy <- tidy(modsva)

est_compare <- tibble(
  LinearModel = mod_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer") %>% select("estimate") %>% unlist,
  ComBat = modcombat_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer") %>% select("estimate") %>% unlist,
  SVA = modsva_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer") %>% select("estimate") %>% unlist)

ggplot(est_compare, aes(x=LinearModel, y=SVA)) +
  geom_point(col="darkgrey", alpha=.5, size=.5) + geom_abline(intercept=0, slope=1, col="darkred") + geom_smooth(method = "lm", se = TRUE)  + theme_bw()

ggplot(est_compare, aes(x=ComBat, y=SVA)) +
  geom_point(col="darkgrey", alpha=.5, size=.5) + geom_abline(intercept=0, slope=1, col="darkred") + geom_smooth(method = "lm", se = TRUE)  + theme_bw()

ggplot(modsva_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer")) + geom_histogram(aes(x=p.value), bins = 100, fill="darkorange")

pvalues <- tibble(
  LinearModel = mod_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer") %>% select("p.value") %>% unlist,
  SVA = modsva_tidy %>% filter(term == "as.factor(pheno$cancer)Cancer") %>% select("p.value") %>% unlist)

ggplot(est_compare, aes(x=LinearModel, y=SVA)) +
  geom_point(col="darkgrey", alpha=.5, size=.5) + geom_abline(intercept=0, slope=1, col="darkred") + geom_smooth(method = "lm", se = TRUE)  + theme_bw()

