#Author: David Lowry

#This R script runs the full model and the additve model and compares
#them to identify eQTLxtreatment interactions. Based on section 7.2
# of "A guide to QTL mapping with R/qtl" by Broman and Sen.

#load r/QTL
library(qtl)

#Load modules that allow for parallization across 
library(snow)
library(rlecuyer)

#Input data file
data <- read.cross("csv", ".", "rQTL_May_SNP_stacked_5_2_13.csv")

#Make backcross formated data into RIL
class(data)[1]="riself"

#Jitter closeby markers
data<-jittermap(data)

#Load data file with covariates
cyto_block<-read.csv(file="block_cyto_treatment.csv", header=TRUE)

#Extract covariates
block1<-as.vector(cyto_block$block1)
block2<-as.vector(cyto_block$block2)
block3<-as.vector(cyto_block$block3)
block4<-as.vector(cyto_block$block4)
block5<-as.vector(cyto_block$block5)
cyto<-as.vector(cyto_block$cyto)
treatment<-as.vector(cyto_block$treatment)

#Combine covariates that will be fit as additive effects in the model
x<-cbind(block1, block2, block3, block4, block5, cyto, treatment)

#Imput genotypes
dataIM<-calc.genoprob(data, step=1, map.function="kosambi")

#Below is the modelling fitting. I had to break up into 2500 expression phenotypes at a time
#because larger sizes cause R to crash due to memory issues

#Fit additive model
out.a=scanone(dataIM, pheno.col=1:2500, addcov=x, model=c("normal"), method=c("hk"))

#Fit interactive model
out.u=scanone(dataIM, pheno.col=1:2500, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"))

#Set seed for permutaitons
set.seed(54955149)

#Permutations for the additive model
operm.a <- scanone(dataIM, pheno.col=1:2500, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)

set.seed(54955149)

#Permutations for the interactive model
operm.u <- scanone(dataIM, pheno.col=1:2500, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)

#Compare LOD profiles across models
out.u <- c(out.u, out.u-out.a, labels=c("f","i"))
operm.u <- cbind(operm.u, operm.u - operm.a, labels=c("f","i"))

#Output results of QTL mapping in tabByCol format
tabbycol_2500<-summary(out.u, perms=operm.u, format="tabByCol", pvalues=TRUE, ci.function="lodint")
write.table(tabbycol_2500, file = "Tab_by_SNP_treatment_2500.txt", quote=FALSE, row.name=FALSE)

#Repeat process for the next 2500 expression phenotypes

out.a=scanone(dataIM, pheno.col=2501:5000, addcov=x, model=c("normal"), method=c("hk"))
out.u=scanone(dataIM, pheno.col=2501:5000, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"))
set.seed(54955149)
operm.a <- scanone(dataIM, pheno.col=2501:5000, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
set.seed(54955149)
operm.u <- scanone(dataIM, pheno.col=2501:5000, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
out.u <- c(out.u, out.u-out.a, labels=c("f","i"))
operm.u <- cbind(operm.u, operm.u - operm.a, labels=c("f","i"))
tabbycol_5000<-summary(out.u, perms=operm.u, format="tabByCol", pvalues=TRUE, ci.function="lodint")

write.table(tabbycol_5000, file = "Tab_by_SNP_treatment_5000.txt", quote=FALSE, row.name=FALSE)

out.a=scanone(dataIM, pheno.col=5001:7500, addcov=x, model=c("normal"), method=c("hk"))
out.u=scanone(dataIM, pheno.col=5001:7500, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"))
set.seed(54955149)
operm.a <- scanone(dataIM, pheno.col=5001:7500, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
set.seed(54955149)
operm.u <- scanone(dataIM, pheno.col=5001:7500, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
out.u <- c(out.u, out.u-out.a, labels=c("f","i"))
operm.u <- cbind(operm.u, operm.u - operm.a, labels=c("f","i"))
tabbycol_7500<-summary(out.u, perms=operm.u, format="tabByCol", pvalues=TRUE, ci.function="lodint")

write.table(tabbycol_7500, file = "Tab_by_SNP_treatment_7500.txt", quote=FALSE, row.name=FALSE)

out.a=scanone(dataIM, pheno.col=7501:10000, addcov=x, model=c("normal"), method=c("hk"))
out.u=scanone(dataIM, pheno.col=7501:10000, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"))
set.seed(54955149)
operm.a <- scanone(dataIM, pheno.col=7501:10000, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
set.seed(54955149)
operm.u <- scanone(dataIM, pheno.col=7501:10000, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
out.u <- c(out.u, out.u-out.a, labels=c("f","i"))
operm.u <- cbind(operm.u, operm.u - operm.a, labels=c("f","i"))
tabbycol_10000<-summary(out.u, perms=operm.u, format="tabByCol", pvalues=TRUE, ci.function="lodint")

write.table(tabbycol_10000, file = "Tab_by_SNP_treatment_10000.txt", quote=FALSE, row.name=FALSE)

out.a=scanone(dataIM, pheno.col=10001:12500, addcov=x, model=c("normal"), method=c("hk"))
out.u=scanone(dataIM, pheno.col=10001:12500, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"))
set.seed(54955149)
operm.a <- scanone(dataIM, pheno.col=10001:12500, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
set.seed(54955149)
operm.u <- scanone(dataIM, pheno.col=10001:12500, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
out.u <- c(out.u, out.u-out.a, labels=c("f","i"))
operm.u <- cbind(operm.u, operm.u - operm.a, labels=c("f","i"))
tabbycol_12500<-summary(out.u, perms=operm.u, format="tabByCol", pvalues=TRUE, ci.function="lodint")

write.table(tabbycol_12500, file = "Tab_by_SNP_treatment_12500.txt", quote=FALSE, row.name=FALSE)

out.a=scanone(dataIM, pheno.col=12501:15000, addcov=x, model=c("normal"), method=c("hk"))
out.u=scanone(dataIM, pheno.col=12501:15000, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"))
set.seed(54955149)
operm.a <- scanone(dataIM, pheno.col=12501:15000, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
set.seed(54955149)
operm.u <- scanone(dataIM, pheno.col=12501:15000, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
out.u <- c(out.u, out.u-out.a, labels=c("f","i"))
operm.u <- cbind(operm.u, operm.u - operm.a, labels=c("f","i"))
tabbycol_15000<-summary(out.u, perms=operm.u, format="tabByCol", pvalues=TRUE, ci.function="lodint")

write.table(tabbycol_15000, file = "Tab_by_SNP_treatment_15000.txt", quote=FALSE, row.name=FALSE)

out.a=scanone(dataIM, pheno.col=15001:17500, addcov=x, model=c("normal"), method=c("hk"))
out.u=scanone(dataIM, pheno.col=15001:17500, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"))
set.seed(54955149)
operm.a <- scanone(dataIM, pheno.col=15001:17500, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
set.seed(54955149)
operm.u <- scanone(dataIM, pheno.col=15001:17500, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
out.u <- c(out.u, out.u-out.a, labels=c("f","i"))
operm.u <- cbind(operm.u, operm.u - operm.a, labels=c("f","i"))
tabbycol_17500<-summary(out.u, perms=operm.u, format="tabByCol", pvalues=TRUE, ci.function="lodint")

write.table(tabbycol_17500, file = "Tab_by_SNP_treatment_17500.txt", quote=FALSE, row.name=FALSE)

out.a=scanone(dataIM, pheno.col=17501:20000, addcov=x, model=c("normal"), method=c("hk"))
out.u=scanone(dataIM, pheno.col=17501:20000, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"))
set.seed(54955149)
operm.a <- scanone(dataIM, pheno.col=17501:20000, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
set.seed(54955149)
operm.u <- scanone(dataIM, pheno.col=17501:20000, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
out.u <- c(out.u, out.u-out.a, labels=c("f","i"))
operm.u <- cbind(operm.u, operm.u - operm.a, labels=c("f","i"))
tabbycol_20000<-summary(out.u, perms=operm.u, format="tabByCol", pvalues=TRUE, ci.function="lodint")

write.table(tabbycol_20000, file = "Tab_by_SNP_treatment_20000.txt", quote=FALSE, row.name=FALSE)

out.a=scanone(dataIM, pheno.col=20001:22500, addcov=x, model=c("normal"), method=c("hk"))
out.u=scanone(dataIM, pheno.col=20001:22500, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"))
set.seed(54955149)
operm.a <- scanone(dataIM, pheno.col=20001:22500, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
set.seed(54955149)
operm.u <- scanone(dataIM, pheno.col=20001:22500, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
out.u <- c(out.u, out.u-out.a, labels=c("f","i"))
operm.u <- cbind(operm.u, operm.u - operm.a, labels=c("f","i"))
tabbycol_22500<-summary(out.u, perms=operm.u, format="tabByCol", pvalues=TRUE, ci.function="lodint")

write.table(tabbycol_22500, file = "Tab_by_SNP_treatment_22500.txt", quote=FALSE, row.name=FALSE)

out.a=scanone(dataIM, pheno.col=22501:25662, addcov=x, model=c("normal"), method=c("hk"))
out.u=scanone(dataIM, pheno.col=22501:25662, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"))
set.seed(54955149)
operm.a <- scanone(dataIM, pheno.col=22501:25662, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
set.seed(54955149)
operm.u <- scanone(dataIM, pheno.col=22501:25662, addcov=x, intcovar=treatment, model=c("normal"), method=c("hk"), n.perm=1000, n.cluster=6)
out.u <- c(out.u, out.u-out.a, labels=c("f","i"))
operm.u <- cbind(operm.u, operm.u - operm.a, labels=c("f","i"))
tabbycol_25662<-summary(out.u, perms=operm.u, format="tabByCol", pvalues=TRUE, ci.function="lodint")

write.table(tabbycol_25662, file = "Tab_by_SNP_treatment_25662.txt", quote=FALSE, row.name=FALSE)
