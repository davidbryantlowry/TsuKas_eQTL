#Author: David Lowry

#This R script runs the additive model only for eQTL mapping



library(qtl)
data <- read.cross("csv", ".", "rQTL_May_SNP_sum_5_13_13.csv")
class(data)[1]="riself"
data<-jittermap(data)
cyto_block<-read.csv(file="block_cyto.csv", header=TRUE)
block1<-as.vector(cyto_block$block1)
block2<-as.vector(cyto_block$block2)
block3<-as.vector(cyto_block$block3)
block4<-as.vector(cyto_block$block4)
block5<-as.vector(cyto_block$block5)
cyto<-as.vector(cyto_block$cyto)
x<-cbind(block1, block2, block3, block4, block5, cyto)
dataIM<-calc.genoprob(data, step=1, map.function="kosambi")

transcriptLOD1=scanone(dataIM, pheno.col=1:5000, addcov=x, model=c("normal"), method=c("hk"))
transcriptLOD2=scanone(dataIM, pheno.col=5001:10000, addcov=x, model=c("normal"), method=c("hk"))
transcriptLOD3=scanone(dataIM, pheno.col=10001:15000, addcov=x, model=c("normal"), method=c("hk"))
transcriptLOD4=scanone(dataIM, pheno.col=15001:20000, addcov=x, model=c("normal"), method=c("hk"))
transcriptLOD5=scanone(dataIM, pheno.col=20001:25662, addcov=x, model=c("normal"), method=c("hk"))
transcriptLODall=cbind(transcriptLOD1, transcriptLOD2, transcriptLOD3, transcriptLOD4, transcriptLOD5)


set.seed(123456)
maxlod1=scanone(dataIM, pheno.col=1:5000, model=c("normal"), method=c("hk"), n.perm=1000, addcov=x)

tabbycol_5000<-summary(transcriptLOD1, format="tabByCol", pvalues=TRUE, perm=maxlod1, ci.function="lodint")

write.table(tabbycol_5000, file = "Tab_by_5000_Mckay_sum_SNP.txt", quote=FALSE, row.name=FALSE)

set.seed(123456)
maxlod2=scanone(dataIM, pheno.col=5001:10000, model=c("normal"), method=c("hk"), n.perm=1000, addcov=x)

tabbycol_10000<-summary(transcriptLOD2, format="tabByCol", pvalues=TRUE, perm=maxlod2, ci.function="lodint")

write.table(tabbycol_10000, file = "Tab_by_10000_Mckay_sum_SNP.txt", quote=FALSE, row.name=FALSE)

set.seed(123456)
maxlod3=scanone(dataIM, pheno.col=10001:15000, model=c("normal"), method=c("hk"), n.perm=1000, addcov=x)

tabbycol_15000<-summary(transcriptLOD3, format="tabByCol", pvalues=TRUE, perm=maxlod3, ci.function="lodint")

write.table(tabbycol_15000, file = "Tab_by_15000_Mckay_sum_SNP.txt", quote=FALSE, row.name=FALSE)


set.seed(123456)
maxlod4=scanone(dataIM, pheno.col=15001:20000, model=c("normal"), method=c("hk"), n.perm=1000, addcov=x)

tabbycol_20000<-summary(transcriptLOD4, format="tabByCol", pvalues=TRUE, perm=maxlod4, ci.function="lodint")

write.table(tabbycol_20000, file = "Tab_by_20000_Mckay_sum_SNP.txt", quote=FALSE, row.name=FALSE)

set.seed(123456)
maxlod5=scanone(dataIM, pheno.col=20001:25662, model=c("normal"), method=c("hk"), n.perm=1000, addcov=x)

tabbycol_25662<-summary(transcriptLOD5, format="tabByCol", pvalues=TRUE, perm=maxlod5, ci.function="lodint")

write.table(tabbycol_24756, file = "Tab_by_25662_Mckay_sum_SNP.txt", quote=FALSE,       row.name=FALSE)
