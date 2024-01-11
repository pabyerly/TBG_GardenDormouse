library(adegenet)
library(hierfstat)
library(pegas)
library(vcfR)
library(reshape2)
library(ggplot2)
library(poppr)
library(ape)
library(igraph)
library(dplyr)
library(vegan)
library(qvalue)
library(OutFLANK)  
library(pcadapt)


#palette
#myCol=c("darkgreen",  "slateblue2",  "goldenrod1", "violetred", "plum1", "palegreen2", "tomato", "cadetblue1", "yellow")
myCol=c("mediumseagreen",  "lightsteelblue",  "wheat", "darkolivegreen", "tomato", "rosybrown3", "gray53", "lavender", "sienna3")
#################################################################################################

#adegenet
#convert vcf to genlight 

#input vcf file (vcf: related individuals removed--n = 86, thinned, sorted by sample name)
vcf<- read.vcfR("populations.snps_indfilter_related_names_sort.vcf.gz")
#variant count: 41775

#to genind format
mouse_gen<-vcfR2genind(vcf)
#to genlight
mouse_snp <- vcfR2genlight(vcf)

#read in population file
snp_pop <- read.table("pop_relate.txt", header=T)

#assign population, ie put it in the population slot 
pop(mouse_gen)=snp_pop$pop
#for genlight
pop(mouse_snp)=snp_pop$pop

#####################################################################
#DACP
#NO a priori grouping 
#identify number of groups with k means 
grp=find.clusters(mouse_gen, max.n.clust=5)

#look at group assignment 
table(pop(mouse_gen), grp$grp)

#cross-validation to determine number of principal components to use
set.seed(999)
mouse1 <- xvalDapc(tab(mouse_gen, NA.method = "mean"),grp$grp)

#check results 
mouse1[-1]

#n.da is number of populations - 1 
dapc1 <- dapc(mouse_gen, var.contrib = TRUE, n.pca=10, n.da=3, grp$grp)

#print contents of the object
print.dapc(dapc1)
#summary/useful info 
summary.dapc(dapc1)
#predict individual assignment 
predict.dapc(dapc1)

scatter(dapc1,scree.da=FALSE, scree.pca=TRUE, bg="white", posi.pca="topleft", 
        legend=TRUE, col=myCol, clab=0, cstar=0, cex=2, pch=c(15, 16, 17, 18), solid=0.8)

compoplot(dapc1, col=myCol,lab="", ncol=2)

loadingplot(dapc1$var.contr, threshold=quantile(dapc1$var.contr,0.75))

###############################################################
#DAPC: 
#a priori grouping 
#cross-validation to determine number of principal components to use
set.seed(999)
mouse <- xvalDapc(tab(mouse_gen, NA.method = "mean"), pop(mouse_gen))
#check results 
mouse[-1]

#n.da is number of populations - 1 
dapc <- dapc(mouse_gen, var.contrib = TRUE, n.pca=50, n.da=7, pop(mouse_gen))

#print contents of the object
print.dapc(dapc)
#summary/useful info 
summary.dapc(dapc)
#predict individual assignment 
predict.dapc(dapc)

scatter(dapc,scree.da=FALSE, scree.pca=TRUE, bg="white", posi.pca="bottomright", 
        legend=TRUE, col=myCol, clab=0, cstar=0, cex=2, pch=c(15, 16, 17, 18), solid=0.5)

compoplot(dapc, col=myCol,lab="", ncol=2)

loadingplot(dapc$var.contr, threshold=quantile(dapc2$var.contr,0.75))

#plot of population assignment predict.dapc(dapc1)
assignplot(dapc, subset=1:33)

#################################################################
#PCA

#look at eigenvalues to pick best # PCs
mouse.pca <- glPca(mouse_snp, nf = 3)
#summary
mouse.pca 

#eigenvalues
barplot(100*mouse.pca$eig/sum(mouse.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance explained", line = 2)
title(xlab="Eigenvalues", line = 1)

mouse.pca.scores <- as.data.frame(mouse.pca$scores)
mouse.pca.scores$pop <- pop(mouse_gen)

#plot on first 2 PCs
p <- ggplot(mouse.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=5)
p <- p + stat_ellipse(level = 0, size = 1)
p <- p + scale_color_manual(values = myCol) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_classic()
p

ggsave("PCA.jpg", width= 6, height = 6, dpi = 600 )

#################################################################

#OUTFLANK OUTLIER ANALYSIS 

data<- read.vcfR("populations.snps_indfilter_related_names_sort.vcf.gz")

#extract just the genotypes
geno <- extract.gt(data)
dim(geno)

#change to OutFLANk format
G <- geno  
G[geno %in% c("0/0")] <- 0
G[geno  %in% c("0/1")] <- 1
G[geno %in% c("1/1")] <- 2
G[is.na(G)] <- 9
tG <- t(G)
dim(tG)

#Bring in the metadata
pop <- read.table("pop_relate.txt", header=T) 

#calculate FST (note that this is FST between all populations!)
OF_SNPs <- MakeDiploidFSTMat(tG, locusNames=seq(1, 41775, by=1), popNames=pop$predict)
head(OF_SNPs)
hist(OF_SNPs$FST,breaks=50)
summary(OF_SNPs$FST) 

OF <- OutFLANK(FstDataFrame=OF_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.1, Hmin=0.1, NumberOfSamples=86, qthreshold=0.1)

#plot
OutFLANKResultsPlotter(OF, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.005, titletext=NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(OF, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.40, titletext = NULL)

#Identify outliers 
P1 <- pOutlierFinderChiSqNoCorr(OF_SNPs,Fstbar=OF$FSTNoCorrbar,
                                dfInferred=OF$dfInferred,qthreshold=0.05,Hmin=0.1)
outliers <- P1$OutlierFlag==TRUE #which of the SNPs are outliers?
table(outliers)

#Alt method
outliers <- which(OF_SNPs$results$OutlierFlag=="TRUE")
print(outliers)

######################################################################################

#outliers with PCAdapt

genos <- read.pcadapt("populations.snps_indfilter_related_names_sort.bed",type=c("bed"))
#bring in pop data
pop <- read.table("pop_relate.txt", header=T) 

#find best number of PCs (k)
x <- pcadapt(input=genos,K=20, min.maf = 0.01)
plot(x,option="screeplot")

#2 is best number of K (PCs)
x <- pcadapt(input=genos, K = 2, min.maf = 0.01)

#pca by group (here, predicted population assignment from DAPC)
plot(x,option="scores",pop=pop$pop, col=myCol)
plot(x,option="scores",pop=pop$predict, col=myCol)

#plot outliers on first 2 PCs
par(mfrow = c(2, 1))
for (i in 1:2)
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
#plot outlier distribution
plot(x, option = "qqplot")
#Manhattan plot
plot(x,option="manhattan")

#isolate number of outliers (P-value 0.05)
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)
#3018 outliers identified 
#record which SNPs are outliers
snp_pc <- get.pc(x, outliers)

#SAME analysis but with LD clumping, sliding window 500 and threshhold 0.2
res <- pcadapt(input=genos, K = 10, LD.clumping = list(size = 500, thr = 0.1), min.maf = 0.01)
plot(res, option = "screeplot")

#choose k = 2 *could also be 5* 
res <- pcadapt(input=genos, K = 2, LD.clumping = list(size = 500, thr = 0.1), min.maf = 0.01)
par(mfrow = c(2, 1))
for (i in 1:2)
  plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
plot(res)
plot(res, option = "qqplot", threshold = 0.1)
# Distribution of Mahalanobis distances.
plot(res, option = "stat.distribution")

#get outliers based on p < 0.01
qval <- qvalue(res$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qval < alpha)
length(outliers)

#315 outliers identified 
#record which SNPs are outliers
snp_pc_ld <- get.pc(res, outliers)
write.csv(snp_pc_ld, "outliers_ld.csv")

#pca by group (here, predicted population assignment from DAPC)
plot(res,option="scores",pop=pop$pop, col=myCol)

par(mfrow = c(1, 1))
#are outliers correlated to latitude along PC1
plot(res$scores[,1]~pop$lat,pch=19,col="mediumseagreen")
cor.test(res$scores[,1], pop$lat)

#are outliers correlated to longitude?
plot(res$scores[,1]~pop$lon,pch=19,col="mediumseagreen")
cor.test(res$scores[,1], pop$lon)

##########################################################
#Pop gen parameters
#snps into genind for pop stats

#create hierfstat object
dm <- genind2hierfstat(mouse_gen)

#genetic diversity 
div <- summary(mouse_gen)
div

#plot observed het
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
#plot observed vs expected het
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")

#W&C FST and FIS 
wc(mouse_gen)

#pairwise Fst
genet.dist(mouse_gen, method = "WC84")
#bootstrapped confidence intervals
boot.ppfst(dat=mouse_gen,nboot=1000,quant=c(0.025,0.975),diploid=TRUE)

#use hierfstat to get basic stats
basicstat <- basic.stats(mouse_gen, diploid = TRUE, digits = 1) 
names(basicstat)

boot.ppfis(dat=mouse_gen,nboot=1000,quant=c(0,1.0),diploid=TRUE)

#compiled as basic.stats: 
ho <- basicstat$Ho
write.csv(x=ho, file = "ho")
hs <- basicstat$Hs
write.csv(x=hs, file = "hs")
fis <- basicstat$Fis
write.csv(x=fis, file = "fis")

#allelic richness per locus and population
ar <- allelic.richness(mouse_gen,min.n=NULL,diploid=TRUE)
write.csv(x=ar, file = "allelic_richness")
