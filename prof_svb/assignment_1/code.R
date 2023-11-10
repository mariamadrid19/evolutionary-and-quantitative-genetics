library(adegenet)
setwd("/Users/mariamadrid/Documents/Evo-Quant Genetics/svb-exercises")

library(vcfR)
data <- read.vcfR("tssb.vcf")
data <- vcfR2genind(data)
data_pop <- read.table("tssb_pop.txt", header=F)
data@pop <- as.factor(data_pop$V2)
data
str(data)
data@pop
data@ploidy
popNames(data)
Hexp <- Hs(data)
n.pop <- seppop(data)
Hobs <- do.call("c",lapply(n.pop,function(x)mean(summary(x)$Hobs)))
Hobs[is.nan(Hobs)] <- NA
Expected_H <- barplot(Hexp,xlab="Population", ylab="Hexp")
Observed_H <- barplot(Hobs,xlab="Population", ylab="Hobs")
H_combined <- rbind(Hexp,Hobs)
barplot(H_combined,beside=T,xlab="Population",col=c("lightblue","blue"),legend=T,ylim = c(0,0.3)) 

library(hierfstat)
data2 <- genind2hierfstat(data)
pp <- genet.dist(data2,diploid=T,method="WC84")
pp
write.table(as.matrix(pp),file="pairwiseWCfst.txt")

dapcstickle <- dapc(data,n.pca=40,n.da=5)
dapcstickle
scatter(dapcstickle)

clust <- find.clusters(data,n.pca=100)
as.data.frame(clust$grp)
table(clust$grp,data$pop)
table.value(table(clust$grp,data$pop),col.labels = popNames(data))

library(corrplot)
#load the data
save(pp, file = "pairwiseWCfst.rda")

#prepare a heat map of the pair- wise FST values
corrplot(as.matrix(pp), is.cor = F, type= 'lower', method = "color", col=colorRampPalette(c("skyblue2","royalblue4"))(200))

#create a data frame that contains the hierarchical structure of the data
dat_hier = data.frame(c(rep('inland', 9), rep('coastal', 4), 
                        rep('inland', 6), rep('coastal', 48), 
                        rep('inland', 32), rep('coastal', 20), 
                        rep('inland', 17), rep('coastal', 18), 
                        rep('inland', 20), rep('coastal', 1)), data@pop)

#give a name to the columns
names(dat_hier) = c('group', 'pop')

#add the data frame to data and make it a stratification factor
data@other = list(population_hierarchy = dat_hier) 
strata(data) = data.frame(other(data)$population_hierarchy)
strata(data)

#make sure not to have too many missing data
#markers with 5% or more missing values would be dropped in the AMOVA
library(poppr) 
info_table(data, plot = TRUE)

#remove missing data
dat2 = missingno(data, type = 'genotype')

#run the AMOVA
stickamova = poppr.amova(dat2, ~group/pop, within = F)
stickamova

#run a permutation test and see where results feature in the output
set.seed(1999)
sticksignif = randtest(stickamova, nrepet = 999)
sticksignif

#plot the permutation test results 
library(mdatools) 
plot(sticksignif)

#install the R packages ‘fsthet’, ‘outflank’ and ‘pcadapt’ to detect selection in SNP data

#Method 1 = FSTHET
library(fsthet)
#Read the file tssb.gen, which contains the reduced data set of Raeymaekers et al. (2017)
gpop <- my.read.genepop("tssb(1).gen", ncode = 2L, quiet = F)

#calculate the actual FST between the markers
fsts<-calc.actual.fst(gpop, fst.choice = "fst")

#plot the FST on the Ht: what do you see? 
plot(fsts$Ht, fsts$Fst,xlab="Ht",ylab="Fst",pch=19)

#run 100 sampling of the loci. This is a procedure to get the empirical confidence intervals for 
#outlier FST values (specified in this function by bootstrap=FALSE). We set the number of loci in each 
#quantile bin to be 10 and we run each sampling for all 4 FST options
quant.out1 = as.data.frame(t(replicate(100, fst.boot(gpop,bootstrap=FALSE, fst.choice="betahat", min.per.bin=10))))
fsts_theta <- calc.actual.fst(gpop, fst.choice = "theta")
quant.out2 = as.data.frame(t(replicate(100, fst.boot(gpop,bootstrap=FALSE, fst.choice="theta", min.per.bin=10))))
quant.out3 = as.data.frame(t(replicate(100, fst.boot(gpop,bootstrap=FALSE, fst.choice="var", min.per.bin=10))))
quant.out4 = as.data.frame(t(replicate(100, fst.boot(gpop,bootstrap=FALSE, fst.choice="fst", min.per.bin=10))))

# look for outliers
outliers1<-find.outliers(fsts,boot.out=quant.out1)
outliers2<-find.outliers(fsts_theta,boot.out=quant.out2)
outliers3<-find.outliers(fsts,boot.out=quant.out3)
outliers4<-find.outliers(fsts,boot.out=quant.out4)
outliers1 
outliers2 
outliers3 
outliers4



#colour graph output with points
plot(fsts_theta$Ht, fsts_theta$Fst,xlab="Ht",ylab="Fst",pch=19,col="black")
points(outliers2$Ht,outliers2$Fst,xlab="Ht",ylab="Fst",pch=19,col="red")

plot(fsts$Ht, fsts$Fst,xlab="Ht",ylab="Fst",pch=19,col="black")
points(outliers1$Ht,outliers1$Fst,xlab="Ht",ylab="Fst",pch=19,col="pink")
points(outliers2$Ht,outliers2$Fst,xlab="Ht",ylab="Fst",pch=19,col="red")
points(outliers3$Ht,outliers3$Fst,xlab="Ht",ylab="Fst",pch=19,col="darkblue") 
points(outliers4$Ht,outliers4$Fst,xlab="Ht",ylab="Fst",pch=19,col="lightgreen")

#THETA
#fsts_theta <- calc.actual.fst(gpop, fst.choice = "theta")
#quant.out2 = as.data.frame(t(replicate(100, fst.boot(gpop,bootstrap=FALSE, fst.choice="theta", min.per.bin=10))))
#outliers2 <- find.outliers(fsts_theta,boot.out=quant.out2)
plot(fsts_theta$Ht, fsts_theta$Fst,xlab="Ht",ylab="Fst",pch=19,col="black")
points(outliers2$Ht,outliers2$Fst,xlab="Ht",ylab="Fst",pch=19,col="pink")
points(outliers2$Ht,outliers2$Fst,xlab="Ht",ylab="Fst",pch=19,col="red")

#Method 2 = OutFLANK
#Install OutFLANK
library(reshape2)
library(qvalue)
library(vcfR)
library(plyr)
library(processx)
library(OutFLANK)
source("read-plink-bed.R")

#read in the SNP, marker, and population data
snps = read.plink('pl1')# SNPS
markers = read.table("pl1.bim", h = F)# Markers
pops = read.table("tssb_pop.txt", h = F)# populations
pops <- pops[-c(2, 3, 5, 8, 9, 12, 13, 14, 23, 25, 26, 28, 34, 35, 54, 55, 61, 68, 69, 70, 71, 72, 73, 74, 75, 76, 78, 79, 84, 86, 87, 88, 90, 92, 93, 94, 104, 114, 127, 128, 131, 141, 143, 148, 150, 151, 152, 159, 163, 168, 170, 171, 172), ] # remove individuals with too much missing values
#replace the missing data coding from NA to 9
snps[is.na(snps)] = 9

#Remove low quality SNPs with lots of missing data
t2 = apply(snps, 2, function(x){sum(x == 9)/length(x)})#we are calculating what percentage of the genotypes in each SNP is missing
snps2 = snps[, t2 <= .25]# we are pruning the data so that we keep SNPs missing in no more than 25% of the data

#run OutFLANK
ofl = MakeDiploidFSTMat(snps2, as.character(markers[t2 <= .25,2]), as.character(pops[,2]))# we are making a Fst matrix
outR = OutFLANK(ofl, NumberOfSamples = 4)
#look at the OutFLANK results
OutFLANKResultsPlotter(outR)
OutFLANKResultsPlotter(outR, Hmin = 0.01, Zoom = T)
head(outR$results, n = 10)

#sort the loci by FST value
outR$results[order(outR$results$FST, decreasing=TRUE),]

#Plot the FST on the Ht and highlight the outliers
plot(outR$results$He, outR$results$FST, pch=20, col="grey") 
points(outR$results$He[outR$results$qvalues<0.1],
       outR$results$FST[outR$results$qvalues<0.1], pch=21, col="red")


#Method 3 = PCAdapt
#Install PCAdapt
library(pcadapt)

#Read the BED file
pcadapt <- read.pcadapt("pl1.bed", type = "bed")

#Run a PCA with a large number of principal components
pca <- pcadapt(pcadapt, K = 20)

#Choose the best number of principal components (K) based on plots
plot(pca, option = "screeplot") # screeplot = eigenvalues of each PCs

#Create a population file
poplist.names <-c(rep("L11",2), rep("U01",2),rep("L02",2), rep("L10",5), rep("L01",3), rep("L05",1), 
                  rep("L06",2), rep("L05",17), rep("L06",16), rep("L12",1), rep("L11",12), 
                  rep("L02",18), rep("L10",14), rep("L01",12), rep("U01",14), rep("L06",1)) 
print(poplist.names)

#Create a score plot that diplays population structure (based on PC1 and PC2)
plot(pca, option = "scores", pop = poplist.names)

#Redo the PCA with the chosen value of K (number of clusters)
pca2 <- pcadapt(pcadapt, K = 3)

#Create a Manhattan plot to look at the SNPs with significant p-values
plot(pca2, option = "manhattan")

#Create a Q-Q plot to look at the difference between expected and observed p- values 
#in order to confirm the presence of outliers
plot(pca2, option = "qqplot")

#Create a histogram of the p-values
hist(pca2$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

#try 3 different methods to choose the cutoff for the outlier detection

#Using the q-value
library(qvalue)
qval <- qvalue(pca2$pvalues)$qvalues 
alpha <- 0.1
outliers <- which(qval < alpha) 
length(outliers)

#Using the Benjamini-Hochberg procedure
padj <- p.adjust(pca2$pvalues, method = "BH") 
alpha <- 0.1
outliers <- which(padj < alpha) 
length(outliers)

#Using the Bonferroni correction
padj <- p.adjust(pca2$pvalues, method = "bonferroni") 
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)

#To extract sign outliers from PCAdapt
sign_pca <- markers[outliers,] 
pca_outliers_results <- cbind(sign_pca, padj[outliers])
pca_outliers_results <- as.data.frame(pca_outliers_results)
str(pca_outliers_results)
install.packages("xlsx")
library(xlsx)
write.xlsx(pca_outliers_results, file = "/Users/mariamadrid/Documents/Evo-Quant Genetics/svb-exercises/pca_outliers_results.xlsx")

#write outliers3 (fsthet) to excel since it's the largest file 
outliers3 <- as.data.frame(outliers3)
write.xlsx(outliers3, file = "/Users/mariamadrid/Documents/Evo-Quant Genetics/svb-exercises/outliers3.xlsx")
