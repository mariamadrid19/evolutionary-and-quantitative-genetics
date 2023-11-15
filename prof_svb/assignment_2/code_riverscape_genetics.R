## THE DATA
## --------

# data from:

# Raeymaekers JAM, Maes GE, Geldof S, Hontis I, Nackaerts K & Volckaert FAM (2008). 
# Modeling genetic connectivity in sticklebacks as a guideline for river restoration. 
# Evolutionary Applications 1, 475-488.
# Available for download on http://onlinelibrary.wiley.com/doi/10.1111/j.1752-4571.2008.00019.x/abstract

browseURL("http://www.researchgate.net/publication/229560553_Modeling_genetic_connectivity_in_sticklebacks_as_a_guide_for_river_restoration")

# THIS PUBLICATION CONTAINS THE SAME TABLES AND FIGURES AS THE ONE THIS SCRIPT GENERATES. SO, WHEN YOU SWITCH BACK
# AND FORTH BETWEEN THE SCRIPT AND THE PUBLICATION, YOU CAN COMPARE YOUR RESULTS DIRECTLY.
# TO HELP YOU WITH THIS, NOTE THAT THE SCRIPT REFERS TO THE CORRESPONDING FIGURE OR TABLE.

# THE EXERCISE HAS TWO PARTS:
# PART I: IMPORTING AND FIRST LOOK AT THE GENETIC DATA (AS BEFORE, BUT NOW STARTING FROM AN EXCEL FILE)
# PART II: LINKING GENETIC DATA WITH GEOGRAPHIC DATA (THE NEW STUFF)
# IT IS NOT NEEDED TO COMPLETE PART I BEFORE STARTING WITH PART II

# PART I: importing and first look at the genetic data ####

# Download the data (Raeymaekers et al-EVA-2008.xlsx) (Click on "View" and "Download" at the right)

browseURL("http://www.researchgate.net/publication/256520299_Raeymaekers_et_al-EVA-2008") 

# set your working directory to the destination where you downloaded the file, for instance:

setwd("/Users/mariamadrid/Documents/Evo-Quant Genetics/svb-assignment-2")

# import the data

library(readxl)

genotypes0 <- read_excel("Raeymaekersetal-EVA-2008.xlsx", sheet = "genotypes Demer 2002") # Specify sheet with a number or name
genotypes0 # inspect the data. 
genotypes0 <- as.data.frame(genotypes0)

# Note that some genotypes have 5 digits and some have 6 digit format! This is because excel stores genotypes
# such as 096096 (i.e. categorically) as 96096 (i.e. numerically). Furthermore, note that missing values are marked as 0
# So, we need to convert all data back to the 6-digit format. This can be done with the converse function below:

colnames(genotypes0)
loci <- colnames(genotypes0)[3:8];loci # define the loci
genotypes <- genotypes0 # make a copy of the genotypes object

converse <- function(x) {
  x0 <- as.vector(x)
  x0[x %in% c(1:99999)] <- paste0(0,x0[x %in% c(1:99999)],sep="")
  x0[x %in% c(0)] <- "000000"
  x0
}

genotypes[,loci] <- apply(genotypes[,loci],2,converse)

View(genotypes0) # old data
View(genotypes) # converted data. Was it successful?

# Now we can convert the data also to a genind object so we can analyse the data with adegenet:

library(adegenet)
stickle <- df2genind(genotypes[,loci],sep=NULL,ncode=3,NA.char="000",ind.names=genotypes$ID,pop=genotypes$Population)

summary(stickle)
indNames(stickle)
dim(stickle@tab)
alleles(stickle)

# So, remember that it is always possible to correctly import genotype data with R for use with adegenet, 
# even starting from a random data format. We can now also convert it to the hierfstat format, 
# for instance to calculate allelic richness:

library(devtools)
install_github("jgx65/hierfstat")
library("hierfstat")
?hierfstat

stickle3 <- hierfstat::genind2hierfstat(stickle)
AR <- allelic.richness(stickle3)
meanAR <- colMeans(AR$Ar)
meanAR

# Finally, let's quantify pairwise fst for these 21 populations:

matFst <- pairwise.WCfst(stickle)
matFst 
matFst.ordered <- matFst[c(6:21,1:5),c(6:21,1:5)] # same, but arranged in a more logical order (downstream to upstream)
matFst.ordered

# Note: you can also go back to the script "3-Conservation Genetics - stickle.R" and obtain the entire VARFINAL object
# with "stickle" and "stickle3" as the input (see +/- line 330 onwards) - for those that want to repeat the first class. 
# It gives you detailed information for genetic diversity across the 21 populations, similar to Table 1 in the paper.
# Likewise, you can test for deviations of HWE, and calculate global Fst, make a dapc plot, etcetera.


# PART II: linking genetic data with geographic data ####

# The objects meanAR and matFst.ordered you have obtained in PART I contain the genetic information we need for the rest of the 
# exercise. They correspond to the 2nd and 3rd sheet of the excel file you have downloaded.
# However, the excel sheets also contain geographical information which we want to use to see which geographic features drive
# genetic variation (allelic richness) and genetic differentiation (pairwise fst). 

geogenvar0 <- read_excel("Raeymaekersetal-EVA-2008.xlsx", sheet = "geo & genvar")
geogendist0 <- read_excel("Raeymaekersetal-EVA-2008.xlsx", sheet = "geo & gendist")

geogenvar <- data.frame(geogenvar0)
geogendist <- data.frame(geogendist0)

# Compare geogenvar with VARFINAL and geogendist with matFst.ordered. We continue with geogenvar and geogendist, since they
# allow us to reproduce the results of the source publication (fst in geogendist object is estimated differently than in matFst.ordered)

# 1) genetic variability and geographical features:

geogenvar
colnames(geogenvar)
geogenvar[,c("population","AR")] # allelic richness (AR) for each population

# 2) genetic differentiation and pairwise geographical features:

geogendist
colnames(geogendist)

# Note that each row in geogendist corresponds to a population pair: S4a-S5a, S4a-S5b, ..., S4a-S14, S5a-S5bd, etc...

# pairwise genetic differentiation (fst) for each population pair:

geogendist[,"fst"] 

# Note that with 21 (n) populations there are 210 possible population pairs (n*(n-1)/2):

length(geogendist[,"fst"])
length(geogendist$fst)

# An easier way to present each of the variables in the geogendist matrix
# can be obtained after running this function:

column.to.matrix <- function(a){nr.obj <- round((1+(1+8*length(a))^0.5)/2);m <- matrix(0,nr.obj,nr.obj);m[lower.tri(m)] <- a;m + t(m)}

# this can be applied now and the result is the following symmetrical matrix (distance matrix):

column.to.matrix(geogendist$fst)

# This is a 21 by 21 population matrix, note the zeros on the diagonal.
# We will sometimes use the column format, and sometimes the matrix format.

# We are now ready to start building landscape genetic models. This means we will choose a genetic parameter Y (allelic richness or
# pairwise fst), and a set of geographical variables X1, X2, X3, .... (for instance, geographic distance and number of barriers)

## PACKAGES TO DOWNLOAD FOR THIS SESSION:

# websites:
# http://www.r-project.org/
# http://www.freestatistics.org/cran/
# http://probability.ca/cran/

# call packages:

library(MASS)
library(vegan)
library(Hmisc)

## CORRELATIONS AND PLOTS (compare the results with the results in the paper)

# It is crucial to verify how environmental variables (X1, X2, X3) are correlated before building landscape genetic models for Y.
# We start with the geogenvar data (Y = allelic richness; X1, X2, X3 = geographical features)

# As a general rule, first visualise the relationships graphically before performing statistical tests:

pairs(geogenvar[,c(6,9,10)]) # This plots all pairwise combinations for the three X variables we are planning to use

# Next, lets have a look at the correlation matrix between all X values.

cor(geogenvar[,c(6,9,10)]) # This shows the correlation matrix between the X variables (Pearson correlation)
cor.test(geogenvar$log10.upstream.distance,geogenvar$log10.habitat.width) # This function also calculates the correlation, and tests if it is significantly different from zero
rcorr(as.matrix(geogenvar[,c(6,9,10)])) # This function gives the correlation matrix (upper output) as well as the P-values (lower output)

# !!! Check if the results of the rcorr funcation are the same as in Table 4A in the publication.
# Note that if you use rcorr (library Hmisc) then you don't need cor and cor.test (library stats)
# What do you conclude based on the plots and the statistics? Do you think we have a problem of multicollinearity here? 

# ANSWER: multicolinearity is the problem that when a set of explanatory variables (such as X1 and X2 and X2) are too similar
# in their relationship with a dependent variable (Y), and hence will "compete" strongly when trying to explain Y. Especially
# when the correlation between X1, X2 and X3 is very strong, multicollinearity can be an issue. In this case, all
# the correlations are significant, but the correlations are not extremely high. The pairs() function shows that 
# the relationships show quite some scatter. This suggests we will not have big multicollinearity issues here.

# However, we want a statistical confirmation of this,  one that also takes into account how X1, X2 and X3 co-vary with Y. 
# This can be done with the "variance inflation factor".
# Search in the stickleback paper for "variance inflation factor" (vif) and find out how it has been used.

library(car)
lm1 <- lm(AR~all_barriers + log10.habitat.width. + log10.upstream.distance., data=geogenvar)
vif(lm1) # no values above 10 for this model

# No vif higher than 10 is observed, so according to this rule of thumb, all X variables can stay in the model
# without multicolinearity issues. 

# The next step is to compare the relationship between Y (allelic richness) and all X variables.
# For instance, to generate figure 2A,C,E,G of the paper, run the following code:

par(mfrow = c(2,2))
plot(geogenvar$distance,geogenvar$AR, xlab="River distance (km)", ylab="Allelic richness")
abline(lm(AR ~ distance, data = geogenvar), col = "blue")
fit_distance <- lm(AR ~ distance, data = geogenvar)
legend("bottomleft",legend=paste("R2=", format(summary(fit_distance)$r.squared,digits=3)))
plot(geogenvar$all_barriers,geogenvar$AR, xlab="All barriers", ylab="Allelic richness")
abline(lm(AR ~ all_barriers, data = geogenvar), col = "blue")
fit_all_barriers <- lm(AR ~ all_barriers, data = geogenvar)
legend("bottomleft",legend=paste("R2=", format(summary(fit_all_barriers)$r.squared,digits=3)))
plot(geogenvar$log10.habitat.width.,geogenvar$AR, xlab="Log(Habitat width)", ylab="Allelic richness")
abline(lm(AR ~ log10.habitat.width., data = geogenvar), col = "blue")
fit_habitat <- lm(AR ~ log10.habitat.width., data = geogenvar)
legend("topleft",legend=paste("R2=", format(summary(fit_habitat)$r.squared,digits=3)))
plot(geogenvar$log10.upstream.distance.,geogenvar$AR, xlab="Log(Upstream distance)", ylab="Allelic richness")
abline(lm(AR ~ log10.upstream.distance., data = geogenvar), col = "blue")
fit_upstream <- lm(AR ~ log10.upstream.distance., data = geogenvar)
legend("topleft",legend=paste("R2=", format(summary(fit_upstream)$r.squared,digits=3)))

# Carefully inspect these plots. What do you conclude?
# ANSWER: at first sight, barriers seem to explain allelic richness better than any other variable.
# Further down in the script we will run a statistical model to see if this can be confirmed when comparing several linear models

# Now let's run all of the steps above again, but now for the geogendist data (Y = pairwise Fst; X = pairwise geographic features)
# So, first let's check the correlation between the pairwise geographical features:

pairs(as.matrix(geogendist[,c(7,9,10,11)])) # Always first make graphs! They tell you much more than statistics.
cor(geogendist[,c(7,9,10,11)])
rcorr(as.matrix(geogendist[,c(7,9,10,11)]))

# !!! Check if the results of the rcorr funcation are the same as in Table 4B in the publication.
# Did you obtain here the same R and P-values as in the paper? 

# ANSWER: you should obtain the same R values, but different P-values. This is because
# the P-values in table 4B are based on the so-called Mantel test.

library(vegan)
?mantel # Find out that there is a simple and a partial mantel test

# As an example of the use of the simple and partial mantel tests, let's consider the relationship between Y (pairwise fst) and X1 or X2 (any pairwise geographical feature)

# simple mantel tests:

mantel(column.to.matrix(geogendist$distance), column.to.matrix(geogendist$fst), method="pearson", permutations=999)
mantel(column.to.matrix(geogendist$all_barriers), column.to.matrix(geogendist$fst), method="pearson", permutations=999)

# The simple mantel test is testing the correlation between Y (pairwise fst) and a single X (either distance or barriers)
# We can conclude from the simple mantel test that fst is correlated with both distance and barriers. 

# partial mantel tests:
mantel.partial(column.to.matrix(geogendist$fst), 
               column.to.matrix(geogendist$distance), 
               column.to.matrix(geogendist$all_barriers), method = "pearson", permutations = 999)
mantel.partial(column.to.matrix(geogendist$fst), 
               column.to.matrix(geogendist$all_barriers), 
               column.to.matrix(geogendist$distance), method = "pearson", permutations = 999)

# The partial mantel test is testing the correlation between Y (first entry in the function) 
# and X1 (second entry in the function) after correcting for X2 (third entry in the function). 
# Such correction is important for us, since we want to know in our case that 
# the increase in genetic divergence with barriers is not just due to distance.
# We conclude that only the effect of barriers on fst is substantial, since this effect remains significant even when we account for distance (second mantel.partial test)
# The effect of distance on fst is not substantial, since this effect is no longer significant when we account for barriers (first mantel.partial test)

# But why do we need the non-parametric Mantel correlation test rather than a Pearson correlation test?

mantel(column.to.matrix(geogendist$fst), column.to.matrix(geogendist$distance), method="pearson") # Mantel correlation: correct R-values and correct P-values
cor.test(geogendist$fst, geogendist$distance, method="pearson") # Pearson correlation: correct R-values but wrong P-values

# ANSWER: because the Pearson correlation assumes independent data points, and pairwise distances
# are not independent. So, while the Pearson correlation and the Mantel correlation are identical, 
# the P-value associated with the Pearson correlation (rcorr function) is not valid because
# we would be violating a basic statistical assumption. The Mantel test provides a valid 
# alternative since it is a permutation test that only relies on the data and not on underlying parametric theories.

# Independent data points: the age of student A, B and C
# Dependent data points: the distance between Antwerp, Ghent and Leuven
# Any pairwise measure such as pairwise fst and the pairwise number of barriers between sites are dependent data points.

# Note that to run the mantel tests, the data have to be converted from column format to a distance matrix.
# So, the Mantel test is expecting indeed a data structure compatible to the idea of dependent data points.

# The next step is to compare the relationship between Y (pairwise fst) and all X variables.
# For instance, to generate figure 2B,D,F,H of the paper, run the following code:

par(mfrow = c(2,2))
plot(geogendist$distance,geogendist$fst, xlab="River distance (km)", ylab="Fst")
abline(lm(fst ~ distance, data = geogendist), col = "blue")
fit_fst_distance <- lm(fst ~ distance, data = geogendist)
legend("topleft",legend=paste("R2=", format(summary(fit_fst_distance)$r.squared,digits=3)))
plot(geogendist$all_barriers,geogendist$fst, xlab="All barriers", ylab="Fst")
abline(lm(fst ~ all_barriers, data = geogendist), col = "blue")
fit_fst_barriers <- lm(fst ~ all_barriers, data = geogendist)
legend("topleft",legend=paste("R2=", format(summary(fit_fst_barriers)$r.squared,digits=3)))
plot(geogendist$log10.habitat.width.,geogendist$fst, xlab="Log(Habitat width)", ylab="Fst")
abline(lm(fst ~ log10.habitat.width., data = geogendist), col = "blue")
fit_fst_width <- lm(fst ~ log10.habitat.width., data = geogendist)
legend("topright",legend=paste("R2=", format(summary(fit_fst_width)$r.squared,digits=3)))
plot(geogendist$log10.upstream.distance.,geogendist$fst, xlab="Log(Upstream distance)", ylab="Fst")
abline(lm(fst ~ log10.upstream.distance., data = geogendist), col = "blue")
fit_fst_upstream <- lm(fst ~ log10.upstream.distance., data = geogendist)
legend("topright",legend=paste("R2=", format(summary(fit_fst_upstream)$r.squared,digits=3)))

# The plots confirm the results of the simple and partial mantel tests: the relationship between fst and barriers is stronger than between fst and distance.
# Further down in the script we will run a statistical model to see if this can be confirmed comparing several linear models.

# Now, compare these plots very carefully. Pay attention to the distribution of 
# the data points, including the outliers. What differences do you observe?
# Which evolutionary forces are acting here?

# ANSWER: isolation-by-geography plots might reflect the combined effect of gene flow (decreasing with geographical isolation, hence influencing the slope) and 
# genetic drift (inducing stochasticity which increases with geographical isolation, hence inducing heteroscedasticity)
# Both effects are stronger for the isolation-by-barrier plot than for the isolation-by-distance plot, suggesting that
# barriers, and not distance, control the balance between gene flow (~migration) and drift.

## NON-METRIC MULTIDIMENSIONAL SCALING (compare the results with the results in the paper)

# Before we do more statistical tests, we can first use multidimensional scaling on pairwise fst to visualise genetic differentiation
# This graphical method is an alternative to dapc.

# Figure 3A:
x <- isoMDS(column.to.matrix(geogendist$fst),y = cmdscale(column.to.matrix(geogendist$fst), 2), k = 2)
rownames(x$points) <- geogenvar$population
x

par(mfrow=c(1,1))
ordiplot(x,type="p") # you can ignore the warning here
text(x$points[,1],x$points[,2],dimnames(x$points)[[1]], cex=1)
title("Non-metric MDS plot of pairwise Fst among 21 populations")

## SPATIAL MODELS FOR ALLELIC RICHNESS

# As the final step, we make landscape genetic models for allelic richness (this section) and pairwise fst (next section)

# For allelic richness, We compare the following models:

lm1 <- lm(AR~all_barriers + log10.habitat.width. + log10.upstream.distance., data=geogenvar)
lm2 <- lm(AR~all_barriers + log10.habitat.width., data=geogenvar)
lm3 <- lm(AR~all_barriers + log10.upstream.distance., data=geogenvar)
lm4 <- lm(AR~log10.habitat.width. + log10.upstream.distance., data=geogenvar)
lm5 <- lm(AR~all_barriers, data=geogenvar)
lm6 <- lm(AR~log10.habitat.width., data=geogenvar)
lm7 <- lm(AR~log10.upstream.distance., data=geogenvar)

# we compare the models above using AIC, the Akaike's Information Criterium.
# AIC is calculated from least-squares regressions as
# AIC = 2K+ n ln(RSS/n), where K is the number of parameters, n
# is the number of populations, and RSS is the residual sum of squares.

# Actually here it was better to use AICc, the Akaike's Information Criterium corrected for small sample size:
# AICc= AIC + 2K(K+1)/(n-K-1)

# AICc for the models above:

2*4 + 21*log(anova(lm1)[4,2]/21) + 2*4*(4+1)/(21-4-1)
2*3 + 21*log(anova(lm2)[3,2]/21) + 2*3*(3+1)/(21-3-1)
2*3 + 21*log(anova(lm3)[3,2]/21) + 2*3*(3+1)/(21-3-1)
2*3 + 21*log(anova(lm4)[3,2]/21) + 2*3*(3+1)/(21-3-1)
2*2 + 21*log(anova(lm5)[2,2]/21) + 2*2*(2+1)/(21-2-1)
2*2 + 21*log(anova(lm6)[2,2]/21) + 2*2*(2+1)/(21-2-1)
2*2 + 21*log(anova(lm7)[2,2]/21) + 2*2*(2+1)/(21-2-1)

# conclusion: model lm3 is best because it has lowest AICc, closely followed by model 5.
# Hence, the smaller AICc, the better!
# The rule of thumb for AIC or AICc says that models with AIC values smaller than 2 are equivalent.

# The results for the "full model", shown in Table 3A of the publication, can be generated as follows:

library(car)
Anova(lm1)

# The P-values of this model indicate that barriers are the only significant predictor of allelic richness!


## MODELS FOR PAIRWISE FST

# As explained, here we work with pairwise distance matrices, but we can follow exactly the same reasoning as for allellic richness:

lm1 <- lm(fst~all_barriers + distance + log10.habitat.width. + log10.upstream.distance., data=geogendist)

lm2 <- lm(fst~all_barriers + distance + log10.habitat.width., data=geogendist)
lm3 <- lm(fst~all_barriers + distance + log10.upstream.distance., data=geogendist)
lm4 <- lm(fst~all_barriers + log10.habitat.width. + log10.upstream.distance., data=geogendist)
lm5 <- lm(fst~distance + log10.habitat.width. + log10.upstream.distance., data=geogendist)

lm6 <- lm(fst~all_barriers + distance, data=geogendist)
lm7 <- lm(fst~all_barriers + log10.habitat.width., data=geogendist)
lm8 <- lm(fst~all_barriers + log10.upstream.distance., data=geogendist)
lm9 <- lm(fst~distance + log10.habitat.width., data=geogendist)
lm10 <- lm(fst~distance + log10.upstream.distance., data=geogendist)
lm11 <- lm(fst~log10.habitat.width. + log10.upstream.distance., data=geogendist)

lm12 <- lm(fst~all_barriers, data=geogendist)
lm13 <- lm(fst~distance, data=geogendist)
lm14 <- lm(fst~log10.upstream.distance., data=geogendist)
lm15 <- lm(fst~log10.habitat.width., data=geogendist)

#AICc= AIC + 2K(K+1)/(n-K-1)

2*5 + 21*log(anova(lm1)[5,2]/21) + 2*5*(5+1)/(21-5-1)
2*4 + 21*log(anova(lm2)[4,2]/21) + 2*4*(4+1)/(21-4-1)
2*4 + 21*log(anova(lm3)[4,2]/21) + 2*4*(4+1)/(21-4-1)
2*4 + 21*log(anova(lm4)[4,2]/21) + 2*4*(4+1)/(21-4-1)
2*4 + 21*log(anova(lm5)[4,2]/21) + 2*4*(4+1)/(21-4-1)
2*3 + 21*log(anova(lm6)[3,2]/21) + 2*3*(3+1)/(21-3-1)
2*3 + 21*log(anova(lm7)[3,2]/21) + 2*3*(3+1)/(21-3-1)
2*3 + 21*log(anova(lm8)[3,2]/21) + 2*3*(3+1)/(21-3-1)
2*3 + 21*log(anova(lm9)[3,2]/21) + 2*3*(3+1)/(21-3-1)
2*3 + 21*log(anova(lm10)[3,2]/21) + 2*3*(3+1)/(21-3-1)
2*3 + 21*log(anova(lm11)[3,2]/21) + 2*3*(3+1)/(21-3-1)
2*2 + 21*log(anova(lm12)[2,2]/21) + 2*2*(2+1)/(21-2-1)
2*2 + 21*log(anova(lm13)[2,2]/21) + 2*2*(2+1)/(21-2-1)
2*2 + 21*log(anova(lm14)[2,2]/21) + 2*2*(2+1)/(21-2-1)
2*2 + 21*log(anova(lm15)[2,2]/21) + 2*2*(2+1)/(21-2-1)

# model 12 (including only barriers!!!) is best, and model 7 is fine as well.
# So barriers alone are sufficient to explain variation in pairwise fst!

# The results for the "full model", shown in Table 3B of the publication, can be APPROXIMATELY generated as follows
# (P-values differ from table 3B, since the method used is not completely identical):

library(ecodist)
?MRM
LM1 <- MRM(fst ~ all_barriers + log10.habitat.width. + log10.upstream.distance.+ distance , data=geogendist,nperm=10000)
LM1

# The P-values of this model indicate that barriers are the only significant predictor of pairwise fst!
# MRM is equivalent to a classical multiple regression, but uses a permutation test
# to generate P-values (since classical P-values are not valid). So, it is a sort extension of the Mantel test.


# Note that we can now also PREDICT pairwise Fst based on geographical features only:

predicted.fst <- predict(lm1)

# Here we use lm1 rather than lm12, because although barriers alone sufficiently explains fst, adding the three
# other geographical features will allow us to make more precise predictions. We are now no longer interested
# in finding the variables with the strongest effects, but in makeing good predictions.

# Some students questioned the use of "pairwise habitat width" and "pairwise upstream distance" in these models,
# because they felt it is kind of meaningless information. For instance, two populations with an intermediate
# habitat width will have an intermediate average habitat width, while a population from a small habitat width
# and a large habitat width would also have an average habitat width. Still, we can argue that these average 
# geographical features are  indicative for the average position of the population pair in the river network,
# and we know from figure 2F and 2G that there is a relationship between these properties and pairwise fst. So,
# including them in the model may lead to better predictions.

# With the following code, we can now make Figure 3A and B. It is a visualisation of OBSERVED pairwise fst
# (figure 3A) and PREDICTED pairwise fst (figure 3B). Note that both figures are indeed to some extent
# similar (for instance the position of the four satellite populations), 
# but figure A is ONLY based on genotypes, and figure B is ONLY based on geographical information.

par(mfrow=c(1,2))
x <- isoMDS(column.to.matrix(geogendist$fst),y = cmdscale(column.to.matrix(geogendist$fst), 2), k = 2)
rownames(x$points) <- geogenvar$population
ordiplot(x,type="p")
title("Observed pairwise Fst among 21 populations")
text(x$points[,1],x$points[,2],dimnames(x$points)[[1]], cex=1)

x <- isoMDS(column.to.matrix(predicted.fst),y = cmdscale(column.to.matrix(predicted.fst), 2), k = 2)
rownames(x$points) <- geogenvar$population
ordiplot(x,type="p")
title("Predicted pairwise Fst among 21 populations")
text(x$points[,1],x$points[,2],dimnames(x$points)[[1]], cex=1)

# For the correlation coefficient in figure 3B, use the following script:

cor(geogendist$fst,predicted.fst)^2 

# This R-squared value reflects the proportion of variation in pairwise fst which we can explain by geography!
# R-squared = 0.54
# This is relatively high - don't forget that patterns in natural populations (e.g. ecology) are usually harder to predict.

# The proportion we cannot account for is due to 1) genetic drift (=coincidence) and 
# 2) additional geographical features not accounted for in the model (e.g. barrier height or type).

# More details on the analysis can be found in the publication!



# A related study can be found here:

browseURL("http://www.researchgate.net/publication/227993944_Guidelines_for_restoring_connectivity_around_water_mills_A_population_genetic_approach_to_the_management_of_riverine_fish")

# How does the sampling design differs from the previous study? What is the advantage/disadvantage of this design?

# -- 
# Joost Raeymaekers
# 
# Ghent University
# Biology Department
# K. L. Ledeganckstraat 35, B-9000 Ghent, Belgium
# Web: http://www.ecology.ugent.be/terec/personal.php?pers=jr
# ResearchGate: http://www.researchgate.net/profile/Joost_Raeymaekers

# University of Leuven
# Laboratory of Biodiversity and Evolutionary Genomics
# Ch. de Beriotstraat 32, B-3000 Leuven, Belgium
# Phone: +32 (0) 16 37 36 49
