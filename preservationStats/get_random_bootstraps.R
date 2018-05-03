
######################
#  Set Up Work Space #
######################
rm(list=ls())
set.seed(2020)
setwd("C:/Users/kiemelel/Documents/Thesis/WGCNA BXD RI/Laura Saba Batch Effects/Robustness/GeneNetworkCompare/VTA/bootstrapMethod")

#########################
#  Make List of Strains #
#########################
load(file="C:/Users/kiemelel/Documents/Thesis/WGCNA BXD RI/Laura Saba Batch Effects/Robustness/GeneNetworkCompare/VTA/VTA_4robust.Rdata")
#each row is a strain and each column is a probeset

VTA_strains <- as.matrix(rownames(VTA_4robust))
length(VTA_strains)
dim(VTA_4robust)

#################################
# Get 100 different subsets of  #
# Strains for Bootstrapping     #
#################################
nBoots = 100
nStrains = length(VTA_strains)
nProbesets = ncol(VTA_4robust)


resamples <- lapply(1:nBoots, function(i)
sample(VTA_strains, size=nStrains, replace = T))

boots_all = do.call(cbind, resamples)
dim(boots_all)

bootNames = c(paste(rep("boot", each=nBoots), rep(c(1:nBoots), sep="")))
pullNames = c(paste(rep("pull", each=nStrains), rep(c(1:nStrains), sep="")))

colnames(boots_all) = bootNames
rownames(boots_all) = pullNames

save(boots_all, file="bootstrap_strains.Rdata")

####################################################
####################################################

