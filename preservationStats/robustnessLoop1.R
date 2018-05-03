setwd("/home/kiemele/Thesis/WGCNA/Robustness/GeneNetwork/VTA/bootstrapMethod")

#setwd("C:/Users/kiemelel/Documents/Thesis/WGCNA BXD RI/Laura Saba Batch Effects/Robustness/GeneNetworkCompare/VTA/bootstrapMethod")
getwd()

#######Set Up the R session#####
rm(list=ls())
# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
set.seed(2020)

#######Data Input###########
 
load("bootstrap_strains.Rdata")
#this loads in object boots_all which is a matrix where each row is a bootstrap pull and each column is a bootstrap


load("/home/kiemele/Thesis/WGCNA/Robustness/GeneNetwork/VTA/VTA_4robust.Rdata")
VTA_labels = as.matrix(read.csv(file = "vta_assign_labelNEW.csv", header=TRUE))
#this is a module membership matrix where each row is a probeset and then there are columns: "Label" and "Color" which correspond to the module assigned to and another column called "labelNew" is for the modules we want to test (i.e. candidate modules).  This reduces computation time when you test less modules.  So for this case non-candidate modules are labeled 0 and then all else labeled 1-5 depending on the module (1 = aquamarine4, 2 = snow2, 3 = tomato2, 4=mediumturquoise and 5 = indianred3)

VTA_labels = as.matrix(VTA_labels[,"labelNEW"])


##DUE TO SLOW COMPUTATION TIME, DO 10 BOOTSTRAPS AT A TIME... NEED TO FIGURE OUT HOW TO SPEED UP 
boots_all = boots_all[,c(1:10)]

nBoots = ncol(boots_all)
nStrain = nrow(boots_all)
nProbesets = ncol(VTA_4robust)

boot_data = vector("list",nBoots)
mp = vector("list", nBoots)

nBoots=2
for(i in c(1:nBoots)){

boot_data[[i]] = matrix(nrow=nStrain, ncol = nProbesets)
boot_data[[i]] = VTA_4robust[boots_all[,i],]
rownames(boot_data[[i]]) = c(1:nStrain)

###Calculate Module Preservation Stats###
setLabels = c("VTA", "bootstrap");
multiExpr = list(VTA = list(data=VTA_4robust), bootstrap = list(data=boot_data[[i]]));
multiColor = list(VTA = VTA_labels);
set.seed(2020)

system.time( {
mp[[i]] = modulePreservation(multiExpr, multiColor,
referenceNetworks = 1,
nPermutations = 200, 
randomSeed = 1,
maxModuleSize = 300,
maxGoldModuleSize = 300,
quickCor = 0,
verbose = 3)
} );
# Save the results
save(mp, file = "modulePreservation_bootstraps_group1.RData");
}
