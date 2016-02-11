##Get a ME file with the Strains as rownames;

setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141113")
rm(list=ls())
library(WGCNA)
options(strainsAsFactors=FALSE)

#load in MEs;
load("LXSbrain.RNAseqTC.WGCNAcontruct.20141112.Rdata")

#load in original data used for WGCNA;
load("heritExprs_strainMeans.Rdata")
write.csv(colnames(heritMeans), file="origStrainToStrainKey.csv")
#in excel, put the new strain names next to the old ones;
#load strain key back into R;

key = read.csv(file="origStrainToStrainKey.csv")

#get new ME file;
#attach the original strain names (becuase it is based on this order);
rownames(MEs) = key[,"origStrain"]

MEs_new = merge(key, MEs, by.x="origStrain", by.y=0)
rownames(MEs_new) = MEs_new[,"strain"]
MEs_new = MEs_new[order(MEs_new[,"strain"]),]

#just keep the rownames as an identifier;
MEs_new  = MEs_new[,-c(1:2)]

MEs_noGray = MEs_new[,-which(colnames(MEs_new)=="ME0")]


write.csv(MEs_noGray, file="MEs_noGray.csv")
