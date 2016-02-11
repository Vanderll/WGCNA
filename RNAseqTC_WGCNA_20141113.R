########################
#   Set up Workspace   #
########################
setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141113")
rm(list=ls())
library(WGCNA)
options(strainsAsFactors=FALSE)

#################
#   Load Data   #
#################
load("heritExprs_strainMeans.Rdata")
ls()
dim(heritMeans)
heritMeans[1:5, 1:5]

datExpr = t(heritMeans)
dim(datExpr)

#########################################
# Check out the Soft Thresholding Power #
#########################################

# Choose a set of soft-thresholding powers
powers = c(c(1:15), seq(from = 16, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power

pdf(file="SoftThresholdingPower.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

save(sft, file="SoftThresholdingPower.Rdata")


##############################################################
##  WGCNA in LXS Brain Expression Data (RNAseq Derived TC)  ##
##############################################################
dim(datExpr)
net = blockwiseModules(datExpr,power=9, minModuleSize=5, deepSplit=4, numericLabels=TRUE, 
	pamRespectsDendro=FALSE,saveTOMs=TRUE, saveTOMFileBase="brain",verbose=3,networkType="unsigned", corType="bicor")

table(net$colors)

pdf(file="WGCNAdendrograms.pdf")

mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
plotDendroAndColors(net$dendrograms[[2]],mergedColors[net$blockGenes[[2]]], "Module colors", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
plotDendroAndColors(net$dendrograms[[3]],mergedColors[net$blockGenes[[3]]], "Module colors", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
TC = colnames(datExpr)

moduleMembership = cbind(TC, moduleLabels, moduleColors)
length(unique(moduleMembership[,"moduleLabels"]))
#length = 930 (includes the grey module so 929 modules);

write.csv(moduleMembership, file="moduleMembership.csv",row.names=FALSE)

MEs = net$MEs
geneTree1 = net$dendrograms[[1]]
geneTree2 = net$dendrograms[[2]]
geneTree3 = net$dendrograms[[3]]

save(net,moduleLabels,moduleColors,MEs,file="LXSbrain.RNAseqTC.WGCNAcontruct.20141112.Rdata")




objects()
length(moduleColors)
length(MEs)
length(moduleLabels)
dim(datExpr)
probeorder <- colnames(datExpr)
length(probeorder)
write.csv(probeorder, file="ProbeOrder.csv")
write.csv(moduleColors, file="moduleColors.csv")
write.csv(moduleLabels, file="moduleLabels.csv")

summary()
sizes = as.numeric(table(moduleMembership[-which(moduleMembership[,"moduleLabels"]=="0"),"moduleColors"]))
summary(sizes)

MEs_noGray = MEs[,-which(colnames(MEs)=="ME0")]
write.csv(MEs_noGray, file="MEs.csv")

#the gray module (tc unassigned to a module) is how big:
dim(moduleMembership[which(moduleMembership[,"moduleLabels"]=="0"),]))