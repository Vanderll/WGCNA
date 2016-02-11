#Set up workspace;
rm(list=ls())
setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141113")
library(WGCNA)
options(stringsAsFactors=FALSE)

#####################################################################
# Step1: filter on probesets present in at least 5% of the samples  #
#####################################################################
#load in DABG values
dabg = read.table(file="/data2/saba/ForPhenoGen/HXB.BXH.Brain.fullPS/dabg.fullPS.HXB_BXH.brain.PhenoGen.txt",sep="\t",header=TRUE)
#load in expression data
expr = read.table(file="/data2/saba/ForPhenoGen/HXB.BXH.Brain.fullPS/rma.fullPS.HXB_BXH.brain.PhenoGen.txt",sep="\t",header=TRUE)

table(dim(dabg)==dim(expr))

presentDABG <- dabg[rowSums(dabg<0.0001)>(ncol(expr)*0.05),]
presentExpr <- expr[rowSums(dabg<0.0001)>(ncol(expr)*0.05),]

#176,230 present probesets

table(nrow(presentDABG)==nrow(presentExpr))
table(ncol(presentDABG)==ncol(presentExpr))

save(presentDABG, presentExpr, file="presentDABGandExpr.Rdata")

######################################################
# Step2: find those present probesets within an exon #
######################################################
#load in file with probeset to exon/transcriptID info;
load("presentDABGandExpr.Rdata")
exonInfo = read.table(file="HXB.Brain.PolyA.psToRNAGene.txt",sep="\t", header=FALSE)
colnames(exonInfo) = c("probeset_id", "transcript_cluster")

length(unique(exonInfo[,"probeset_id"]))
length(unique(exonInfo[,"transcript_cluster"]))
#11,273 unique transcript clusters


presentExonExpr = merge(exonInfo, presentExpr, by.x="probeset_id", by.y="probeset_id")
rownames(presentExonExpr) = presentExonExpr[,"probeset_id"]
length(unique(presentExonExpr[,"transcript_cluster"]))
#8,715 unique genes

presentExonExpr = presentExonExpr[,-c(1:2)]

presentKey = merge(exonInfo, as.matrix(presentExpr[,1]), by.x="probeset_id", by.y=1)
#############################################
# Step 3: Look for outlier samples (arrays) #
#############################################

sampleTree = flashClust(as.dist(1-cor(presentExonExpr)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

##Looks like we have 3 sample outliers (height > 0.03):
toRemove = c("PD_8604_02_Brain.CEL")

sig=c()
for(i in toRemove){
	sig = c(sig, which(colnames(presentExonExpr)==i))
}

length(sig)

PresentExonExpr_clean = presentExonExpr[,-sig] 

####################################################################################
# Step 4: Correlate strain means & determine what PS should be collapsed together  #
####################################################################################
strain = sapply(strsplit(colnames(PresentExonExpr_clean), split="_", fixed=TRUE), "[[", 1)
means = t(apply(PresentExonExpr_clean, 1, function(a) lm(a~as.factor(strain) -1)$coefficients))
colnames(means) = sapply(strsplit(colnames(means), split=")", fixed=TRUE), "[[", 2)

#remove parents and the other inbreds (PD, WKY.Lj) from correlation analysis
toDelete = c("BN.LX", "PD", "SHR.H", "SHR.lj", "SHR.Lx", "WKY.Lj")
sig = c()
for(i in toDelete){
	sig = c(sig, which(colnames(means)==i))
}

length(sig)
RImeans = means[,-sig]


#make a list of matrices, where each matrix represents the expression matrix for a specific TC
TCmatrices = function(x, key){
	xList = list()
	sig = list()
	TC = as.matrix(unique(key[,"transcript_cluster"]))
	for(i in TC){
		sig[[i]] = key[which(key[,"transcript_cluster"]==i),"probeset_id"]
		xList[[i]] = merge(x, as.matrix(sig[[i]]), by.x=0, by.y=1)
		rownames(xList[[i]]) = xList[[i]]$Row.names
		xList[[i]] = xList[[i]][,-1]
	}
	return(xList)
}

TCmeansMatrices = TCmatrices(RImeans, presentKey)

#now we want to know what ps can be collapsed down together
#then we can make a function and get the 1st principal component when for those with the individual array data;

cluster = function(x){
  if(nrow(x)>1){	
  tree = hclust(as.dist(1-cor(t(x))), method = "average")
  cutTree = cutree(tree, h=0.75)
   }
 if(nrow(x)==1){
  cutTree = 1
  names(cutTree) = rownames(x)
  }
 return(cutTree)		
}

sampleTree = hclust(as.dist(1-cor(t(TCmeansMatrices[[1]]))), method = "average")
sizeGrWindow(12,9)
pdf(file = "exampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Clustering Example: gene XLOC_022227", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 0.75, col="red")
abline(h = 0.5, col="blue")

dev.off()

#5 if height = 0.75
#9 if height = 0.5

cluster(TCmeansMatrices[[1]])


clustersToMake = lapply(TCmeansMatrices, cluster)

getClustNames = function(x){
	names = list()
	PS = list()
	for(i in 1:length(x)){
	  names[[i]] = as.matrix(paste(names(x)[[i]], ".clust", x[[i]], sep=""))
 	  PS[[i]] = names(x[[i]])	
	}
	collapsedResults = cbind(unlist(PS), unlist(names)) 
	colnames(collapsedResults) = c("probeset", "transcriptClusterCollapsed")
	return(collapsedResults)
}

PSclustKey = getClustNames(clustersToMake)
dim(PSclustKey)

length(unique(PSclustKey[,"transcriptClusterCollapsed"]))

save(PSclustKey, clustersToMake, TCmeansMatrices, file="findingCollapsedTC.Rdata")
write.csv(PSclustKey, file="PStoCollapsedTransciptClusterKey.csv")

################################################################################################
# Step 5: Generate the 1st principal component scores for those that need to be collapsed down #
################################################################################################
#want to use the PresentExonExpr_clean... But need to remove the parentals;

strain = sapply(strsplit(colnames(PresentExonExpr_clean), split="_", fixed=TRUE), "[[", 1)
#remove parents and the other inbreds (PD, WKY.Lj) when generating the PCA analyses;
toDelete = c("BN.LX", "PD", "SHR.H", "SHR.lj", "SHR.Lx", "WKY.Lj")
sig = c()
for(i in toDelete){
	sig = c(sig, which(strain==i))
}

length(sig)
RIforCollapse = PresentExonExpr_clean[,-sig]

TCcollapsed = collapseRows(RIforCollapse, PSclustKey[,"transcriptClusterCollapsed"], PSclustKey[,"probeset"], method="ME", thresholdCombine=NA)
length(TCcollapsed)
save(TCcollapsed, file="TCcollapsed.datExpr.Rdata")

##################################################################################################
# Step 6: Calculate the heritability and apply a filter of > 0.25 to stay in dataset for network #
##################################################################################################
expr = TCcollapsed$datETcollapsed
strain = sapply(strsplit(colnames(expr), split="_", fixed=TRUE), "[[", 1)

test = apply(expr, 1, function(a) summary(lm(a~as.factor(strain)))$r.squared)

herits = cbind(rownames(expr), test)
colnames(herits) = c("transcript", "heritability")
write.table(herits, file="herits.clustersTry20141113.HXB.brain.txt", sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

heritable <- herits[herits[,2]>0.25,1]
heritablePS <- match(heritable,rownames(expr))
heritablePS <- heritablePS[!is.na(heritablePS)]
heritExprs <- expr[heritablePS,]
save(heritExprs, file="heritExprs.Rdata")

pdf(file = "heritabilities.pdf", width = 12, height = 9);
hist(as.numeric(herits[,2]), breaks=100, main = "Histogram of Heritability", xlim=c(0,1), axes=TRUE)
abline(v = 0.25, col="red")
dev.off()
#51,654 after heritanility filter;

heritMeans = t(apply(heritExprs, 1, function(a) lm(a~as.factor(strain) -1)$coefficients))
colnames(heritMeans) = sapply(strsplit(colnames(heritMeans), split=")", fixed=TRUE), "[[", 2)
save(heritMeans, file="heritExprs_strainMeans.Rdata")

####################################
# Step 7: Look for outlier strains #
####################################
sampleTree = flashClust(as.dist(1-cor(heritMeans)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
pdf(file = "strainClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Strain clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
#keep all strains;


