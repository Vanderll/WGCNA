####################################
# Get adjacency matrix for Spencer #
####################################

rm(list=ls())
setwd("/data2/kiemele/WGCNA/HXBbrain_forPhenogen/20141113")
library(WGCNA)
options(strainsAsFactors=FALSE)

##########################################
# Load in Dataset Used to Create Modules #
##########################################

load("heritExprs_strainMeans.Rdata")
#rows = transcript clusters, columns = strains

###############################
# Create a Correlation Matrix #
###############################

corMatrix = bicor(t(heritMeans))
dim(corMatrix)

###############################
# Create the Adjacency Matrix #
###############################

adjMatrix_all = corMatrix^9
dim(adjMatrix_all)

save(adjMatrix_all, file="adjMatrix_all.Rdata")

###################################################
# Create a List of Adjacency Matrices (by Module) #
###################################################
key = read.csv("moduleMembership.csv")

table(rownames(adjMatrix_all)==key[,"TC"])

TCwant = list()
adjMatricesByModule = list()
for(i in names(table(key[,"moduleColors"]))){	
	TCwant[[i]] = as.vector(key[which(key[,"moduleColors"]==i),"TC"])
	adjMatricesByModule[[i]] = adjMatrix_all[TCwant[[i]], TCwant[[i]]]
}

save(adjMatricesByModule, file="adjMatricesByModule.Rdata")

##Try creating this list so the names and numbers correspond to the module label;
length(names(table(key[,"moduleColors"])))

forAdj = names(table(key[,"moduleLabels"]))
forAdj = forAdj[-which(forAdj=="0")]

TCwant = list()
adjMatricesByModule = list()
for(i in forAdj){	
	TCwant[[i]] = as.vector(key[which(key[,"moduleLabels"]==i),"TC"])
	adjMatricesByModule[[i]] = adjMatrix_all[TCwant[[i]], TCwant[[i]]]
}


labelColorKey = key[order(key[,"moduleLabels"]), c(2:3)]

sig = c()
for(i in 2:nrow(labelColorKey)){
	sig = c(sig, labelColorKey[i,2]==labelColorKey[(i-1),2])
}

sig2 = c("FALSE", sig)

uniqueLabelColorKey = labelColorKey[-which(sig2=="TRUE"),]
uniqueLabelColorKey = uniqueLabelColorKey[-which(uniqueLabelColorKey[,"moduleColors"]=="grey"),]

table(uniqueLabelColorKey[,"moduleLabels"]==names(adjMatricesByModule))
names(adjMatricesByModule) = uniqueLabelColorKey[,"moduleColors"]

save(adjMatricesByModule, file="adjMatricesByModule.Rdata")
