setwd("/home/kiemele/Thesis/WGCNA/Robustness/GeneNetwork/VTA/bootstrapMethod")
#setwd("C:/Users/kiemelel/Documents/Thesis/WGCNA BXD RI/Laura Saba Batch Effects/Robustness/GeneNetworkCompare/VTA/bootstrapMethod")
rm(list=ls())

library(WGCNA)
options(stringsAsFactors = FALSE)

##Load Results 
load(file = "modulePreservation_bootstraps_group1.RData")
mp1 = mp
load(file = "modulePreservation_bootstraps_group2.RData")
mp2 = mp
load(file = "modulePreservation_bootstraps_group3.RData")
mp3 = mp
load(file = "modulePreservation_bootstraps_group4.RData")
mp4 = mp
load(file = "modulePreservation_bootstraps_group5.RData")
mp5 = mp
load(file = "modulePreservation_bootstraps_group6.RData")
mp6 = mp
load(file = "modulePreservation_bootstraps_group7.RData")
mp7 = mp
load(file = "modulePreservation_bootstraps_group8.RData")
mp8 = mp
load(file = "modulePreservation_bootstraps_group9.RData")
mp9 = mp
load(file = "modulePreservation_bootstraps_group10.RData")
mp10 = mp

##Load module assignments
assign = read.csv(file="vta_assign_labelNEW.csv", header=TRUE)
assign = as.matrix(assign[,"labelNEW"])
nModules = max(assign) + 2  #add 2 for the grey and test modules the modulePreservation calculates
nModules
nBoots = length(mp1)

##Make empty matrices for Zsummary & p-values
Zsummary_group1 = matrix(nrow=nModules, ncol = nBoots)
log_pval_group1 = matrix(nrow=nModules, ncol = nBoots)
Zsummary_group2 = matrix(nrow=nModules, ncol = nBoots)
log_pval_group2 = matrix(nrow=nModules, ncol = nBoots)
Zsummary_group3 = matrix(nrow=nModules, ncol = nBoots)
log_pval_group3 = matrix(nrow=nModules, ncol = nBoots)
Zsummary_group4 = matrix(nrow=nModules, ncol = nBoots)
log_pval_group4 = matrix(nrow=nModules, ncol = nBoots)
Zsummary_group5 = matrix(nrow=nModules, ncol = nBoots)
log_pval_group5 = matrix(nrow=nModules, ncol = nBoots)
Zsummary_group6 = matrix(nrow=nModules, ncol = nBoots)
log_pval_group6 = matrix(nrow=nModules, ncol = nBoots)
Zsummary_group7 = matrix(nrow=nModules, ncol = nBoots)
log_pval_group7 = matrix(nrow=nModules, ncol = nBoots)
Zsummary_group8 = matrix(nrow=nModules, ncol = nBoots)
log_pval_group8 = matrix(nrow=nModules, ncol = nBoots)
Zsummary_group9 = matrix(nrow=nModules, ncol = nBoots)
log_pval_group9 = matrix(nrow=nModules, ncol = nBoots)
Zsummary_group10 = matrix(nrow=nModules, ncol = nBoots)
log_pval_group10 = matrix(nrow=nModules, ncol = nBoots)

#######################################################
# Extract Just the preservation Z summary and p-value #
#######################################################
for(i in c(1:nBoots)){
	Zsummary_group1[,i] = mp1[[i]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
	log_pval_group1[,i] = mp1[[i]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
}
for(i in c(1:nBoots)){
	Zsummary_group2[,i] = mp2[[i]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
	log_pval_group2[,i] = mp2[[i]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
}
for(i in c(1:nBoots)){
	Zsummary_group3[,i] = mp3[[i]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
	log_pval_group3[,i] = mp3[[i]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
}
for(i in c(1:nBoots)){
	Zsummary_group4[,i] = mp4[[i]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
	log_pval_group4[,i] = mp4[[i]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
}
for(i in c(1:nBoots)){
	Zsummary_group5[,i] = mp5[[i]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
	log_pval_group5[,i] = mp5[[i]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
}
for(i in c(1:nBoots)){
	Zsummary_group6[,i] = mp6[[i]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
	log_pval_group6[,i] = mp6[[i]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
}
for(i in c(1:nBoots)){
	Zsummary_group7[,i] = mp7[[i]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
	log_pval_group7[,i] = mp7[[i]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
}
for(i in c(1:nBoots)){
	Zsummary_group8[,i] = mp8[[i]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
	log_pval_group8[,i] = mp8[[i]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
}
for(i in c(1:nBoots)){
	Zsummary_group9[,i] = mp9[[i]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
	log_pval_group9[,i] = mp9[[i]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
}
for(i in c(1:nBoots)){
	Zsummary_group10[,i] = mp10[[i]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
	log_pval_group10[,i] = mp10[[i]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
}

rownames(Zsummary_group1) = rownames(mp1[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)
rownames(Zsummary_group2) = rownames(mp2[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)
rownames(Zsummary_group3) = rownames(mp3[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)
rownames(Zsummary_group4) = rownames(mp4[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)
rownames(Zsummary_group5) = rownames(mp5[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)
rownames(Zsummary_group6) = rownames(mp6[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)
rownames(Zsummary_group7) = rownames(mp7[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)
rownames(Zsummary_group8) = rownames(mp8[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)
rownames(Zsummary_group9) = rownames(mp9[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)
rownames(Zsummary_group10) = rownames(mp10[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)

colnames(Zsummary_group1) = c(paste(rep("boot", each=nBoots), rep(c(1:nBoots), sep=""))) ##THIS NEED TO CHANGE FOR EACH GROUP!!!
colnames(Zsummary_group2) = c(paste(rep("boot", each=nBoots), rep(c(11:(10+nBoots)), sep=""))) ##THIS NEED TO CHANGE FOR EACH GROUP!!!
colnames(Zsummary_group3) = c(paste(rep("boot", each=nBoots), rep(c(21:(20+nBoots)), sep=""))) ##THIS NEED TO CHANGE FOR EACH GROUP!!!
colnames(Zsummary_group4) = c(paste(rep("boot", each=nBoots), rep(c(31:(30+nBoots)), sep=""))) ##THIS NEED TO CHANGE FOR EACH GROUP!!!
colnames(Zsummary_group5) = c(paste(rep("boot", each=nBoots), rep(c(41:(40+nBoots)), sep=""))) ##THIS NEED TO CHANGE FOR EACH GROUP!!!
colnames(Zsummary_group6) = c(paste(rep("boot", each=nBoots), rep(c(51:(50+nBoots)), sep=""))) ##THIS NEED TO CHANGE FOR EACH GROUP!!!
colnames(Zsummary_group7) = c(paste(rep("boot", each=nBoots), rep(c(61:(60+nBoots)), sep=""))) ##THIS NEED TO CHANGE FOR EACH GROUP!!!
colnames(Zsummary_group8) = c(paste(rep("boot", each=nBoots), rep(c(71:(70+nBoots)), sep=""))) ##THIS NEED TO CHANGE FOR EACH GROUP!!!
colnames(Zsummary_group9) = c(paste(rep("boot", each=nBoots), rep(c(81:(80+nBoots)), sep=""))) ##THIS NEED TO CHANGE FOR EACH GROUP!!!
colnames(Zsummary_group10) = c(paste(rep("boot", each=nBoots), rep(c(91:(90+nBoots)), sep=""))) ##THIS NEED TO CHANGE FOR EACH GROUP!!!

rownames(log_pval_group1) = rownames(Zsummary_group1)
colnames(log_pval_group1) = colnames(Zsummary_group1)
rownames(log_pval_group2) = rownames(Zsummary_group2)
colnames(log_pval_group2) = colnames(Zsummary_group2)
rownames(log_pval_group3) = rownames(Zsummary_group3)
colnames(log_pval_group3) = colnames(Zsummary_group3)
rownames(log_pval_group4) = rownames(Zsummary_group4)
colnames(log_pval_group4) = colnames(Zsummary_group4)
rownames(log_pval_group5) = rownames(Zsummary_group5)
colnames(log_pval_group5) = colnames(Zsummary_group5)
rownames(log_pval_group6) = rownames(Zsummary_group6)
colnames(log_pval_group6) = colnames(Zsummary_group6)
rownames(log_pval_group7) = rownames(Zsummary_group7)
colnames(log_pval_group7) = colnames(Zsummary_group7)
rownames(log_pval_group8) = rownames(Zsummary_group8)
colnames(log_pval_group8) = colnames(Zsummary_group8)
rownames(log_pval_group9) = rownames(Zsummary_group9)
colnames(log_pval_group9) = colnames(Zsummary_group9)
rownames(log_pval_group10) = rownames(Zsummary_group10)
colnames(log_pval_group10) = colnames(Zsummary_group10)

Zsummary_all = cbind(Zsummary_group1, Zsummary_group2, Zsummary_group3, Zsummary_group4, Zsummary_group5, Zsummary_group6,
	Zsummary_group7, Zsummary_group8, Zsummary_group9, Zsummary_group10)

log_pval_all = cbind(log_pval_group1, log_pval_group2, log_pval_group3, log_pval_group4, log_pval_group5, log_pval_group6,
	log_pval_group7, log_pval_group8, log_pval_group9, log_pval_group10)

pval_all = exp(log_pval_all)

save(Zsummary_all, pval_all, file="ZsummaryResults_AllInfo.Rdata")

#########
Zsummary_means = apply(Zsummary_all, 1, mean)
Zsummary_sd = apply(Zsummary_all, 1, sd)
Zsummary_4view = cbind(Zsummary_means, Zsummary_sd)

key = read.table(file="moduleLabelKey.txt", sep="\t", header=TRUE)
test = merge(key, Zsummary_4view, by.x = "label", by.y = 0)

write.csv(test, file="modulePreservation_ZsummaryResults.csv")

####################
###Original Code ###
####################
moduleSize = mp[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"moduleSize"]
Zsummary = mp[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"Zsummary.pres"]
log_pval = mp[[1]]$preservation$log.p$ref.VTA$inColumnsAlsoPresentIn.bootstrap[,"log.psummary.pres"]
pval = exp(log_pval)

test = cbind(moduleSize, Zsummary, pval)
rownames(test) = rownames(mp[[1]]$preservation$Z$ref.VTA$inColumnsAlsoPresentIn.bootstrap)

write.csv(test, file="Zsummary_boots.csv")
