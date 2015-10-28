library(foreach)
library(doMC)
registerDoMC()
library(proxy)
library(WGCNA)
library("org.Mm.eg.db")
library(biomaRt)
library(GOstats)
library("org.Mm.eg.db")
library("edgeR")


#source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
 #biocLite("GOstats")

getDoParWorkers()
options(cores=27)
getDoParWorkers()

#enableWGCNAThreads(nThreads = 27)

setwd("/home/dan/workDir/HDID2")

source("functionDefinitions.R")
try(dir.create("resultsCoexpr"), silent = F)
try( dir.create("figuresCoexpr"), silent = F)


load("data/selectedData.RData")
sample_info=read.csv("data/sampleInfo.csv", row.names=1)
colnames(selectedGeneCounts)= colnames(sample_info)
geneNames=rownames(selectedGeneCounts)
load("data/gene_biotype.RData")

# divide the data in different groups
HDID_M=selectedGeneCounts[,sample_info["group",]=="HDID2" & sample_info["gender",]=="M"]
HDID_F=selectedGeneCounts[,sample_info["group",]=="HDID2" & sample_info["gender",]=="F"]
HSNPT_M=selectedGeneCounts[,sample_info["group",]=="HSNpt" & sample_info["gender",]=="M"]
HSNPT_F=selectedGeneCounts[,sample_info["group",]=="HSNpt" & sample_info["gender",]=="F"]

countsMales=selectedGeneCounts[,sample_info["gender",]=="M"]
countsFemales=selectedGeneCounts[,sample_info["gender",]=="F"]

adjConsensusMales=adjacency(t(countsMales), power=1)
adjConsensusFemales=adjacency(t(countsMales), power=1)

save(adjConsensusMales, adjConsensusFemales,  file="data/adjConsensus.RData")

adj_HDID_M=adjacency(t(HDID_M), power=1, corFnc="bicor")
save(adj_HDID_M,   file="data/adj_HDID_M.RData")
rm(adj_HDID_M)

adj_HDID_F=adjacency(t(HDID_F), power=1, corFnc="bicor")
save(adj_HDID_F,   file="data/adj_HDID_F.RData")
rm(adj_HDID_F)

adj_HSNPT_M=adjacency(t(HSNPT_M), power=1, corFnc="bicor")
save(adj_HSNPT_M,   file="data/adj_HSNPT_M.RData")
rm(adj_HSNPT_M)

adj_HSNPT_F=adjacency(t(HSNPT_F), power=1, corFnc="bicor")
save(adj_HSNPT_F,   file="data/adj_HSNPT_F.RData")
rm(adj_HSNPT_F)

load("data/adjConsensus.RData")

########################################################################################################################################
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(adjConsensus, powerVector = powers, verbose = 5, moreNetworkConcepts=T)

plotNetConstruction(sft)
quartz.save("figuresCoexpr/netConstructionCoexpr.tif", type="tif", bg="white", dpi=300)
quartz.save("figuresCoexpr/netConstructionCoexpr.jpg", type="jpg", bg="white")

softPowerCoexpr=6
adjCoexpr=adjConsensus^softPowerCoexpr
adjCoexpr[is.na(adjCoexpr)]=0
hierADJCoexpr = hclust(as.dist(1-adjCoexpr),method="average");

# code below might need modifications to make the number of modules between ~ 15 and 50, and to make number of grey genes less than ~ 2-3000
hybridCoexpr=cutreeHybrid(dendro = hierADJCoexpr, distM=1-adjCoexpr, cutHeight = 0.9995, minClusterSize = 100, deepSplit = 4, maxCoreScatter = NULL, minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL, pamStage = TRUE, pamRespectsDendro = F, useMedoids = FALSE,  respectSmallClusters = TRUE, verbose = 2, indent = 0)

colorsCoexpr = labels2colors(hybridCoexpr$labels)
names(colorsCoexpr)=geneNames
table(colorsCoexpr)
length(table(colorsCoexpr))
modulesCoexpr=names(table(colorsCoexpr))
sum(colorsCoexpr=="grey")

# length(table(colorsCoexpr))
# [1] 32
# > modulesCoexpr=names(table(colorsCoexpr))
# > sum(colorsCoexpr=="grey")
# [1] 53

fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")

writeLines(paste("Number modules ", length(table(colorsCoexpr)), sep=','), fileConnSummary)
writeLines(paste("Number grey genes  ",  sum(colorsCoexpr=="grey"), sep=','), fileConnSummary)

close(fileConnSummary)


save(adjCoexpr, colorsCoexpr, file="data/adjCoexprModules.RData")
load("data/adjCoexprModules.RData")

########################################################################################################
# save the results below for use with enrinchR
coexprConn=intramodularConnectivity(adjCoexpr, colorsCoexpr, scaleByMax=T)
totalScaledConnectivity=coexprConn[,"kTotal"]/max(coexprConn[,"kTotal"])

coexprResultsTable=cbind(colorsCoexpr, round(coexprConn[,"kWithin"],3), round(totalScaledConnectivity,3))
colnames(coexprResultsTable)=c("module", "moduleScaledConn", "networkScaledConn")

#write.csv(coexprConn, file=paste("resultsCoexpr/coexprModulesConnInfo",".csv", sep=""))  

for (module in modulesCoexpr){
  print(module)
    
  currModuleInfo=cbind(rownames(coexprConn)[colorsCoexpr==module],round(coexprConn[colorsCoexpr==module,"kWithin"],2))
  write.csv(currModuleInfo, file=paste("resultsCoexpr/module_", module, ".csv", sep=""), row.names=F, col.names=F)  

  fileConn<-file(paste("resultsCoexpr/moduleConn_", module, ".txt", sep=""),  open="at")
      
  for (i in 1:dim(currModuleInfo)[1]){
    writeLines(paste(currModuleInfo[i,1],currModuleInfo[i,2],sep=','), fileConn)
  }
  close(fileConn)
  

}




#############################################################################
# GO annotations

load("data/transcriptInfoMouse.RData")
annotateMouseModulesGO(colorsCoexpr, transcriptInfoMouse, type="Coexpr")

##############################################################################
# evaluate differential expression
geneReadsRaw=read.csv("data/RNASeq016_mm10_gene_reads_not_normalized.csv", header=F, row.names=1)
geneNamesRaw=rownames(geneReadsRaw)

sample_info=read.csv("data/sampleInfo.csv", row.names=1)
colnames(geneReadsRaw)= colnames(sample_info)

# divide the data in different groups
HDID_M=geneReadsRaw[,sample_info["group",]=="HDID2" & sample_info["gender",]=="M"]
HDID_F=geneReadsRaw[,sample_info["group",]=="HDID2" & sample_info["gender",]=="F"]
HSNPT_M=geneReadsRaw[,sample_info["group",]=="HSNpt" & sample_info["gender",]=="M"]
HSNPT_F=geneReadsRaw[,sample_info["group",]=="HSNpt" & sample_info["gender",]=="F"]

groupSelectionM=c(rep("HSNPT_M",dim(HSNPT_M)[2]),rep("HDID_M",dim(HDID_M)[2]))
groupSelectionM =factor(groupSelectionM)

#deMales=differentialExpression(cbind(HSNPT_M, HDID_M), groupSelectionM)
d=DGEList(counts= cbind(HSNPT_M, HDID_M), group= groupSelectionM)
d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d, dispersion="tagwise") 
de.calls <- decideTestsDGE(de.tgw, p=0.05, adjust="BH")

resultsDEMales=cbind(gene_biotype[rownames(de.tgw$table)], de.tgw$table, de.calls)
write.csv(resultsDEMales, file="results/DEmales.csv")


groupSelectionF=c(rep("HSNPT_F",dim(HSNPT_F)[2]),rep("HDID_F",dim(HDID_F)[2]))
groupSelectionF =factor(groupSelectionF)

#deFemales=differentialExpression(cbind(HSNPT_F, HDID_F), groupSelectionF)
d=DGEList(counts= cbind(HSNPT_F, HDID_F), group= groupSelectionF)
d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d, dispersion="tagwise") 
de.calls <- decideTestsDGE(de.tgw, p=0.05, adjust="BH")

resultsDEFemales=cbind(gene_biotype[rownames(de.tgw$table)],de.tgw$table, de.calls)
#write.csv(resultsDEFemales, file="results/DEFemales.csv")

selectedResDEF=resultsDEFemales[rownames(coexprResultsTable),]
selectedResDEM=resultsDEMales[rownames(coexprResultsTable),]

connResultsFemales=cbind(coexprResultsTable, selectedResDEF)
connResultsMales=cbind(coexprResultsTable, selectedResDEM)

write.csv(connResultsFemales, file="resultsCoexpr/ConnectivityResultsFemales.csv")
write.csv(connResultsMales, file="resultsCoexpr/ConnectivityResultsMales.csv")


# load("data/gene_biotype.RData")
# 
# gene_biotype_net=gene_biotype[names(colorsCoexpr)]
# connResultsFemales=cbind(gene_biotype_net, connResultsFemales)
# connResultsMales=cbind(gene_biotype_net, connResultsMales)
# 
# write.csv(connResultsFemales, file=paste("resultsCoexpr/coexprModulesFemales",".csv", sep=""))
# write.csv(connResultsMales, file=paste("resultsCoexpr/coexprModulesMales",".csv", sep=""))


##############################################################################
# evaluate module enrichment in differentially expressed genes
# these modules are next utilized as inputs for enrichR

pValuesFemales=resultsDEFemales[ ,"PValue"]
names(pValuesFemales)=rownames(resultsDEFemales)
femaleDEModuleEnrich=moduleEnrichment(colorsCoexpr, pValuesFemales)

#femaleDEModuleEnrich=moduleEnrichment(colorsCoexpr, pValuesFemales)
# Bonferroni corrected module enrichments
modulesCoexpr[femaleDEModuleEnrich[1,]<0.05]
#[1] "brown"       "cyan"        "pink"        "red"         "royalblue"   "saddlebrown" "tan"         "yellow"     
fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
writeLines(paste("Female Modules enriched in DE genes ", modulesCoexpr[femaleDEModuleEnrich[1,]<0.05], sep=','), fileConnSummary)
close(fileConnSummary)

pValuesMales=resultsDEMales[ ,"PValue"]
names(pValuesMales)=rownames(resultsDEMales)

maleDEModuleEnrich=moduleEnrichment(colorsCoexpr, pValuesMales)
modulesCoexpr[maleDEModuleEnrich[1,]<0.05]

fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")

writeLines(paste("Male Modules enriched in DE genes ", modulesCoexpr[maleDEModuleEnrich[1,]<0.05], sep=','), fileConnSummary)

close(fileConnSummary)


#black"         "cyan"          "darkgreen"     "darkgrey"      "darkturquoise" "green"         "royalblue"     "salmon"        "steelblue"    
#[10] "yellow"   
##############################################################################
# evaluate module einrichment in non coding RNA - do not run for polyA datasets


##########################################################################################################
biotypeCategories=names(table(gene_biotype_net))

biotypeModuleTable = matrix(0, nrow = length(biotypeCategories), ncol = length(modulesCoexpr));
rownames(biotypeModuleTable)= biotypeCategories
colnames(biotypeModuleTable)= modulesCoexpr

logPtable = biotypeModuleTable 
textTable= matrix(data="", nrow = length(biotypeCategories), ncol = length(modulesCoexpr))

for (module in modulesCoexpr){
  print(module)
  moduleGenes=geneNames[which(colorsCoexpr==module)]
  for (biotype in biotypeCategories){
    print(biotype)
    biotypeGenes=geneNames[which(gene_biotype_net==biotype)]
    logPtable[biotype, module]=-log10(fisher.test(geneNames %in% moduleGenes, geneNames %in% biotypeGenes, alt="g")$p.value)+runif(1,0,0.2)
    
  }
}

logPtable[logPtable > 10]=10
logPtable[logPtable < 3]=0

for (module in modulesCoexpr){
  
  for (biotype in biotypeCategories){
    
    logPtable[biotype, module]=logPtable[biotype, module]+runif(1,0.7,0.9)
    
  }
}
#################################################################

par(mar=c(10,10,4,4))
labeledHeatmap(logPtable, xLabels = biotypeCategories, yLabels = modulesCoexpr, 
               colors = greenWhiteRed(100)[50:100],
               setStdMargins = FALSE, 
               textMatrix = textTable,cex.text=0.8,
               main = paste("Module enrichment in ncRNA "));

quartz.save("figuresCoexpr/NCRNA_modules.tif", type="tif", bg="white", dpi=300)
quartz.save("figuresCoexpr/NCRNA_modules.jpg", type="jpg")

#############################################################################

#diffConnM = evaluateDiffConn(colorsCoexpr, HSNPT_F, HDID_F, softPower=6, nPerm=4)
  
