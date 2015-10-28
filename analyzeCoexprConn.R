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
disableWGCNAThreads()

registerDoMC()

getDoParWorkers()
options(cores=27)
getDoParWorkers()


setwd("workDir/HDID2")

source("functionDefinitions.R")
# dir.create("resultsCoexpr")
# dir.create("figuresCoexpr")


load("data/selectedData.RData")

#save(selectedGeneCounts, file="data/selectedGeneCounts.RData")

sample_info=read.csv("data/sampleInfo.csv", row.names=1)
colnames(selectedGeneCounts)= colnames(sample_info)
geneNames=rownames(selectedGeneCounts)

# divide the data in different groups
HDID_M=selectedGeneCounts[,sample_info["group",]=="HDID2" & sample_info["gender",]=="M"]
HDID_F=selectedGeneCounts[,sample_info["group",]=="HDID2" & sample_info["gender",]=="F"]
HSNPT_M=selectedGeneCounts[,sample_info["group",]=="HSNpt" & sample_info["gender",]=="M"]
HSNPT_F=selectedGeneCounts[,sample_info["group",]=="HSNpt" & sample_info["gender",]=="F"]


load("data/adjCoexprModules.RData")

########################################################################################################

#############################################################################
start.time <- Sys.time()

diffConnF = evaluateDiffConn(colorsCoexpr, HSNPT_F, HDID_F, softPower=6, nPerm=1000, nCores=10)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



diffConnM = evaluateDiffConn(colorsCoexpr, HSNPT_M, HDID_M, softPower=6, nPerm=1000, nCores=10)

save(diffConnM, diffConnF, file="data/diffConn.RData")
##############################################################################################
start.time <- Sys.time()

diffConnNPT = evaluateDiffConn(colorsCoexpr, HSNPT_F, HSNPT_M, softPower=6, nPerm=1000, nCores=10)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



diffConnHDID = evaluateDiffConn(colorsCoexpr, HDID_F, HDID_M, softPower=6, nPerm=1000, nCores=10)

save(diffConnNPT, diffConnHDID, file="data/diffConnNPTvsHDID.RData")

