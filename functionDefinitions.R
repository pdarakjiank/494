plotNetConstruction = function (sft)

{
sizeGrWindow(15, 9)
par(mfrow = c(2,3));
par(mar=c(5, 5, 3, 3))
par(lab=c(5, 7, 3))

cex1 = 1.5;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Fit",
     main = paste("Scale independence"), font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");

x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y.at <- c(seq(from = 0, to=1, by=.1))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)


abline(h=0.75,col="black")

mtext("A", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="l",
     main = paste("Mean Connectivity"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");


x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,5]))

y_step=(y_range[2]-y_range[1])/7

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("B", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Density as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,8],
     xlab="Soft Threshold (power)",ylab="Density",type="l",
     main = paste("Density"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");


x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,8]))

y_step=(y_range[2]-y_range[1])/5

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("C", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Centralization as a function of the soft-thresholding power
y_range=c(0,max(sft$fitIndices[,9]))

plot(sft$fitIndices[,1], sft$fitIndices[,9],
     xlab="Soft Threshold (power)",ylab="Centralization",type="l", ylim=y_range,
     main = paste("Centralization"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");


x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,9]))

y_step=(y_range[2]-y_range[1])/5

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("D", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Heterogeneity as a function of the soft-thresholding power
y_range=c(0,max(sft$fitIndices[,10]))

plot(sft$fitIndices[,1], sft$fitIndices[,10],
     xlab="Soft Threshold (power)",ylab="Centralization",type="l", ylim=y_range,
     main = paste("Heterogeneity"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");

x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,10]))

y_step=(y_range[2]-y_range[1])/5

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("E", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)

}

####################################################################################

annotateMouseModulesGO = function (colorsCoexpr, transcriptInfoMouse, type)
  # this assumes the gene identifiers are gene symbols, whichr are the names of the colorsCoexpr
  
{
  library("GOstats")
  
  fileName=paste("results", type,"/modules", type, "GO.csv", sep="")
  modules=names(table(colorsCoexpr))
  geneNames=names(colorsCoexpr)    
  networkTranscriptInfo= transcriptInfoMouse[which(transcriptInfoMouse[,"external_gene_name"] %in% geneNames),]
  
  pValue=0.01/length(modules)
  univ= unique(as.character(networkTranscriptInfo[,"entrezgene"]))
  
  univ=unique(univ)
  
  for (module in modules){
    print(module)
    geneNamesModule=geneNames[colorsCoexpr==module]
    
    moduleEntrez = unique(as.character(networkTranscriptInfo[which(networkTranscriptInfo[,"external_gene_name"] %in% geneNamesModule),"entrezgene"]))
    BPparams <- new("GOHyperGParams", geneIds = moduleEntrez, universeGeneIds = univ,  annotation = "org.Mm.eg.db",ontology = "BP", pvalueCutoff = pValue, conditional = T, testDirection = "over")
        
    BPres <- hyperGTest(BPparams)
    res=summary(BPres)
    print(res)
    res=cbind(rep(paste(module,"module BP"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
    
    
    CCparams <- new("GOHyperGParams", geneIds = moduleEntrez,
                    universeGeneIds = univ, annotation = "org.Mm.eg.db", 
                    ontology = "CC", pvalueCutoff = pValue, conditional = F,
                    testDirection = "over")
    
    CCres <- hyperGTest(CCparams)
    res=summary(CCres)
    res=cbind(rep.int(paste(module,"module CC"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
    
    
    MFparams <- new("GOHyperGParams", geneIds = moduleEntrez,
                    universeGeneIds = univ, annotation = "org.Mm.eg.db", 
                    ontology = "MF", pvalueCutoff = pValue, conditional = F,
                    testDirection = "over")
    
    MFres <- hyperGTest(MFparams)
    res=summary(MFres)
    res=cbind(rep.int(paste(module,"module MF"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
       
  }
    
}

##########################################################################################
differentialExpression = function (geneReads, groupFactor)
  
{
  d=DGEList(counts= geneReads, group= groupFactor)
  
  d <- estimateTagwiseDisp(d)
  de.tgw <- exactTest(d, dispersion="tagwise") 

  de.calls <- decideTestsDGE(de.tgw, p=0.05, adjust="BH")


  return(as.list(c(de.tgw, de.calls)))
}

##########################################################################################
#cbind(HDID_M, HSNPT_M),exonGeneNameSelected, geneNames, phenotypeVectorM, nPerm=0
# exonReads=cbind(HDID_M, HSNPT_M)
# exonGeneNames=exonGeneNameSelected
# phenotypeVector=phenotypeVectorM
# nPerm=0

splicingSignificance = function (exonReads, exonGeneNames, geneNames, phenotypeVector, nPerm)
  
{
  mantelGSS=rep(0,length(geneNames))
  names(mantelGSS)=geneNames
  mantelGSS_P= mantelGSS
  
  phenotypeDistances =stats::dist(phenotypeVector, method="manhattan")
  nSamples=dim(exonReads)[2] 
  sampleNames=colnames(exonReads)
  mantelList=foreach (gene_name = geneNames, .inorder=T) %dopar% {
    
    totalNumberReads=rep(0, length(geneNames))
    names(totalNumberReads)=geneNames
    for (gene_name in geneNames[1:100]){
      currGeneExonCounts= t(exonReads[which(exonGeneNames==gene_name),]) # distance is computed between rows of the matrix
      
    print(geneName)
    print(paste("Numer exons is :", as.character(dim(currGeneExonCounts)[2]), sep=""))
    totalNumberReads[gene_name]=sum(currGeneExonCounts)
    print(sum(currGeneExonCounts)) 
    }
        
    currGeneExonCounts= t(exonReads[which(exonGeneNames==gene_name),]) # distance is computed between rows of the matrix
    if (sum(currGeneExonCounts) > 0){
      exon_numbers =dim(currGeneExonCounts)[2]
      #print(exon_numbers)
      
      # compute significance of each exon individually
      exonSignif=0*(1:exon_numbers)
      for (exonIdx in 1:exon_numbers){
        currExonCounts= currGeneExonCounts[,exonIdx]
        if (sum(currExonCounts)>0){
          currExonDist=dist(currExonCounts, method="canberra")
          exonSignif[exonIdx]=vegan::mantel(phenotypeDistances, currExonDist, permutations=0, na.rm=T)$statistic
        }else{
          
          exonSignif[exonIdx]=0
        }
      }
      exonSignif[is.na(exonSignif)]=0
      
      # rank the exons in decreased order of significance
      rankedExons=length(exonSignif)-rank(exonSignif,ties.method="first")+1
      
      # record the cumulative significance of exon groups
      groupExonSignif=rep(0, exon_numbers)
      
      sequenceDistances=array(data=0,dim=c(exon_numbers,nSamples,nSamples),dimnames=list(1:exon_numbers, sampleNames, sampleNames))
      
      for (idx in 1:exon_numbers){
        exonIndexes=which(rankedExons <=idx)
        cumulativeDist= (1/idx)*as.matrix(dist(currGeneExonCounts[,exonIndexes], method="canberra"))
        cumulativeDist[is.na(cumulativeDist)]=0		
        groupExonSignif[idx]=vegan::mantel(phenotypeDistances, as.dist(cumulativeDist), permutations=0, na.rm = T)$statistic	
        sequenceDistances[idx,,]= cumulativeDist
      }
      
      mantelGSS[gene_name]=max(groupExonSignif)	
      max_idx=which(groupExonSignif==max(groupExonSignif, na.rm=T))
      max_idx= max_idx[1]
      
      mantelResults=vegan::mantel(phenotypeDistances, as.dist(sequenceDistances[max_idx,,]), permutations=nPerm, na.rm=TRUE)
      mantelResults
    } 
    # finished selecting the group of exons that maximizes splicing significance           
  }
  
  for (i in 1:length(geneNames)){
    if (!is.null(mantelList[[i]])){
    mantelGSS[i]=mantelList[[i]]$statistic
    mantelGSS_P[i]=mantelList[[i]]$signif	
    }
  }
  
  mantelGSS_P[is.na(mantelGSS_P)]=1
  mantelGSS[is.na(mantelGSS)]=0
  mantelGSS_Q=qvalue(mantelGSS_P)
  
  
  return(rbind(mantelGSS, mantelGSS_P, mantelGSS_Q$qvalues))
}

##########################################################################################


moduleEnrichment = function (colorsCoexpr, pValues)
  
{ 
  modulesCoexprNames=names(table(colorsCoexpr))
  modulesEnrichPvalues=rep(1, length(modulesCoexprNames))
  names(modulesEnrichPvalues)=modulesCoexprNames   
  modulesEnrichZscores=rep(0, length(modulesCoexprNames))
  names(modulesEnrichZscores)=modulesCoexprNames   
  
  netGenes=names(colorsCoexpr)
  pValueGenes=names(pValues)
    
  presentGenes=intersect(netGenes, pValueGenes)
  presentGeneColors=colorsCoexpr[presentGenes]
  presentPValues=pValues[presentGenes]
  
  for (module in modulesCoexprNames){
    print(module)
    moduleP=presentPValues[presentGeneColors==module]
    moduleStatistic=wilcox.test(moduleP, presentPValues, alternative="l")$statistic 
    
    randomStatistic=rep(0,1000)
    for (perm in 1:1000){
      randomGenes=sample(presentGenes, length(moduleP), replace=F)
      randomPValues=presentPValues[randomGenes]
      randomStatistic[perm]=wilcox.test(randomPValues, presentPValues, alternative="l")$statistic    
    }
     modulesEnrichZscores[module]=(moduleStatistic-mean(randomStatistic))/sd(randomStatistic)  
    modulesEnrichPvalues[module]=sum(moduleStatistic > randomStatistic)/length(randomStatistic)  
    
  }  
  
  return(rbind(length(modulesCoexprNames)*modulesEnrichPvalues, modulesEnrichZscores))
}
                    
############################################################

##########################################################################################

# this function computes z scores corresponding to differences in connectivity
# both modular and whole network

#colorsCoexpr, HSNPT_F, HDID_F, softPower=6, nPerm=10
# exprData1=HSNPT_F
# exprData2=HDID_F
# softPower=6
# nPerm=2
# nCores=14
evaluateDiffConn = function (colorsCoexpr, exprData1, exprData2, softPower, nPerm, nCores)
{
  #t1 <- Sys.time()
  
  options(cores=nCores)

  geneNames=rownames(exprData1)
    
  coexprConn1=intramodularConnectivity(adjacency(t(exprData1), power=softPower), colorsCoexpr, scaleByMax=F)
  coexprConn2=intramodularConnectivity(adjacency(t(exprData2), power=softPower), colorsCoexpr, scaleByMax=F)
  diffConnTotal=coexprConn2[, "kTotal"]-coexprConn1[, "kTotal"]
  diffConnWithin=coexprConn2[, "kWithin"]-coexprConn1[, "kWithin"]

  
  bootMatrixTotal=matrix(data=0, nrow=dim(exprData1)[1], ncol=nPerm)
  rownames(bootMatrixTotal)=geneNames
  bootMatrixWithin=bootMatrixTotal

  combinedData=cbind(exprData1, exprData2)
  
    
#   enableWGCNAThreads(nThreads=32)
#   options(cores=32)
#   
  #rm(resultListTotal)
  resultListTotal <- foreach (i=1:nPerm, .inorder=F, .errorhandling = "remove", .verbose = T) %dopar% {
 # for (perm in 1:nPerm){
   # t1 <- Sys.time()
    
    #print(perm)
    random_samples1 <- sample(1:dim(combinedData)[2],size= dim(exprData1)[2],replace=T)
    random_samples2 <- sample(1:dim(combinedData)[2],size= dim(exprData1)[2],replace=T)
    
   coexprConn1Perm=intramodularConnectivity(adjacency(t(combinedData[,random_samples1]), power=softPower), colorsCoexpr, scaleByMax=F)
    coexprConn2Perm=intramodularConnectivity(adjacency(t(combinedData[,random_samples2]), power=softPower), colorsCoexpr, scaleByMax=F)
    diffConnTotalPerm=coexprConn2Perm[, "kTotal"]-coexprConn1Perm[, "kTotal"]
    diffConnWithinPerm=coexprConn2Perm[, "kWithin"]-coexprConn1Perm[, "kWithin"]
    
#     bootMatrixTotal[,perm]=diffConnTotalPerm
#     bootMatrixWithin[,perm]=diffConnWithinPerm
    
#     t2 <- Sys.time()
#     print(difftime(t2,t1))
    
    returnVal=cbind(diffConnTotalPerm, diffConnWithinPerm)
   returnVal
    #adjacency(t(combinedData[,random_samples1]), power=softPower)
    #intramodularConnectivity(adjacency(t(combinedData[,random_samples1]), power=softPower), colorsCoexpr, scaleByMax=F)
    }
 
 
 
 for (perm in 1:nPerm){
   bootMatrixTotal[,perm]<-as.vector(resultListTotal[[perm]][,1])
   bootMatrixWithin[,perm]<-as.vector(resultListTotal[[perm]][,2])  
 }
 
 geneZTotal=rep(0, length(geneNames))
 names(geneZTotal)=geneNames
 for (i in 1:length(geneNames)) {
   geneMean=mean(bootMatrixTotal[i,])
   geneSd=sd(bootMatrixTotal[i,])
   geneZTotal[i]=(diffConnTotal[i]-geneMean)/geneSd
 }
 
 
 geneZWithin=rep(0, length(geneNames))
 names(geneZWithin)=geneNames
 for (i in 1:length(geneNames)){
   geneMean=mean(bootMatrixWithin[i,])
   geneSd=sd(bootMatrixWithin[i,])
   geneZWithin[i]=(diffConnWithin[i]-geneMean)/geneSd
 }
#  t2 <- Sys.time()
#  print(difftime(t2,t1))
 
 return(cbind(geneZTotal, geneZWithin))
 
}


#######################################################################

##########################################################################################


moduleEnrichmentConn = function (colorsCoexpr, zScores)
  
{ 
  zScores=abs(zScores)
  modulesCoexprNames=names(table(colorsCoexpr))
  modulesEnrichPvalues=rep(1, length(modulesCoexprNames))
  names(modulesEnrichPvalues)=modulesCoexprNames   
  
  
  for (module in modulesCoexprNames){
    print(module)
    moduleZ=zScores[colorsCoexpr==module]
    moduleStatistic=wilcox.test(moduleZ, zScores, alternative="g")$statistic 
    
    randomStatistic=rep(0,1000)
    for (perm in 1:1000){
      randomZ=sample(zScores, length(moduleZ), replace=F)
      randomStatistic[perm]=wilcox.test(randomZ, zScores, alternative="g")$statistic    
    }
    modulesEnrichZscores[module]=(moduleStatistic-mean(randomStatistic))/sd(randomStatistic)  
    modulesEnrichPvalues[module]=sum(moduleStatistic < randomStatistic)/length(randomStatistic)  
    
  }  
  
  return(rbind(length(modulesCoexprNames)*modulesEnrichPvalues, modulesEnrichZscores))
}


#############################################################################################

computeROC_zscores = function (goldStandardZ, comparisonZ){
  
  commonIds=intersect(names(goldStandardZ), names(goldStandardZ))  
  goldStandartZ=abs(goldStandardZ[commonIds])
  comparisonZ=abs(comparisonZ[commonIds])
  
  minZ=min(min(goldStandartZ), min(comparisonZ))
  maxZ=max(max(goldStandartZ), max(comparisonZ))
  
  zVals=seq(minZ,maxZ, (maxZ-minZ)/100)
  TP=c()
  FP=c()
  
  for (zvalue in zVals){
    
    goldStandardCurr=goldStandardZ > zvalue
    comparisonPCurr=comparisonP > zvalue
    
    FP=c(FP,sum(comparisonPCurr & (!goldStandardCurr))/sum(!goldStandardCurr))
    TP=c(TPemma,sum(comparisonPCurr & goldStandardCurr)/sum(goldStandardCurr))
  }
  return(cbind(FP, TP))
}




computeROC_pvals = function (goldStandardP, comparisonP){

  commonIds=intersect(names(goldStandardP), names(comparisonP))  
  goldStandartP=goldStandardP[commonIds]
  comparisonP=comparisonP[commonIds]
  
pVals=(1:100)/100

TP=c()
FP=c()


for (pvalue in pVals){
  
  goldStandardCurr=goldStandardP<pvalue
  comparisonPCurr=comparisonP<pvalue
  
  FP=c(FP,sum(comparisonPCurr & (!goldStandardCurr))/sum(!goldStandardCurr))
  TP=c(TPemma,sum(comparisonPCurr & goldStandardCurr)/sum(goldStandardCurr))
  
  
}

return(cbind(FP, TP))

}


plotROC = function (ROCresultList, colors, legendText){

  plot(ROCresultList[[1]][,1], ROCresultList[[1]][,2],type="l",col=colors[1], xlim=c(0,1),ylim=c(0,1), lwd=4,xlab="False Positives",ylab="True Positives",cex.lab=1.7, axes=FALSE)
  text(x=0.7,y=0.1,legendText[1],col=colors[1],cex=1.7)
  
  if (length(ROCresultList) > 1){
    for (i in 2:length(ROCresultList)){
      points(ROCresultList[[i]][,1], ROCresultList[[i]][,2],type="l",col=colors[i], xlim=c(0,1),ylim=c(0,1), lwd=4,xlab="False Positives",ylab="True Positives",cex.lab=1.7, axes=FALSE)
      text(x=0.7,y=0.1+i*0.1,legendText[1],col=colors[i],cex=1.7)
    }
  }


points((1:100)/100,(1:100)/100,type="l", lwd=4)

axis(1,at=c(0,0.25,0.5,0.75,1), c(0,0.25,0.5,0.75,1),cex.axis=1.5, lwd=1.5)
axis(2,at=c(0,0.25,0.5,0.75,1), c(0,0.25,0.5,0.75,1),cex.axis=1.5, lwd=1.5)
box() #- to make it look "as usual"

}
