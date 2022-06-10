Cheng_initialRenvironment <- function(Update){
  if (!hasArg("Update")) {
    Update = F
  }
  if (Update == T) {
    install.packages("xlsx")
    install.packages("XLConnect")
    install.packages("ggplot2")
    install.packages("ggvis")
    install.packages("rgl")
    install.packages("dplyr")
    install.packages("tidyr")
    install.packages("stringr")
    install.packages("lubridate")
    install.packages("gplots")
    install.packages("httr")
    install.packages("XML")
    install.packages("survival")
    install.packages("foreign")
    install.packages("checkmate")
    install.packages("RSQLite")
    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    biocLite("limma")
    biocLite("DESeq2")
    #biocLite("BiocUpgrade")
    biocLite("piano", dependencies=TRUE)
    biocLite("GO.db")
    biocLite("Biobase")
    biocLite("biomaRt")
  }
  
  library("grid")
  library("xlsx")
  library("XLConnect")
  library("ggplot2")
  library("ggvis")
  library("rgl")
  library("dplyr")
  library("tidyr")
  library("stringr")
  library("lubridate")
  require(gplots)
  library("survival")
  require("xlsx")
  library("DESeq2")
  library("biomaRt")
  library("DESeq2")
  library("piano")
  library("Biobase")
}

Cheng_deseq2 <- function(Y, group,filter=FALSE,returnEffect=FALSE){
  
  
  
  dds2 <-DESeqDataSetFromMatrix(countData=Y,colData=as.data.frame(group),design=~group)
  
  dds2 <- estimateSizeFactors(dds2)
  dds2 <- estimateDispersions(dds2)
  dds2 <- nbinomWaldTest(dds2)
  # if (filter) {
  #   res <- results(dds2)
  # } else {
  #   res <- results(dds2,cooksCutoff=FALSE, independentFiltering=FALSE)
  # }
  # 
  # if ( returnEffect) {
  #   return(list(res[,5],2^res[,2]))
  # } else {
  #   return(res[,5])
  # }
}




Cheng_readFile <- function(fileName,sepMark) { 
  if (!hasArg("sepMark")) { 
    sepMark = '\t'
  }
  data = read.table(fileName,sep = sepMark,stringsAsFactors = F)
  rowNames = data[2:length(data[,1]),1]
  colNames = data[1,2:length(data[1,])]
  data = data[,-1]
  data = data[-1,]
  rownames(data) = rowNames
  colnames(data) = colNames
  return(data)
  
  
  
}

Cheng_readMatrixFile <- function(fileName,sepMark) { 
  if (!hasArg("sepMark")) { 
    sepMark = '\t'
  }
  data = read.table(fileName,sep = sepMark,stringsAsFactors = F)
  rowNames = data[2:length(data[,1]),1]
  colNames = data[1,2:length(data[1,])]
  data = data[,-1]
  data = data[-1,]
  data = as.matrix(sapply(data, as.numeric))  
  rownames(data) = rowNames
  colnames(data) = colNames
  return(data)
  
  
  
}

Cheng_generateKMplot <- function(SurvInput,quantCut,outFile,figureTitle) { 
  
  
#  content = read.table("ASS1.txt",sep = '\t',header = TRUE)
  
#  DeadInd = content$Status %in% c('dead')
#  LivingDays = content$LivingDays
#  SurvInput= as.data.frame(LivingDays)
#  SurvInput$DeadInd = DeadInd
#  EXP = content$ENSG00000130707
#  SurvInput$EXP = EXP
  
  
  if (!hasArg("quantCut")) {
    cutoffs = sort(unique(SurvInput$EXP))
    Percentile20 = quantile(SurvInput$EXP,0.2)
    Percentile80 = quantile(SurvInput$EXP,0.8)
    cutoffs = cutoffs[cutoffs>Percentile20&cutoffs<=Percentile80]
    Pvalue = 1
    cutoff = 0
    mini = 0
    Pvalues = cutoffs
    for (i in 1:length(cutoffs)){
      tempdata = SurvInput
      tempdata$EXP[tempdata$EXP<cutoffs[i]]=0
      tempdata$EXP[tempdata$EXP>=cutoffs[i]]=1
      survOut = survfit(SurvObj ~ EXP,tempdata)
      res.km2 <- survdiff(SurvObj~ EXP, tempdata, rho=0)
      icutoff = cutoffs[i]
      iPvalue = pchisq(res.km2$chisq,length(res.km2$n)-1,lower.tail = FALSE)
      Pvalues[i] = iPvalue
      if (iPvalue < Pvalue) {
        cutoff = icutoff
        Pvalue = iPvalue
        mini = i
      }
      
      
      
      
    }
  } else {
    cutoff = quantile(SurvInput$EXP,quantCut)
  }
  SurvInput$EXP[SurvInput$EXP<cutoff]=0
  SurvInput$EXP[SurvInput$EXP>=cutoff]=1
  survOut = survfit(SurvObj ~ EXP,SurvInput)
  
  
  pdf(file = paste(outFile,"pdf",sep = "."),width=5.8,height=6)
  plot(survOut, col=c("red","black"), mark.time=T, cex=1.4,xlab="Time (year)",xscale = 365,lty =1, ylab = "Survival Probability",las=1, cex.lab=1.4)
  group1legend = paste("Low expression",paste("(n =",as.character(sum(SurvInput$EXP==0)),")",sep = ""), sep = " ")
  group2legend = paste("High expression",paste("(n =",as.character(sum(SurvInput$EXP==1)),")",sep = ""), sep = " ")
  if (!hasArg("figureTitle")) {
    #title("ASS1")
  } else {
    title(figureTitle)
  }
  legend(
    "bottomleft",
    legend=c(group1legend,group2legend),
    col=c("red","black"),
    horiz=FALSE,
    lty=1,
    bty='n',
    cex=1.4)
#  title("ASS1") 
  res = survdiff(SurvObj ~ EXP, SurvInput)
  logRankP = 1 - pchisq(res$chisq, length(res$n)-1)
  legend("topright",legend =c(paste0("P=",as.character(format(logRankP,scientific = TRUE,digits = 3)))),text.font=2,bty="n",cex=1.4)
  dev.off()
  sum_result<-summary(coxph(SurvObj ~ EXP, SurvInput))
  coef<-sum_result$coefficients[1]
  result = as.data.frame(logRankP)
  result$EXPcut = cutoff
  result$coef=coef
  return(result)
}


Cheng_generateSurvInputfromTCGA <- function(gene,cancerType,dataDir) { 
  #if (!hasArg("dataDir")) {
  #  dataDir = paste0("C:/Pan_cancer/TCGA_transcript_result_store/",cancerType,"/PKLR/")
  #}
  data = Cheng_readFile(paste0(dataDir,paste0(cancerType,"_trans_exp_TPM_1.txt")))
  LivingDays = as.numeric(data$LivingDays)
  SurvInput= as.data.frame(LivingDays)
  SurvInput$EXP = as.numeric(data[,gene])
  SurvInput$DeadInd = data$Status %in% c('dead')
  SurvInput$SurvObj <- with(SurvInput, Surv(LivingDays, DeadInd))
  return(SurvInput)
}

Cheng_DEanalysis <- function(inputFile,groupInd,outDir,NameCondA,NameCondB,gmtFile) {
  # groupInd example c(rep(0,2),rep(1,5),rep(2,3)) in which 0 means not used, 1 means condition A and 2 means condition B
  # Always Condition B compared to Condition A, that is to say, up-/down-regulation means up-/down-regulated in condition B
  if (!hasArg("outDir")) { 
    outDir = 'C:/work/R scripts/R dump/'
  }
  if (!hasArg("NameCondA")) { 
    NameCondA = 'CondA'
  }
  if (!hasArg("NameCondB")) { 
    NameCondB = 'CondB'
  }
  if (!hasArg("gmtFile")) {
    gmtFile = "C:/work/Database/ENTREZ/c5.bp.v5.0.entrez.gmt"
  }
  data = Cheng_readMatrixFile(inputFile)
  ind = which(groupInd==0)
  data = data[,-ind]
  data = round(data)
  groupInd = groupInd[-ind]
  groupInd = groupInd - 1 
  result = Cheng_deseq2(data,groupInd)
  outFileName = paste(paste(outDir,NameCondB,"vs",NameCondA,sep = ""),"csv",sep = ".")
  write.csv(results(result),file = outFileName)
  
  DEdata = Cheng_readMatrixFile(outFileName,sep = ",")
  
  ## Convert Ensemble to entrez
  mart = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
  #listDatasets(mart)
  ENSB2ETREZ = getBM(attributes = c("ensembl_gene_id","entrezgene"),values = rownames(DEdata),filters = "ensembl_gene_id",mart = mart)
  ETREZ2HGNC = getBM(attributes = c("entrezgene","hgnc_symbol"),values = rownames(DEdata),filters = "ensembl_gene_id",mart = mart)
  res = DEdata[match(  ENSB2ETREZ[,1], rownames(DEdata)),]
  res2 = DEdata[match(  ENSB2HGNC[,1], rownames(DEdata)),]
  
  #GO ANALYSIS
  ## "Download GMT Files": entrez genes ids
  myGsc <- loadGSC(gmtFile) #Your own path
  #head(myGsc)
  
  pval= res[ ,6]
  pval = as.matrix(pval)
  rownames(pval) = paste(as.matrix(ENSB2ETREZ[ ,2] ) )
  pval[is.na(pval)] <- 1
  fc= res[,2]  #extract fold changes
  fc = as.matrix(fc)
  fc[fc>0] <- 1
  fc[fc<0] <- -1
  fc[is.na(fc)] <- 0
  fc = as.numeric(fc)
  fc = as.matrix(fc)
  rownames(fc) = paste(as.matrix(ENSB2ETREZ[ ,2] ) )
  
  #check everything OK:
  head(pval)
  head(fc)
  
  ## Run GSEA
  gsaRes <- runGSA(pval, fc, gsc=myGsc, geneSetStat="reporter", signifMethod="nullDist", gsSizeLim=c(10,10000))
  
  ## Write output to excel file
  outFileNamePIANO = paste(paste(outDir,"PIANO_output_",NameCondB,"vs",NameCondA,sep = ""),"xls",sep = ".")
  GSAsummaryTable(gsaRes, save=TRUE, file=outFileNamePIANO)
  
  ## network plot
  #par(mar = rep(2, 4))
  #networkPlot(gsaRes, class='non', significance=0.0000005,lay=5, ncharLabel=75)
  
  ## heatmap
  outFileNamePIANOpdf = paste(paste(outDir,"PIANO_output_",NameCondB,"vs",NameCondA,"_heatmap",sep = ""),"pdf",sep = ".")
  pdf(outFileNamePIANOpdf ,width=20, height=20)
  GSAheatmap(gsaRes, cutoff=25, adjusted=FALSE, ncharLabel=75, cellnote="none", columnnames="full", colorkey=TRUE, colorgrad=NULL, cex=0.5)
  dev.off()
  
  
}

Cheng_coxForEXP<- function(SurvInput) {
  cox = coxph(SurvObj ~ EXP, SurvInput)
  
}
