require(matrixStats)

require(data.table)

#install.packages('igraph')

require(igraph)

require(Hmisc)

require(compiler)

require(reshape2)



#filterExprMat <- function(exprMat, meanCut, mode="mean") {

#if (mode == "mean") rMeans = rowMeans(exprMat)

#if (mode =="median") rMeans = rowMedians(exprMat)

#return(exprMat[rMeans>meanCut,]) }



makeCorTable <- function(exprMat, cutoff= 0.99, mode="pearson", self=F, debug= F){
  #cor function alwary calculate the the co-expression between column elements 
  minValue = -2
  
  if(mode == "pearson"){
    
    corMat = cor(t(exprMat))
    
  }
  
  if(mode == "incomplete"){
    
    corMat = cor(t(exprMat), use="pairwise.complete.obs")
    
  }
  
  if(mode == "spearman"){
    
    corMat =cor(t(exprMat), method = "spearman")
    
  }
  
  if (mode == "r2") {
    
    corMat = cor(t(exprMat))
    
    corMat = corMat^2
    
  }
  
  if (debug) print("na amount:")
  
  if (debug) print(table(is.na(corMat)))
  
  #corMat[is.na(corMat)] = minValue
  corMat[is.na(corMat)] = 0
  cutoff = as.numeric(quantile(corMat[lower.tri(corMat)], cutoff, na.rm = T))
  #order from the smallest to largest, so get the high 0.99 quantile if quantile=0.99
  corMat[which(corMat<cutoff)]=0
  return(corMat)
  
  #if(debug) print("quantile: ")
  
  #if(debug) print(cutoff)
  
  #resTable = melt(corMat)[,1:3]
  
  #resTableNew = as.data.table(resTable)
  
  #setkey(resTableNew, value)
  
  #resTableNew = resTableNew[order(value, decreasing = T)]
  
  #resTable = as.data.frame(resTableNew)
  
  #if (!self) {
    
   # resTable = resTable[resTable[,1] != resTable[,2],]
    
 # }
  
  #if(debug) print("cut-off before")
  
  #if(debug) print(dim(resTable))
  
  #resTable =resTable[resTable[,3]>cutoff,]
  
  #if(debug) print("cut-off on the way")
  
  #if(debug) print(dim(resTable))
  
  #resTable = resTable[!is.na(resTable[,3]),]
  
  #resTable = resTable[resTable[,3] != minValue,]
  
  #if(debug) print("cutoff after")
  
  #if(debug) print(dim(resTable))
  
  #return(resTable)
  
}



makeCorNet <- function(corMat) {
  
  #  corNet = igraph::graph.data.frame(corTable[,1:2], directed = F)
  corNet = igraph::graph_from_adjacency_matrix(corMat, mode = "undirected", weighted = TRUE,diag = F, add.colnames=T, add.rownames=T)
  corNet = set_vertex_attr(corNet,'name',value = rownames(corMat))
  corNet = delete_vertex_attr(corNet,'TRUE')
  corNet = igraph::simplify(corNet, remove.multiple = T, remove.loops = T)
  
  return(corNet)
  
}



getModulePruned <- function(moduleList, cutCluster, debug=T) {
  
  currGeneClusterLen = unlist(lapply(moduleList, length))
  
  names(currGeneClusterLen) = names(moduleList)
  
  if (debug) print(table(currGeneClusterLen>=cutCluster))
  
  currGeneClusterCut = names(currGeneClusterLen[currGeneClusterLen>=cutCluster])
  
  moduleList = moduleList[currGeneClusterCut]
  
  return(moduleList)  
  
}



getModuleMatrix <- function(corNet, modulePruned, debug=F) {
  
  numEdges = ecount(corNet)
  
  adjMat = as_adj(corNet)
  
  degMat = getDegMat(corNet)
  
  
  
  numModules = length(modulePruned)
  
  if (debug) print(numModules)
  
  moduleMat = matrix(0, numModules, numModules)
  
  for (i in 1:(numModules-1)) {
    
    if (debug) print(i)
    
    for (j in (i+1):numModules) {
      
      moduleMat[i,j] = areModuleInteract(adjMat, degMat, numEdges, modulePruned[[i]], modulePruned[[j]])
      
    }
    
  }  
  
  moduleMat[lower.tri(moduleMat)] = t(moduleMat)[lower.tri(moduleMat)]
  
  rownames(moduleMat) = names(modulePruned)
  
  colnames(moduleMat) = names(modulePruned)
  
  return(moduleMat)
  
}



getDegMat <- function(igraphNet) {
  
  #we regard network as undirected
  
  degs = igraph::degree(igraphNet, mode="in", loops=F)
  
  #if using undirected graph mode argument is ignored and loops removed
  
  #degs = igraph::degree(igraphNet, loops= T)
  
  names(degs) = V(igraphNet)$name
  
  degMat = as.matrix(degs) %*% t(as.matrix(degs))
  
  return(degMat)
  
}



areModuleInteract <- function(adjMat, degMat, numEdges, genesA, genesB) {
  
  observedEdges = observedEdgesBtw(adjMat, genesA, genesB)
  
  expectedEdges = expectedEdgesBtw(degMat, numEdges, genesA, genesB)
  
  #diff=(observedEdges-expectedEdges)/expectedEdges
  
  diff=observedEdges/expectedEdges
  
  return(diff)
  
}







annotateModulesByCC <- function(corNet, moduleList, cutCluster=5, cutCC = 0.5, debug=F){
  
  modulePrunedList = getModulePruned(moduleList, cutCluster)
  
  moduleNames = names(modulePrunedList)
  
  moduleMat = getModuleMatrix(corNet, modulePrunedList)
  
  print(moduleMat)
  
  moduleCC= getClusterScores(corNet, modulePrunedList)
  
  highCCModuleNames = names(moduleCC[moduleCC>=cutCC])
  
  summaryTypes = rep("None", length(modulePrunedList))
  
  names(summaryTypes) = moduleNames
  
  summaryTypes[moduleNames %in% highCCModuleNames] = "HighCC"
  
  summaryTypes[!moduleNames %in% highCCModuleNames] = "LowCC"
  
  
  
  moduleCut=1
  
  moduleMelt = melt(moduleMat)
  
  moduleMelt = moduleMelt[moduleMelt[,3]>moduleCut,]
  
  moduleMelt = lapply(moduleMelt, as.character)
  
  moduleSizes = unlist(lapply(modulePrunedList, length))
  
  moduleAttr = data.frame(node=moduleNames, size=moduleSizes, summary = summaryTypes, cc = moduleCC, stringsAsFactors = F)
  
  modulesNotShown = moduleNames[!moduleNames %in% moduleAttr]
  
  return(list(nodeTable = moduleAttr, edgeTable=moduleMelt, nodesNotShown= modulesNotShown))
  
}



observedEdgesBtw <- function(adjMat, genesA, genesB) {
  
  adjMatSel = adjMat[genesA, genesB]
  
  count = sum(adjMatSel)
  
  return(count)
  
}



expectedEdgesBtw <- function(degMat, numEdges, genesA, genesB) {
  
  degMatSel = degMat[genesA, genesB]/(2*numEdges)
  
  count = sum(degMatSel)
  
  return(count)
  
}



getClusterScores <- function(corNet, geneClusterList) {
  
  ccScores = lapply(geneClusterList, function(currClusterGenes){
    
    #currGraph = igraph::subgraph(corNet, currClusterGenes)
    
    #currGraph = subgraph(corNet, currClusterGenes)
    
    currGraph = induced.subgraph(corNet, currClusterGenes) # use this one if subgraph is giving you an error
    
    cc = transitivity(currGraph, "globalundirected")
    
    return(cc)
    
  })
  
  ccScores = unlist(ccScores)
  
  names(ccScores) = names(geneClusterList)
  
  return(ccScores)
  
}



# write THE FILES FOR CYTOSCPAE

write.edge.cytoscape <- function(edgeTable, nodesNotShown, outFile) {
  
  #lenTable=dim(edgeTable)[1]
  
  lenTable=3
  
  #lenNodesNotShown = length
  
  print(head(edgeTable))
  
  sink(outFile, append=F)
  
  cat("source")
  
  cat("\t")
  
  cat("target")
  
  cat("\n")
  
  
  
  for (i in 1:lenTable) {
    
    # ints = unlist(edgeTable[i,])
    
    # cat(edgeTable[i,1])
    
    cat(edgeTable$Var1)
    
    cat("\t")
    
    # cat(edgeTable[i,2])
    
    cat(edgeTable$Var2)
    
    cat("\n")
    
  }
  
  for (node in nodesNotShown) {
    
    cat(node)
    
    cat("\n")
    
  }
  
  sink() 
  
}



write.node.cytoscape <- function(nodeTable, outFile) {
  
  write.table(nodeTable, outFile, row.names = F, quote = F, sep="\t")
  
}



#makeModuleList <- function(corNet,debug=F) {
  
  #fc = cluster_walktrap(corNet)
  
  #cluster = fc$membership
  
  #geneCluster = data.frame(gene=get.vertex.attribute(corNet)$'TRUE'[as.vector(V(corNet))], cluster=cluster, stringsAsFactors = F)
  
  #geneClusterList = split.data.frame(geneCluster, geneCluster$cluster)
  
  #geneClusterList = lapply(geneClusterList, "[[", "gene")
  
  #geneClusterSizes = do.call(c, lapply(geneClusterList, length))
  
  #if(debug) {
    
  # print("... done")
    
  #  print(paste("... modularity:", as.character(modularity(fc))))
    
  #  print(paste("... no clusters:", as.character(length(geneClusterList))))
    
  #  print(paste("... no of genes max cluster: ", as.character(sort(geneClusterSizes,T)[1])))
    
  #}
  
  #return(geneClusterList)
  
#}

makeModuleList <- function(corNet, debug=F) {
  
  fc = cluster_walktrap(corNet)
  
  cluster = fc$membership
  
  geneCluster = data.frame(gene=V(corNet)$name, cluster=cluster, stringsAsFactors = F)
  
  geneClusterList = split.data.frame(geneCluster, geneCluster$cluster)
  
  geneClusterList = lapply(geneClusterList, "[[", "gene")
  
  geneClusterSizes = do.call(c, lapply(geneClusterList, length))
  
  if(debug) {
    
    print("... done")
    
    print(paste("... modularity:", as.character(modularity(fc))))
    
    print(paste("... no clusters:", as.character(length(geneClusterList))))
    
    print(paste("... no of genes max cluster: ", as.character(sort(geneClusterSizes,T)[1])))
    
  }
  
  return(geneClusterList)
  
}



#exprMat = filterExprMat(tissueTpmMat, meanCut = 1, mode = "mean")


#corTable = makeCorTable(exprMat, cutoff = 0.99, mode = "pearson", self=F, debug =F)

#corNet = makeCorNet(corTable)

#moduleList= makeModuleList(corNet, debug = F)

#cytoscapematerial = annotateModulesByCC(corNet, moduleList, cutCluster = 20, cutCC = 0.5, debug = F)


#aaa=as.data.frame(cytoscapematerial$edgeTable)

#write.table(aaa, file = "cytoEdgefile.txt", sep = "\t",row.names = F, col.names = T)

#write.node.cytoscape(cytoscapematerial$nodeTable, 'cytoNodefile.txt')

