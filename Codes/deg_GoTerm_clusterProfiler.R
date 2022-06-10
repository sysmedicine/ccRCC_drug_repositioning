deg_GoTerm_clusterProfiler<-function(clusterProfile_result){
  #we set the fdr_GoTerm for controling p.adjust
  
  path_result<-clusterProfile_result@result
  GO_ID<-as.matrix(path_result$ID)
  term_name<-as.matrix(path_result$Description)
  pvalue<-as.matrix(path_result$pvalue)
  p.adjust<-as.matrix(path_result$p.adjust)
  qvalue<-as.matrix(path_result$qvalue)
  #Count<-as.matrix(path_result$Count)
  geneID<-as.matrix(path_result$geneID)
  path_out_table<-cbind(GO_ID,term_name,pvalue,p.adjust,qvalue,geneID)
  colnames(path_out_table)<-c("GO_ID","term_name","pvalue","p.adjust","qvalue","geneID")
  
  path_out_table<-as.data.frame(path_out_table)
  mode(path_out_table$pvalue)="numeric"
  mode(path_out_table$p.adjust)="numeric"
  mode(path_out_table$qvalue)="numeric"
  return(path_out_table)
}
