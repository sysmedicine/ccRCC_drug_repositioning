#########################################################################################
#This codes include four main parts as follows:
#Kaplan-Meier analysis
#Go Term enrichment 
#co-expression network
#Drug repositioning
#########################################################################################
#Kaplan-Meier analysis
#TPM expression profiles were used as input.
#Find the example data in the folder "Example_data": TCGA_KIRC_trans_exp_TPM_1.txt
library(survival)
setwd("C:\\Codes\\")
source('cheng_toolbox_beta_1.R')

cancerType<-"TCGA_KIRC"

path_input<-"C:\\Example_data\\"
setwd(path_input)
exp<-as.matrix(read.table(paste0(cancerType,"_trans_exp_TPM_1.txt"),header=T,sep="\t"))
exp<-exp[,-c(1:7)]
gene_list<-as.matrix(colnames(exp))

dataDir=path_input
path_out_pdf<-"C:\\Example_data\\output_pdf\\"
path_out_stat<-"C:\\Example_data\\output_stat\\"

p_list<-NULL
cutoff_list<-NULL
coef_list<-NULL
for (j in 1:length(gene_list)){
  print(j)
  output<-Cheng_generateSurvInputfromTCGA(gene_list[j,1],cancerType,dataDir)
  setwd(path_out_pdf)
  result<-Cheng_generateKMplot(output,outFile=paste0(cancerType,"_",gene_list[j,1]))
  log_rank_p<-as.matrix(result$logRankP)
  cut_off<-as.matrix(result$EXPcut)
  coef<-as.matrix(result$coef)
  
  p_list<-rbind(p_list,log_rank_p)
  cutoff_list<-rbind(cutoff_list,cut_off)
  coef_list<-rbind(coef_list,coef)
  rm(log_rank_p,cut_off,coef)
}
final_result<-cbind(as.matrix(gene_list),coef_list,cutoff_list,p_list)
colnames(final_result)<-c("gene","coef","cutoff","p")

setwd(path_out_stat)
write.table(final_result,file=paste0(cancerType,"_KM_stat_result_1.txt"),sep="\t",row.names=F,col.names=T,quote=F)
#########################################################################################
#Go Term enrichment for marker genes (overlapped prognostic genes)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)

setwd("C:\\Codes\\")
source('deg_GoTerm_clusterProfiler.R')

path_TCGA<-"C:\\Pan_cancer\\TCGA_transcript\\TCGA-KIRC\\tpm\\gene_level_exp\\"
setwd(path_TCGA)
load("exp_tpm_gene_level_rem1_co.Rdata")#load symbol_exp (your expression profiles, rownames is gene symbol)
background<-as.matrix(rownames(symbol_exp))
rm(ensem_of_symbol_exp,clinical_exp,symbol_exp)

path_raw<-"E:\\L1000_2021\\kidney_drug_reposition\\KIRC_marker_gene\\"
setwd(path_raw)
int_gene<-as.matrix(read.csv("KIRC_markers.txt",header=T,sep="\t"))
int_gene<-as.matrix(int_gene[,"symbol"])#get the gene symbol
pathway=enrichGO(gene=int_gene,OrgDb='org.Hs.eg.db',keyType = "SYMBOL",ont="BP",universe = background ,pAdjustMethod = "BH",qvalueCutoff=0.05)
path_table_all<-deg_GoTerm_clusterProfiler(pathway)
setwd(path_raw)
write.table(path_table_all,file="GO_terms_KIRC_markers.txt",sep="\t",row.names=F,col.names=T,quote=F)
#########################################################################################
#co-expression network----------random walk
##extract the module form top 1% network, neworkson-----sunjie code
setwd("C:\\Codes\\")
source('networkson_change.R')

top_quantile=0.99 #extract the co-expressed gene pair with r values within top 1%
num_nodes=10
transitivity=0.5
exp_th=5

setwd("C:\\Pan_cancer\\TCGA_transcript\\TCGA-KIRC\\tpm\\gene_level_exp\\")
load("exp_tpm_gene_level_rem1_co.Rdata")#load symbol_exp (your expression profiles, row names are gene symbol)
mean_value<-as.matrix(rowMeans(symbol_exp))
index<-which(mean_value>exp_th)
symbol_exp<-symbol_exp[index,]

path_out<-"E:\\L1000_2021\\kidney_drug_reposition\\coexp_network\\based_all_genes\\TPM_5\\result_TCGA_top1\\"
setwd(path_out)

corMatrix = makeCorTable(symbol_exp, cutoff = top_quantile, mode = "spearman", self=F, debug =F)##Create co-expression matrix (corMatrix)
save(file="coexp_network_top1.Rdata",corMatrix)

corNet = makeCorNet(corMatrix)#Create the network graph based on the co-expression matrix
moduleList= makeModuleList(corNet, debug = F) #Use random walk method to extract gene modules
save(file="module_list.Rdata",moduleList)

cytoscapematerial = annotateModulesByCC(corNet, moduleList, cutCluster = num_nodes, cutCC = transitivity, debug = F)#Extract the modue list with number of nodes >10; If connectivity (cc) > 0.5, label as "HighCC", otherwise "LowCC".

write.node.cytoscape(cytoscapematerial$nodeTable, 'cytoNodefile.txt')#output the module list
#########################################################################################
#Drug Repositioning
#Spearman correlation between shRNA and drug perturbed signature profiles for cell line HA1E
#We have extracted the shRNA/compound perturbed signature profiles (level 5) from CMap database (https://clue.io/data/CMap2020#LINCS2020), find the data in the folder "Example_data"(HA1E_sh.gctx and HA1E_cp.gctx).
library(cmapR)

int_cell="HA1E"
int_target=as.matrix(c("BUB1B","CCNB2","CEP55","ASF1B","RRM2"))

path_sig<-"E:\\L1000_2021\\metadata\\sig_info_process\\"
setwd(path_sig)
load("siginfo_beta_yuan.Rdata")#load sig_id, find it in the folder "Example_data" ,which same as the "siginfo_beta.txt" downloaded from CMap database 2020 (https://clue.io/data/CMap2020#LINCS2020)

setwd("E:\\L1000_2021\\metadata\\gene_info_process\\")
load("geneinfo_proteinCode_vs_KIRCexpGenes_remDup.Rdata")#load gene_info (not necessary)

path_out<-"E:\\L1000_2021\\kidney_drug_reposition\\coexp_network\\based_all_genes\\TPM_1\\Drug_reposition\\Cor_SHvsCP_HA1E\\"

path_sh_gctx<-"E:\\L1000_2021\\shRNA_processed\\"
setwd(path_sh_gctx)
data_sh<-parse_gctx(paste0(int_cell,"_sh.gctx"))
exp_sh<-data_sh@mat
exp_sh<-exp_sh[gene_info[,"gene_id"],]

path_cp_gctx<-"E:\\L1000_2021\\drug_processed\\"
setwd(path_cp_gctx)
data_cp<-parse_gctx(paste0(int_cell,"_cp.gctx"))
exp_cp<-data_cp@mat
exp_cp<-exp_cp[gene_info[,"gene_id"],]

rm(data_sh,data_cp)
gc()

for (i in 1:length(int_target)){
  print(i)
  index_sh<-which(sig_info[,"pert_type"]=="trt_sh"&sig_info[,"cell_iname"]==int_cell&sig_info[,"cmap_name"]==int_target[i])
  sig_id_sh<-as.matrix(unique(sig_info[index_sh,"sig_id"]))#extract the sig_id for shRNA perturbation for the specific target gene
  int_exp_sh<-exp_sh[,sig_id_sh]
  
  r_matrix<-NULL
  for (j in 1:dim(int_exp_sh)[2]){
    r_each<-t(cor(int_exp_sh[,j],exp_cp,method="spearman"))
    colnames(r_each)<-colnames(int_exp_sh)[j]
    r_matrix<-cbind(r_matrix,r_each)
    rm(r_each)
    gc()
  }
  
  path_out_each<-paste0(path_out,int_target[i],"\\")
  dir.create(path_out_each)
  setwd(path_out_each)
  save(file=paste0("r_cor_SHvsCP_",int_target[i],"_end.Rdata"),r_matrix)
  
  rm(index_sh,sig_id_sh,int_exp_sh,r_matrix,path_out_each)
  gc()
}

#########################################################################################
#Filter the coefficient matrix (obtained from above step) by selecting each drug with best dose and treat time. In the experimental setting in CMap, they used different doses and treated time to treat a cell line. Here we extract the best drug setting with highest correlation coefficient.

rm(list=ls())

int_target=as.matrix(c("BUB1B","CCNB2","CEP55","ASF1B","RRM2"))

path_sig<-"E:\\L1000_2021\\metadata\\sig_info_process\\"
setwd(path_sig)
load("siginfo_beta_yuan.Rdata")#load sig_id, same as the "siginfo_beta.txt" downloaded from CMap database 2020 (https://clue.io/data/CMap2020#LINCS2020)
pertID_to_cmapName<-unique(sig_info[,c("pert_id","cmap_name")])

path_raw<-"E:\\L1000_2021\\kidney_drug_reposition\\coexp_network\\based_all_genes\\TPM_1\\Drug_reposition\\Cor_SHvsCP_HA1E\\"

for (i in 1:length(int_target)){
  print(i)
  path_each<-paste0(path_raw,int_target[i],"\\")
  setwd(path_each)
  load(paste0("r_cor_SHvsCP_",int_target[i],"_end.Rdata"))#load r_matrix
  index_cp<-match(rownames(r_matrix),sig_info[,"sig_id"])
  sig_info_cp<-as.matrix(sig_info[index_cp,])# mapped to r_matrix
  uni_cp_id<-as.matrix(unique(sig_info_cp[,"pert_id"]))
  
  r_matrix_end<-matrix(NA,length(uni_cp_id),dim(r_matrix)[2])
  sig_id_matrix_end<-matrix(NA,length(uni_cp_id),dim(r_matrix)[2])
  for (n in 1:length(uni_cp_id)){
    
    ind<-which(sig_info_cp[,"pert_id"]==uni_cp_id[n])
    if (length(ind)>1){
      int_r_matrix<-r_matrix[ind,]
      int_sig_info_cp<-sig_info_cp[ind,]
      for (m in 1:dim(r_matrix)[2]){
        ind_2<-which(int_r_matrix[,m]==max(int_r_matrix[,m]),arr.ind=T)
        r_matrix_end[n,m]=int_r_matrix[ind_2,m]
        sig_id_matrix_end[n,m]=int_sig_info_cp[ind_2,"sig_id"]   
      }
    }
    
    if(length(ind)==1){
      int_r_matrix<-t(as.matrix(r_matrix[ind,]))
      int_sig_info_cp<-t(as.matrix(sig_info_cp[ind,]))
      r_matrix_end[n,]=int_r_matrix[1,]
      sig_id_matrix_end[n,]=int_sig_info_cp[1,"sig_id"]   
      
    }
    
  }
  
  loc<-match(uni_cp_id,pertID_to_cmapName[,"pert_id"])
  r_matrix_end<-cbind(pertID_to_cmapName[loc,],r_matrix_end)
  colnames(r_matrix_end)<-c(colnames(pertID_to_cmapName),colnames(r_matrix))
  
  sig_id_matrix_end<-cbind(pertID_to_cmapName[loc,],sig_id_matrix_end)
  colnames(sig_id_matrix_end)<-c(colnames(pertID_to_cmapName),colnames(r_matrix))
  
  write.table(r_matrix_end,file=paste0("r_matrix_sel_maxR_",int_target[i],".txt"),sep="\t",row.names=F,col.names=T,quote=F)
  write.table(sig_id_matrix_end,file=paste0("sig_id_cp_sel_maxR_",int_target[i],".txt"),sep="\t",row.names=F,col.names=T,quote=F)
  
  loc_2<-match(colnames(r_matrix),sig_info[,"sig_id"])
  sh_sig_info<-sig_info[loc_2,]
  write.table(sh_sig_info,file=paste0("sig_id_sh_sek_maxR_",int_target[i],".txt"),sep="\t",row.names=F,col.names=T,quote=F)
  
}

#########################################################################################
#Filter the coefficient matrix by selecting the best subset of shRNAs. In the experimental setting in CMap, they used different shRNAs for the knockdown of a gene. Here we extract the a subset of shRNAs which lead to similar effect on cells by using clustering analysis. Based on the clustering figure and median coefficient, choose the subset of shRNAS.
#Heatmap--clustering of shRNAs
rm(list=ls())

library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(ggpubr)

int_target=as.matrix(c("BUB1B","CCNB2","CEP55","ASF1B","RRM2"))

path_raw<-"E:\\L1000_2021\\kidney_drug_reposition\\coexp_network\\based_all_genes\\TPM_1\\Drug_reposition\\Cor_SHvsCP_HA1E\\"

my.breaks <- c(seq(-1, -0.0001, by=0.001), seq(0, 1, by=0.001))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2),
               colorRampPalette(colors = c("white","red"))(length(my.breaks)/2))

for (i in 1:length(int_target)){
  
  path_each<-paste0(path_raw,int_target[i],"\\")
  setwd(path_each)
  
  r_matrix<-as.matrix(read.csv(paste0("r_matrix_sel_maxR_",int_target[i],".txt"),header=F,sep="\t"))
  colnames(r_matrix)<-r_matrix[1,]
  r_matrix<-r_matrix[-1,]
  
  drug_inform<-r_matrix[,1:2]
  r_matrix<-r_matrix[,-c(1,2)]
  mode(r_matrix)="numeric"
  
  cor.matrix=rcorr(r_matrix, type="spearman")
  cor.matrix.R=cor.matrix$r
  
  tiff(file=paste0("pheatmap2_r_matrix_sel_maxR_",int_target[i],".tiff"),height =800 ,width = 800)
  
  cols.cor <- cor(cor.matrix.R,  method = "spearman")
  rows.cor <- cor(t(cor.matrix.R), method = "spearman")
  
  pheatmap(cor.matrix.R,
           cluster_cols = T,
           cluster_rows = T,
           clustering_distance_cols = as.dist(1 - cols.cor),
           clustering_distance_rows = as.dist(1 - rows.cor),
           clustering_method = "ward.D2",
           color =  my.colors, 
           breaks = my.breaks,
           legend_breaks = seq(-1,1,by=1),
           fontsize = 15,
           cellwidth = 26,
           cellheight = 26)
  
  dev.off()
}
#After obtaining the best subset of shRNAs, filter the coefficient matrix by removing other shRNAs. Then you could rank the drugs by different rules, such as top drugs with highest median coefficients.