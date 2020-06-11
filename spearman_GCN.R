# Code to generate the networks NRdS and NRdP from a correlation matrix R obtained using Spearman correlation
# It loads the networks NSdS and NSdP

# It requieres 2 arguments:
# 1) Input file .RData with the expression matrix M (rows=genes, columns=measurements)
# 2) Name (including path) which will be used to load and save files

require("devtools")
install_github("lbozhilova/COGENT")
require(COGENT)
require(preprocessCore)
require(parallel)
require(Matrix)
require(igraph)

args = commandArgs(trailingOnly=TRUE)
file=args[1]
name=args[2]
load(file) #load file with the expression matrix. The expression matrix needs to be called "M".
load(paste(name,"_NPdP.RData",sep=""))
load(paste(name,"_NSdS.RData",sep=""))


preprocessing=function(em){
  # Sets the 20% of the lowest expressed values from each sample to the lowest value in the original dataset, and quantile-normalises the data
  # It first calculates the minumum value in the original dataset and identifies which entries need to be modified
  # Arguments: em (expression matrix where rows are genes and columns different samples)
  min_val=min(em,na.rm = T)
  p=0.2
  genes_remove=apply(em,2,function(x) get_lowest_values(x,p)
                     )
  # It then removes those genes from the dataset and performs the quantile normalisation
  em[which(genes_remove==T)]=NA
  em_qnorm=normalize.quantiles(em)
  # Lastly, it sets the identified values to the minimum gene expression value in the original dataset
  em_qnorm[which(genes_remove==T)]=min_val
  return(em_qnorm)
}

get_lowest_values=function(column,p){
  # Identification of the p lowest values from each sample
  # Arguments: column (vector to be tested), p (proportion of values to be identified)
  v=column<=quantile(column,p,na.rm=T,type=1)
  return(v)
}

srcor_matrix=function(em,prep=T,threads=1){
  # Preprocesses expression matrix (if prep=T), and computes and returns a Spearman correlation matrix
  # Arguments: em (expression matrix in which columns are samples and rows genes), prep (T if preprocessing is needed), threads (number of threads to be used)
 
  em=as.matrix(em[,-1])
  if (prep == T){
    # PREPROCESSING
    M_q=preprocessing(em)
    rownames(M_q)=rownames(em)
  }else{
    M_q=em
  }
  rm(em)
  # CORRELATION MEASUREMENT
  # Calculate Spearman correlation matrix
  R=cor(t(M_q),use = "pairwise.complete.obs",method="spearman")
  rownames(R)=rownames(M_q)
  colnames(R)=rownames(M_q)
  rm(M_q)
  return(R)
}

df=as.data.frame(M/10) #Some functions in COGENT only work with "double" values, we need to transform possible "numeric" values
df=cbind(rownames(M),df) 
colnames(df)[1]="Name"

R=srcor_matrix(df,prep= T,threads =  threads)
diag(R)=NA
A=matrix(0,nrow=nrow(M),ncol=nrow(M))
A[R>=sort(as.vector(R),decreasing = T)[ecount(NSdS)*2]]=1
NRdS=graph_from_adjacency_matrix(A,mode="upper")
V(NRdS)$name=V(NSdS)$name

save(NRdS, file=paste(name,"_NRdS.RData",sep=""))


A=matrix(0,nrow=nrow(M),ncol=nrow(M))
A[R>=sort(as.vector(R),decreasing = T)[ecount(NPdP)*2]]=1
NRdS=graph_from_adjacency_matrix(A,mode="upper")
V(NRdP)$name=V(NPdP)$name

save(NRdP, file=paste(name,"_NRdP.RData",sep=""))
