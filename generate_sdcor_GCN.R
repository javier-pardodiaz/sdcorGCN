# Code to generate an unweighted gene coexpression network using signed distance correlation
# It requieres 3 arguments:
# 1) Input file .RData with the expression matrix M (rows=genes, columns=measurements)
# 2) Name (including path) which will be used to save the generated files
# 3) Number of threads to be used
# There is a commented out section which generates unweighted gene coexpression networks using Pearson and Spearman correlations. This section has been used in the paper in order to compare the network construction methods. 

# This code saves two correlation matrices (A and B) for each COGENT iteration (see paper). Afterward, it loads the correlation matrices to generate and compare gene coexpression networks using different thresholds
# This code saves the unweighted gene coexpression network (igraph object) obtained using distance correlation and the optimal threshold value according to COGENT

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
threads=as.numeric(args[3])
load(file) #load file with the expression matrix. The expression matrix needs to be called "M".

iterations=25 #Number of COGENT iterations

# Threshold values to be tested

thrs=c(0.3,0.4,0.45,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.8,0.9,0.95) #Rhizobium leguminosarum dataset
#thrs=c(0.3,0.5,0.6,0.65,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95) #Yeast Dataset


# FUNCTIONS
# PREPROCESSING FUNCTIONS

preprocessing=function(em){
  # Sets the 20% of the lowest expressed values from each sample to the lowest value in the original dataset, and quantile-normalises the data
  # It first calculates the minumum value in the original dataset and identifies which entries need to be modified
  # Arguments: em (expression matrix where rows are genes and columns different samples)
  min_val=min(em,na.rm = T)
  p=0.2
  genes_remove=apply(em,2,function(x) get_lowest_values(x,p))
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

# CORRELATION MEASUREMENT FUNCTIONS

dif=function(a,b){
  # Function employed in the calculation of signed distance correlation
  # Calculates difference between two numbers
  return(abs(a-b))
}

element_minus_vector_elements=function(element,vector){
  # Function employed in the calculation of signed distance correlation
  # For each the element in the vector, it computes the absolute value of its difference with the resto of the elements
  # It returns a vector with as many zeros at the start as the position of the element in the vector since those values will have been calculated already, in a previous iteration
  j=c(element+1:(length(vector)-element))
  distance_vector=sapply(j,function(j) dif(vector[element],vector[j]))
  distance_vector=c(rep(0,element),distance_vector)
  return(distance_vector)
}

center_matrix_1=function(row_idx,m_distance_2,mean_elements,mean_total){
  # Function employed in the calculation of signed distance correlation
  # It double-centers the expression distances
  j=c(1:row_idx)
  row_norm=sapply(j, function(j) center_matrix_2(row_idx,j,m_distance_2,mean_elements,mean_total))
  row_norm=c(row_norm,rep(0,length(mean_elements)-row_idx))
  return(row_norm)
}

center_matrix_2=function(row_idx,col_idx, m_distance_2,mean_elements,mean_total){
  # Function employed in the calculation of signed distance correlation
  # It double-centers the expression distances
  value_norm=m_distance_2[row_idx,col_idx]-mean_elements[row_idx]-mean_elements[col_idx]+mean_total
  return(value_norm)
}

get_normalised_distances_matrix=function(expression_vector){
  # Function that from a expression vector, with all the sample values for a single gene, returns its double centered distances matrix
  # Arguments: expression_vector (expression values of a gene in different samples)
  i = c(1:(length(expression_vector)-1))
  # Get a triangular distances matrix (absolute values)
  m_distance=sapply(i,function(i) element_minus_vector_elements(i,expression_vector))
  m_distance=cbind(m_distance,0)
  # Get the complete matrix
  m_distance_2=m_distance+t(m_distance)
  mean_elements=apply(m_distance_2,1,mean)
  mean_total=mean(m_distance_2)
  k = c(1:(length(expression_vector)))
  # Double centering
  m_norm=sapply(k,function(k) center_matrix_1(k,m_distance_2,mean_elements,mean_total))
  # Return the matrix as a vector
  m_norm_2=as.vector(m_norm+t(m_norm)-(m_norm*Diagonal(nrow(m_norm))))
  return(m_norm_2)
}

pairwise_correlation=function(indexes){
  # Calculates the Pearson correlation of two distance matrices (as columns in a distances_matrix matrix)
  # Arguments: indexes (columns in the distances_matrix to be used)
  corr_v=(cor(distances_matrix[,indexes[1]],distances_matrix[,indexes[2]],method="pearson"))
  return(corr_v)
}

distance_correlation_simple=function(em,cores=1){
  # For a (quantile normalised) expression matrix, it computes and returns a distance correlation matrix
  # Calculates the doubled centered matrix of distances for each gene. This operation can be done using one or more cores
  # Arguments: em (expression matrix), cores (threads to be used in the computing)
  if (cores>1){
    cl1= makeCluster(cores, "FORK")
    distances_matrix<<-parApply(cl1,em,1,get_normalised_distances_matrix)
    stopCluster(cl1)
  }else{
    distances_matrix<<-apply(em,1,get_normalised_distances_matrix)
  }
  
  # Get all the possible combinations of pairs of genes (each pair of genes is in one of the columns in pairwise_combs)
  pairwise_combs<<-combn(nrow(em),2)
  l=ncol(pairwise_combs)

  # Calculate the correlation of the distances vectors for each pairs of genes. This operation can be done using one or more cores.
  if (cores>1){
    cl2= makeCluster(cores, "FORK")
    correlation_vector_d=parApply(cl2,pairwise_combs, 2, pairwise_correlation)
    stopCluster(cl2)
    
  }else{
    correlation_vector_d=apply(pairwise_combs, 2, pairwise_correlation)
  }
  
  # Square root of the correlation vector to get distance correlation
  correlation_vector_d=sqrt(correlation_vector_d)
  
  # Transform vector into a triangular matrix
  cm=matrix(0,nrow=nrow(em),ncol=nrow(em))
  cm[lower.tri(cm, diag=FALSE)] <- correlation_vector_d
  return(t(cm))
}

sdcor_matrix=function(em,prep=T,threads=1){
  # Preprocesses expression matrix (if prep=T), and computes and returns a signed distance correlation matrix
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
  # Calculate Pearson and Distance correlation matrices
  P=cor(t(M_q),use = "pairwise.complete.obs",method="pearson")
  D=distance_correlation_simple(M_q,threads) #note that D is a triangular matrix
  # Assign the sign of the entries in P to the values in D
  S=D*sign(P) #note that S is a triangular matrix
  rownames(S)=rownames(M_q)
  colnames(S)=rownames(M_q)
  rm(P,D,M_q)
  return(S)
}

pcor_matrix=function(em,prep=T,threads=1){
  # Preprocesses expression matrix (if prep=T), and computes and returns a Pearson correlation matrix
  # Arguments: em (expression matrix in which columns are samples and rows genes), prep (T if preprocessing is needed), threads (number of threads to be used)
  em=as.matrix(em[,-1])
  if (prep == T){
    #PREPROCESSING
    M_q=preprocessing(em)
    rownames(M_q)=rownames(em)
  }else{
    M_q=em
  }
  rm(em)
  # CORRELATION MEASUREMENT
  # Calculate Pearson correlation matrix
  P=cor(t(M_q),use = "pairwise.complete.obs",method="pearson")
  rownames(P)=rownames(M_q)
  colnames(P)=rownames(M_q)
  rm(M_q)
  return(P)
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

# THRESHOLDING FUNCTIONS

cogentParallelCM=function(df,corrFun,repCount=25,propShared=0.5,threadCount = 4,name,preprocessing){
  # It generates and saves repCount pairs of correlation matrices
  # Arguments: df (dataframe containing the gene expression data), corrFun (function employed to generate the network), repCount (number of iterations), propShared (proportion of the total number of samples shared between the two networks generated in each iteration),threadCount (numbers of threads to be used), name (used to save the generated files), preprocessing (whether preprocessing is necessary) 
  j=c(1:repCount)
  cl1= makeCluster(threadCount, "FORK")
  parSapply(cl1,j, function(j) cogentSingleCM(df,corrFun,propShared,name,iter=j,preprocessing))
  stopCluster(cl1)
  return(0)
}

cogentSingleCM=function(df,corrFun, propShared,name,iter,preprocessing){
  # It splits the expression data in two overlaping sets for which it computes and saves their correlation matrices
  # Arguments: df (dataframe containing the gene expression data), corrFun (function employed to generate the network), propShared (proportion of the total number of samples shared between the two networks generated in each iteration),threadCount (numbers of threads to be used), name (used to save the generated files), iter (iteration number), preprocessing (whether preprocessing is necessary) 

  set.seed(iter)
  dfList <- splitExpressionData(df, propShared)
  CM_reps <- lapply(dfList, function(x) do.call(corrFun, list(x),preprocessing))
  # Note that CM_reps contains triangular matrices
  CM=CM_reps[[1]]
  save(CM,file=paste(name,"_cm_",iter,"_A.RData",sep = ""))
  CM=CM_reps[[2]]
  save(CM,file=paste(name,"_cm_",iter,"_B.RData",sep = ""))
  rm(CM_reps)
  return(0)
}

threshold_cm=function(correlation_matrix, thr){
  # It sets all values in the correlation matrix which are lower than thr to 0
  # Arguments: correlation_matrix (square matrix with correlation values for all pairs of genes), thr (threshold)
  adj=Matrix(0L,ncol=ncol(correlation_matrix),nrow=ncol(correlation_matrix),sparse = T)
  adj[correlation_matrix>=thr]=1
  adj[correlation_matrix<thr]=0
  # The correlation matrix is/might be a triangular matrix and we need to correct for this
  if(sum(as.vector(correlation_matrix[lower.tri(correlation_matrix, diag = FALSE)]),na.rm = T)==0){
    adj=adj+t(adj)
  }
  rm(correlation_matrix)
  diag(adj)=0L
  return(adj)
}

get_unweighted_network_adjancency_matrix_and_compare=function(thr,cm_1,cm_2){
  # Generate and compare unweighted networks. It returns a similarity measure obtained using the COGENT package
  # Arguments: thr (threshold value), cm_1 and cm_2 (correlation matrices)
  A_list <- lapply(list(cm_1,cm_2), function(x) do.call(threshold_cm, list(x, thr)))
  rm(cm_1,cm_2)
  similarity_scores=getEdgeSimilarityCorrected(A_list,reduce=TRUE,type="expected")[[2]]
  rm(A_list)
  return(similarity_scores)
}

cogentSingleThr=function(thrs,name,iter,align=FALSE){
  # Load 2 correlation matrices, generate and compare unweighted networks obtained using different thresholds
  # Arguments: thrs (threshold values to test), name (file name used to load files), iter (iteration number)
  load(paste(name,"_cm_",iter,"_A.RData",sep = ""))
  CM[is.na(CM)]=0
  cm_1=Matrix(CM,sparse = T)
  load(paste(name,"_cm_",iter,"_B.RData",sep = ""))
  CM[is.na(CM)]=0
  cm_2=Matrix(CM,sparse = T)
  rm(CM)
  similarity_scores=lapply(thrs,function(thrs) get_unweighted_network_adjancency_matrix_and_compare(thrs, cm_1,cm_2))
  rm(cm_1,cm_2)
  return(similarity_scores)
}

cogentParallelThr=function(repCount=25,threadCount = 4,name,thrs){
  # Generates and compares unweighted networks obtained different thresholds and correlation matrices previously saved
  # Arguments: repCount (number of iterations), threadCount (number of threads to be used), name (name used to load and save files), thrs (thresholds to test)
  j = c(1:repCount)
  cl3= makeCluster(threadCount, "FORK")
  similarity=parSapply(cl3,j, function(j) cogentSingleThr(thrs,name,iter=j))
  stopCluster(cl3)
  return(similarity)
}

# END OF FUNCTIONS

# Assign the appropriate format to the correlation matrix
df=as.data.frame(M/10) #Some functions in COGENT only work with "double" values, we need to transform possible "numeric" values
df=cbind(rownames(M),df) 
colnames(df)[1]="Name"

# Generate 25 (iterations) pairs of signed distance correlation matrices from overlapping datasets
cogentParallelCM(df,"sdcor_matrix",repCount = iterations,propShared =  0.5, threadCount = 4,name = name ,preprocessing =  TRUE)

# Generate and compare 25 (iterations) pairs of unweighted networks using different threshold values 
similarity_scores=cogentParallelThr(repCount = iterations, threadCount = threads,name =name,thrs=thrs)
save(similarity_scores,file=paste(name,"_similarity_results_thrs.RData"))
similarity_scores=matrix(as.vector(similarity_scores,mode="numeric"),nrow=length(thrs))
similarities=apply(similarity_scores,1,mean) # Vector with a similarity value for each tested threshold

# Generate signed distance correlation matrix (with all the samples in the dataset)
S=sdcor_matrix(df,prep= T,threads =  threads)
S=S+t(S)
diag(S)=NA
save(S,file=paste(name,"_cm.RData",sep=""))

densities=c()
for (thr in thrs){
  densities=c(densities,sum(S>=thr,na.rm = T)/(nrow(S)*(nrow(S)-1))) # Vector with an edge density value for each tested threshold
}
save(densities,file=paste(name,"_densities.RData",sep=""))
scores=similarities-densities # Vector with a score value for each tested threshold
save(scores,file=paste(name,"_scores.RData",sep=""))

# Unweighted signed distace correlation network with the optimal threshold for the signed distance correlation matrix
A=matrix(0,nrow=nrow(S),ncol=nrow(S)) # Adjancency matrix
A[S>=thrs[match(max(scores),scores)]]=1
NSdS=graph_from_adjacency_matrix(A,mode="upper") # Network
V(NSdS)$name=rownames(M)

save(NSdS, file=paste(name,"_NSdS.RData",sep=""))




