##Code to generate weighted and thresholded networks from correlation matrices
##Previous to run this code, you will need to run the code on https://github.com/javier-pardodiaz/sdcorGCN/blob/master/generate_sdcor_GCN.R lines 1-321 and 329-333

##The variable "name" should match the one used previously so that the files with the adjacency matrices for the different iterations can be accessed

iterations=25
##Threshold values to test

thrs=c(0.1,0.3,0.4,0.45,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.65,0.7,0.9)

## n is the number of genes in the dataset
n=nrow(M)


library(igraph)
##Similarity matrix that contains the similarity values for each of the 25 iterations and all tested threshold values

sim=matrix(nrow=interations,ncol=length(thrs))
for (i in c(1:iterations)){
  load(paste(name,"_cm_",iter,"_A.RData",sep = ""))
  # load(paste(name,"p_cm_",iter,"_A.RData",sep = "")) # To use Pearson correlation instead of signed distance correlation
  cm_A=CM
  load(paste(name,"_cm_",iter,"_B.RData",sep = ""))
  # load(paste(name,"p_cm_",iter,"_B.RData",sep = "")) # To use Pearson correlation instead of signed distance correlation
  cm_B=CM
  rm(CM)
  cm_A[is.na(cm_A)]=0
  cm_B[is.na(cm_B)]=0
  diag(cm_B)=0
  diag(cm_A)=0
  
  # The correlation matrix is/might be a triangular matrix and we need to correct for this
  if(sum(as.vector(cm_A[lower.tri(cm_A, diag = FALSE)]),na.rm = T)==0){
    cm_A=cm_A+t(cm_A)
  }
  if(sum(as.vector(cm_B[lower.tri(cm_B, diag = FALSE)]),na.rm = T)==0){
    cm_B=cm_B+t(cm_B)
  }
  
  for (thr in thrs){
    
    ##Generate random networks with the same weight distribution as cm_A and cm_B
    
    rm_A_v=sample(c(1:ncol(cm_A)),ncol(cm_A))
    rm_A=cm_A[c(rm_A_v),]
    rm_A=rm_A[,rm_A_v]
    rm_B_v=sample(c(1:ncol(cm_B)),ncol(cm_B))
    rm_B=cm_B[c(rm_B_v),]
    rm_B=rm_B[,rm_B_v]
    
    ##Threshold correlation matrices
  
    cm_1=cm_A
    cm_1[cm_A<thr]=0
    cm_2=cm_B
    cm_2[cm_B<thr]=0
  
    rm_A[rm_A<thr]=0
    rm_B[rm_B<thr]=0
    
    ## Calulate the correction term gamma
  
    int1=sum(rm_B[cm_1>=rm_B],cm_1[rm_B>cm_1])
    int2=sum(rm_A[cm_2>=rm_A],cm_2[rm_A>cm_2])
    cor=(int1+int2)/2
    
    ## Calculate similarity
    
    den=sum(cm_1[cm_1>=cm_2],cm_2[cm_2>cm_1])+cor
    num=sum(cm_2[cm_1>=cm_2],cm_1[cm_2>cm_1])-cor
    
    sim[i,match(thr,thrs)]=(num/den)
  }
}

## Calculate mean similarity for the 25 iterations

sim_mean=apply(sim,2,mean,na.rm=T)

## Calculate the sum of weights for each threshold

sum_weights=c()

for (thr in thrs){
  A_S=S
  A_S[A_S<thr]=0
  diag(A_S)=0
  
  sum_weights=c(sum_weights,sum(A_S))
}

## Calculate the score for each threshold

score=(sim_mean)-(sum_weights/(n^2-n))


##Edit ylim and xlim so that it fits your data

plot(x=sum_weights/2,y=score,type="l",col="blue",lwd=3,ylim=c(0.3,0.5),xlim=c(0,1500000),main="R.")
abline(v=sum_weights[match(max(score),score)]/2, col="blue",lty=2)

abline(h=max(score),col="blue",lty=2)


##The optimal threshold results in the highest score

thr_s=thrs[match(max(score),score)]

w_s=sum_weights[match(max(score),score)]

##Adjancency matricx for the full dataset.
## Need to load the correlation matrix for the full dataset (obtained from https://github.com/javier-pardodiaz/sdcorGCN/blob/master/generate_sdcor_GCN.R lines 329-333 )
load(paste(name,"_cm.RData",sep=""))
A_S=S
A_S[A_S<thr_s]=0
diag(A_S)=0
## Obtain network and save
WSwS=graph_from_adjacency_matrix(A_S,mode="upper",weighted = T)
save(WSwS,file=paste(name,"_NSwS.RData",sep=""))


## To generate the optimal Pearson correlation network (you need to load the Pearson correlation matrices in lines 21 and 23)
## Need to load the correlation matrix for the full dataset (obtained from https://github.com/javier-pardodiaz/sdcorGCN/blob/master/generate_sdcor_GCN.R lines 329-333 )
load(paste(name,"p_cm.RData",sep=""))
A_P=P
A_P[A_P<thr_s]=0
diag(A_P)=0
## Obtain network and save
WPwP=graph_from_adjacency_matrix(A_P,mode="upper",weighted = T)
save(WPwP,file=paste(name,"_NPwP.RData",sep=""))
