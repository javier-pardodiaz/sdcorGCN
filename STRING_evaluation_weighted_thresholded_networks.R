# Code to analyse the networks (previously saved) using STRING
# It requieres 3 arguments:
# 1) Name (including path) which will be used to load the files (it needs to be the same as in generate_sdcor_GCN.R)
# 2) Path to STRING data with the proper format: 
#     A single .RData file with 3 matrices (not dataframes!). 
#     The names of the matrices must be "total", "coexpression" and "no_coexpression"
#     The three matrices are obtained using the combine_subscores.py script asnd following the indications in the Section 6 of the Supplementary information
#     Each matrix should have 3 columns: Interactor1, Interactor2, and Confidence score
#     We suggest you remove those entries (rows) which have genes that not appear in the expression data

# Here we show the STRING evaluation for the network WSdS using all the information from STRING. The evaluation for the other networks and sets of data is similar.
# We also show how we generate 30 random networks with edge weight wS
require(Matrix)
require(igraph)


args = commandArgs(trailingOnly=TRUE)
name=args[1]
string_file=args[2]

args = commandArgs(trailingOnly=TRUE)
name=args[1]
string_file=args[2]

load(string_file) #load file with the matrices with STRING information (proper formatting)
load(paste(name,"_WSwS.RData",sep=""))

edges_WSwS=ends(WSwS,c(1:ecount(WSwS)))

total_paste=cbind(paste(total[,1],total[,2]),total[,3])
coexpression_paste=cbind(paste(coex[,1],coex[,2]),coex[,3])
no_coexpression_paste=cbind(paste(no_coexpression[,1],no_coexpression[,2]),no_coexpression[,3])

edges_WSwS_paste=cbind(paste(edges_WSwS[,1],edges_WSwS[,2]),E(WSwS)$weight)
colnames(edges_WSwS_paste)=c("edge","weight")
sim_WSwS_total=total_paste[(total_paste[,1]%in%edges_WSwS_paste[,1]),]
colnames(sim_WSwS_total)=c("edge","score")
edge_score_weight_WSwS_total=as.matrix(merge(sim_WSwS_total,edges_WSwS_paste))
score_WSwS_total=sum(as.numeric(edge_score_weight_WSwS_total[,2])*as.numeric(edge_score_weight_WSwS_total[,3]))/sum(E(WSwS)$weight)


#Random networks
a=combn(V(WSwS)$name, 2)

random_total_wS=list()

weights=E(WSwS)$weight

for (i in c(1:30)){
  
  edges_random_index=sample(ncol(a),ecount(WSwS))
  edges_paste=paste(a[1,edges_random_index],a[2,edges_random_index])
  edges_random_paste=cbind(edges_paste,weights)
  colnames(edges_random_paste)=c("edge","weight")
  
  sim_WRwS_total=total_paste[(total_paste[,1]%in%edges_random_paste[,1]),]
  colnames(sim_WRwS_total)=c("edge","score")
  edge_score_weight_WRwS_total=as.matrix(merge(sim_WRwS_total,edges_random_paste))
  random_total_wS[[i]]=sum(as.numeric(edge_score_weight_WRwS_total[,2])*as.numeric(edge_score_weight_WRwS_total[,3]))/sum(weights)
  
}
