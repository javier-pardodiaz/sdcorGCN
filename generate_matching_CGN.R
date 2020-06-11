# Code to generate unweighted gene coexpression networks from a correlation matrix that matches de edge density of a given network
# It loads the networks NSd(S) and NPd(P) and generates and saves the networks NPd(S) and NSd(P)
# It shows on screen the summaries of the networks

# It requieres 1 argument:
# 1) Name (including path) which will be used to load and save the files. It needs to be the same as in generate_sdcorGCN.R

require("devtools")
install_github("lbozhilova/COGENT")
require(COGENT)
require(preprocessCore)
require(parallel)
require(Matrix)
require(igraph)

args = commandArgs(trailingOnly=TRUE)
name=args[1]

load(paste(name,"_NPdP.RData",sep=""))
load(paste(name,"_NSdS.RData",sep=""))
load(paste(name,"_p_cm.RData",sep=""))
load(paste(name,"_cm.RData",sep=""))


#S is a triangular matrix
A=matrix(0,nrow=nrow(S),ncol=nrow(S))
A[S>=sort(as.vector(S),decreasing = T)[ecount(NPdP)]]=1
NSdP=graph_from_adjacency_matrix(A,mode="upper")
V(NSdP)$name=V(NSdS)$name

#P is a complete matrix, with diag=1
diag(P)=NA
A=matrix(0,nrow=nrow(S),ncol=nrow(S))
A[P>=sort(as.vector(P),decreasing = T)[ecount(NSdS)*2]]=1
NPdS=graph_from_adjacency_matrix(A,mode="upper")
V(NPdS)$name=V(NSdS)$name

# Print networks summaries

network_summaries=function(network,name){
  edges=ecount(network)
  components_network=components(network)
  lcc=induced_subgraph(network,components_network$membership==names(sort(table(components_network$membership),decreasing=TRUE)[1]))
  lcc_size=vcount(lcc)
  lcc_density=edge_density(lcc)*100
  lcc_gcc=transitivity(lcc)
  cat(name,"\nNumber Edges: ", edges,"\nSize Lcc: ",lcc_size,"\nEdge density LCC*100: ",lcc_density,"\nGlobal clustering coefficient LCC: ", lcc_gcc,"\n\n")
  return(list(edges,lcc_size,lcc_density,lcc_gcc))
}

networks=c("NSdS","NPdS","NSdP","NPdP")
summaries=lapply(networks,function(networks) network_summaries(eval(parse(text = networks)),networks))
