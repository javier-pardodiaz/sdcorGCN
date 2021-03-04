## Code to generate an unweighted network NSwS with the same edges as the optimal weighted network WSwS
## And a weighted network WSdS with the same edges as the optimal unweighted network NSdS

load(paste(name,"_WSwS.RData",sep=""))
load(paste(name,"_NSdS.RData",sep=""))
load(paste(name,"_cm.RData",sep=""))

A=S
thr=sort(S,decreasing = T)[2*ecount(NSdS)]
A[A<thr]=0
diag(A)=0
WSdS=graph_from_adjacency_matrix(A,mode="upper",weighted = T)

NSwS=WSwS
E(NSwS)$weight=mean(E(WSwS)$weight)

E(NSdS)$weight=mean(E(WSdS)$weight)

