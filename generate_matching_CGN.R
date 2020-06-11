# Code to generate unweighted gene coexpression networks from a correlation matrix that matches de edge density of a given network
# It loads the networks NSd(S) and NPd(P) and generates and saves the networks NPd(S) and NSd(P)

# It requieres 1 argument:
# 1) Name (including path) which will be used to load and save the files


#S is a triangular matrix
##A=matrix(0,nrow=nrow(S),ncol=nrow(S))
##A[S>=sort(as.vector(S),decreasing = T)[ecount(NPdP)]]=1
##NSdP=graph_from_adjacency_matrix(A,mode="upper")
##V(NSdP)$name=V(NSdS)$name

#P is a complete matrix, with diag=1
##diag(P)=NA
##A=matrix(0,nrow=nrow(S),ncol=nrow(S))
##A[P>=sort(as.vector(P),decreasing = T)[ecount(NSdS)*2]]=1
##NPdS=graph_from_adjacency_matrix(A,mode="upper")
##V(NPdS)$name=V(NSdS)$name
