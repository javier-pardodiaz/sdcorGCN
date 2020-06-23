# Code to analyse the networks (previously saved) using STRING
# It requieres 3 arguments:
# 1) Name (including path) which will be used to load the files (it needs to be the same as in generate_sdcor_GCN.R)
# 2) Path to STRING data with the proper format: 
#     A single .RData file with 3 matrices (not dataframes!). 
#     The names of the matrices must be "total", "coexpression" and "no_coexpression"
#     The three matrices are obtained using the combine_subscores.py script asnd following the indications in the Section 6 of the Supplementary information
#     Each matrix should have 3 columns: Interactor1, Interactor2, and Confidence score
#     We suggest you remove those entries (rows) which have genes that not appear in the expression data

require(Matrix)
require(igraph)


args = commandArgs(trailingOnly=TRUE)
name=args[1]
string_file=args[2]

load(string_file) #load file with the matrices with STRING information (proper formatting)
load(paste(name,"_NSdS.RData",sep=""))
load(paste(name,"_NSdP.RData",sep=""))
load(paste(name,"_NPdS.RData",sep=""))
load(paste(name,"_NPdP.RData",sep=""))
load(paste(name,"_NRdP.RData",sep=""))
load(paste(name,"_NRdP.RData",sep=""))


edges_NPdP=ends(NPdP,c(1:ecount(NPdP)))
edges_NPdS=ends(NPdS,c(1:ecount(NPdS)))
edges_NSdP=ends(NSdP,c(1:ecount(NSdP)))
edges_NSdS=ends(NSdS,c(1:ecount(NSdS)))
edges_NRdP=ends(NRdP,c(1:ecount(NRdP)))
edges_NRdS=ends(NRdS,c(1:ecount(NRdS)))

total_paste=cbind(paste(total[,1],total[,2]),total[,3])
coexpression_paste=cbind(paste(coex[,1],coex[,2]),coex[,3])
no_coexpression_paste=cbind(paste(no_coexpression[,1],no_coexpression[,2]),no_coexpression[,3])

edges_NPdP_paste=paste(edges_NPdP[,1],edges_NPdP[,2])
edges_NPdS_paste=paste(edges_NPdS[,1],edges_NPdS[,2])
edges_NSdS_paste=paste(edges_NSdS[,1],edges_NSdS[,2])
edges_NSdP_paste=paste(edges_NSdP[,1],edges_NSdP[,2])
edges_NRdS_paste=paste(edges_NRdS[,1],edges_NRdS[,2])
edges_NRdP_paste=paste(edges_NRdP[,1],edges_NRdP[,2])

sim_NPdP_total=total_paste[(total_paste[,1]%in%edges_NPdP_paste),]
sim_NSdP_total=total_paste[(total_paste[,1]%in%edges_NSdP_paste),]
sim_NPdS_total=total_paste[(total_paste[,1]%in%edges_NPdS_paste),]
sim_NSdS_total=total_paste[(total_paste[,1]%in%edges_NSdS_paste),]
sim_NRdS_total=total_paste[(total_paste[,1]%in%edges_NRdS_paste),]
sim_NRdP_total=total_paste[(total_paste[,1]%in%edges_NRdP_paste),]

sim_NPdP_coexpression=coexpression_paste[(coexpression_paste[,1]%in%edges_NPdP_paste),]
sim_NSdP_coexpression=coexpression_paste[(coexpression_paste[,1]%in%edges_NSdP_paste),]
sim_NPdS_coexpression=coexpression_paste[(coexpression_paste[,1]%in%edges_NPdS_paste),]
sim_NSdS_coexpression=coexpression_paste[(coexpression_paste[,1]%in%edges_NSdS_paste),]
sim_NRdS_coexpression=coexpression_paste[(coexpression_paste[,1]%in%edges_NRdS_paste),]
sim_NRdP_coexpression=coexpression_paste[(coexpression_paste[,1]%in%edges_NRdP_paste),]

sim_NPdP_no_coexpression=no_coexpression_paste[(no_coexpression_paste[,1]%in%edges_NPdP_paste),]
sim_NSdP_no_coexpression=no_coexpression_paste[(no_coexpression_paste[,1]%in%edges_NSdP_paste),]
sim_NPdS_no_coexpression=no_coexpression_paste[(no_coexpression_paste[,1]%in%edges_NPdS_paste),]
sim_NSdS_no_coexpression=no_coexpression_paste[(no_coexpression_paste[,1]%in%edges_NSdS_paste),]
sim_NRdS_no_coexpression=no_coexpression_paste[(no_coexpression_paste[,1]%in%edges_NRdS_paste),]
sim_NRdP_no_coexpression=no_coexpression_paste[(no_coexpression_paste[,1]%in%edges_NRdP_paste),]

# Generation of a full graph from which we sample random edges
A=matrix(1,nrow=vcount(NPdP),ncol=vcount(NPdP))
diag(A)=0
complete_graph=graph_from_adjacency_matrix(A,mode="upper")
V(complete_graph)$name=V(NSdP)$name
edges_all=ends(complete_graph,c(1:ecount(complete_graph)))

random_all_dS=list()
random_coex_dS=list()
random_no_coex_dS=list()
sum_values_all_dS=c()
sum_values_coex_dS=c()
sum_values_no_coex_dS=c()

for (i in c(1:30)){
  edges_random=edges_all[sample(c(1:ecount(complete_graph)),ecount(NSdS)),]
  edges_random_paste=paste(edges_random[,1],edges_random[,2])
  random_all_dS[[i]]=as.numeric(RL_total_paste[(RL_total_paste[,1]%in%edges_random_paste),2])
  random_coex_dS[[i]]=as.numeric(RL_coexpression_paste[(RL_coexpression_paste[,1]%in%edges_random_paste),2])
  random_no_coex_dS[[i]]=as.numeric(RL_no_coexpression_paste[(RL_no_coexpression_paste[,1]%in%edges_random_paste),2])
  sum_values_all_dS=c(sum_values_all_dS,sum(random_all_dS[[i]]))
  sum_values_coex_dS=c(sum_values_coex_dS,sum(random_coex_dS[[i]]))
  sum_values_no_coex_dS=c(sum_values_no_coex_dS,sum(random_no_coex_dS[[i]])) 
}

random_all_dP=list()
random_coex_dP=list()
random_no_coex_dP=list()
sum_values_all_dP=c()
sum_values_coex_dP=c()
sum_values_no_coex_dP=c()
for (i in c(1:30)){
  edges_random=edges_all[sample(c(1:ecount(complete_graph)),ecount(NSdP)),]
  edges_random_paste=paste(edges_random[,1],edges_random[,2])
  random_all_dP[[i]]=as.numeric(RL_total_paste[(RL_total_paste[,1]%in%edges_random_paste),2])
  random_coex_dP[[i]]=as.numeric(RL_coexpression_paste[(RL_coexpression_paste[,1]%in%edges_random_paste),2])
  random_no_coex_dP[[i]]=as.numeric(RL_no_coexpression_paste[(RL_no_coexpression_paste[,1]%in%edges_random_paste),2])
  sum_values_all_dP=c(sum_values_all_dP,sum(random_all_dP[[i]]))
  sum_values_coex_dP=c(sum_values_coex_dP,sum(random_coex_dP[[i]]))
  sum_values_no_coex_dP=c(sum_values_no_coex_dP,sum(random_no_coex_dP[[i]]))
}

#Table with results

results_STRING_table=matrix(nrow=8,ncol=3)
rownames(results_STRING_table)=c("NSdS","NPdS","NRdS","REdS","NSdP","NPdP","NRdP","REdP")
colnames(results_STRING_table)=c("All","Coexpression","No coexpression")
results_STRING_table[1,]=c(sum(as.numeric(sim_NSdS_total[,2])),sum(as.numeric(sim_NSdS_coexpression[,2])),sum(as.numeric(sim_NSdS_no_coexpression[,2])))
results_STRING_table[2,]=c(sum(as.numeric(sim_NPdS_total[,2])),sum(as.numeric(sim_NPdS_coexpression[,2])),sum(as.numeric(sim_NPdS_no_coexpression[,2])))
results_STRING_table[3,]=c(sum(as.numeric(sim_NRdS_total[,2])),sum(as.numeric(sim_NRdS_coexpression[,2])),sum(as.numeric(sim_NRdS_no_coexpression[,2])))
results_STRING_table[4,]=c(mean(sum_values_all_dS),mean(sum_values_coex_dS),mean(sum_values_no_coex_dS))
results_STRING_table[5,]=c(sum(as.numeric(sim_NSdP_total[,2])),sum(as.numeric(sim_NSdP_coexpression[,2])),sum(as.numeric(sim_NSdP_no_coexpression[,2])))
results_STRING_table[6,]=c(sum(as.numeric(sim_NPdP_total[,2])),sum(as.numeric(sim_NPdP_coexpression[,2])),sum(as.numeric(sim_NPdP_no_coexpression[,2])))
results_STRING_table[7,]=c(sum(as.numeric(sim_NRdP_total[,2])),sum(as.numeric(sim_NRdP_coexpression[,2])),sum(as.numeric(sim_NRdP_no_coexpression[,2])))
results_STRING_table[8,]=c(mean(sum_values_all_dP),mean(sum_values_coex_dP),mean(sum_values_no_coex_dP))

print(results_STRING_table)
save(results_STRING_table,file=paste(name,"_STRING_results.RData",sep=""))

#Below the code to plot the results as in the paper
par(mfrow=c(1,4))
boxplot(x=list(sum_values_all_dS,sum_values_all_dP),
        main="All evidence",
        ylab="Score",
        xlab="Network edge density",
        names=c("d(S)","d(P)"),
        ylim=c(0,sum(as.numeric(sim_NSdP_total[,2])+10)),
        yaxs="i",cex.lab=1.5, cex.axis=1.5,cex.main=1.5
)

points(1,sum(as.numeric(sim_NPdS_total[,2])),bg="red",pch=24,cex=2,col="black")
points(2,sum(as.numeric(sim_NPdP_total[,2])),bg="red",pch=24,cex=2,col="black")
points(1,sum(as.numeric(sim_NRdS_total[,2])),bg="yellow",pch=22,cex=2,col="black")
points(2,sum(as.numeric(sim_NRdP_total[,2])),bg="yellow",pch=22,cex=2,col="black")
points(1,sum(as.numeric(sim_NSdS_total[,2])),bg="blue",pch=21,cex=2,col="black")
points(2,sum(as.numeric(sim_NSdP_total[,2])),bg="blue",pch=21,cex=2,col="black")

boxplot(x=list(sum_values_coex_dS,sum_values_coex_dP),
        main="Only coexpression evidence",
        ylab="Score",
        xlab="Network edge density",
        names=c("d(S)","d(P)"),
        ylim=c(0,sum(as.numeric(sim_NSdP_coexpression[,2])+10)),
        yaxs="i",cex.lab=1.5, cex.axis=1.5,cex.main=1.5
)

points(1,sum(as.numeric(sim_NPdS_coexpression[,2])),bg="red",pch=24,cex=2,col="black")
points(2,sum(as.numeric(sim_NPdP_coexpression[,2])),bg="red",pch=24,cex=2,col="black")
points(1,sum(as.numeric(sim_NRdS_coexpression[,2])),bg="yellow",pch=22,cex=2,col="black")
points(2,sum(as.numeric(sim_NRdP_coexpression[,2])),bg="yellow",pch=22,cex=2,col="black")
points(1,sum(as.numeric(sim_NSdS_coexpression[,2])),bg="blue",pch=21,cex=2,col="black")
points(2,sum(as.numeric(sim_NSdP_coexpression[,2])),bg="blue",pch=21,cex=2,col="black")



boxplot(x=list(sum_values_no_coex_dS,sum_values_no_coex_dP),
        main="No coexpression evidence",
        ylab="Score",
        xlab="Network edge density",
        names=c("d(S)","d(P)"),
        ylim=c(0,sum(as.numeric(sim_NSdP_no_coexpression[,2])+10)),
        yaxs="i",cex.lab=1.5, cex.axis=1.5,cex.main=1.5
)

points(1,sum(as.numeric(sim_NPdS_no_coexpression[,2])),bg="red",pch=24,cex=2,col="black")
points(2,sum(as.numeric(sim_NPdP_no_coexpression[,2])),bg="red",pch=24,cex=2,col="black")
points(1,sum(as.numeric(sim_NRdS_no_coexpression[,2])),bg="yellow",pch=22,cex=2,col="black")
points(2,sum(as.numeric(sim_NRdP_no_coexpression[,2])),bg="yellow",pch=22,cex=2,col="black")
points(2,sum(as.numeric(sim_NSdP_no_coexpression[,2])),bg="blue",pch=21,cex=2,col="black")
points(1,sum(as.numeric(sim_NSdS_no_coexpression[,2])),bg="blue",pch=21,cex=2,col="black")

plot(10, 10, axes=FALSE, frame.plot=FALSE,ylim=c(-1,1),xlim=c(-1,1),ylab = "",xlab="",xaxs="i")
legend(x=-1, y=0.5,c("Signed distance cor.","Pearson cor.","Spearman cor."),pch=c(21,24,22),pt.bg=c("blue","red","yellow"),cex=1.5)


