#Code to generate weighted and thresholded networks that match the optimal networks for Pearson and signed distance correlation
load(name,"_WPwP.RData",sep="")
load(name,"_WSwS.RData",sep="")

load(paste(name,"_cm.RData",sep=""))
load(paste(name,"p_cm.RData",sep=""))

#Networks with matching sums of edge weights
limits=sort(thrs_p[match(sort(abs(sum_weights_p-w_s),decreasing = F)[c(1,2)],abs(sum_weights_p-w_s))])

new_w=sum(P[P>=limits[1]],na.rm = T)
new_limit=sum(limits)/2
while (abs(new_w-w_s)>new_limit){
  new_limit=sum(limits)/2
  print(new_limit)
  new_w=sum(P[P>=new_limit],na.rm = T)
  print(new_w)
  if (new_w>w_s){
    limits[1]=new_limit
  }else{
    limits[2]=new_limit
  }
  print(limits)
}

A_P=P
A_P[A_P<new_limit]=0
diag(A_P)=0
WPwS=graph_from_adjacency_matrix(A_P,mode="upper",weighted = T)

save(WPwS,file=paste(name,"_WPwS.RData",sep=""))


limits=sort(thrs[match(sort(abs(sum_weights_s-w_p),decreasing = F)[c(1,2)],abs(sum_weights_s-w_p))])

new_w=sum(S[S>=limits[1]],na.rm = T)
new_limit=sum(limits)/2

while (abs(new_w-w_p)>new_limit){
  new_limit=sum(limits)/2
  print(new_limit)
  new_w=sum(S[S>=new_limit],na.rm = T)
  print(new_w)
  if (new_w>w_p){
    limits[1]=new_limit
  }else{
    limits[2]=new_limit
  }
  print(limits)
}

A_S=S
A_S[A_S<new_limit]=0
diag(A_S)=0
WSwP=graph_from_adjacency_matrix(A_S,mode="upper",weighted = T)

save(WSwP,file=paste(name,"_WSwP.RData",sep=""))


#Networks with matching mean edge weights

limits=c(0.7,0.5)

m_s=mean(E(WSwS)$weight)
new_m=mean(P[P>=limits[1]],na.rm = T)
new_limit=sum(limits)/2
while (abs(new_m-m_s)>(new_limit/sum(P>new_limit,na.rm = T))){
  new_limit=sum(limits)/2
  print(new_limit)
  new_m=mean(P[P>=new_limit],na.rm = T)
  print(new_m)
  if (new_m>m_s){
    limits[1]=new_limit
  }else{
    limits[2]=new_limit
  }
  print(limits)
}

A_P=P
A_P[A_P<new_limit]=0
diag(A_P)=0
WPaS=graph_from_adjacency_matrix(A_P,mode="upper",weighted = T)

save(WPaS,file=paste(name,"_WPaS.RData",sep=""))


limits=c(0.7,0.5)

m_p=mean(E(WPwP)$weight)
new_m=mean(S[S>=limits[1]],na.rm = T)
new_limit=sum(limits)/2
while (abs(new_m-m_p)>(new_limit/sum(S>new_limit,na.rm = T))){
  new_limit=sum(limits)/2
  print(new_limit)
  new_m=mean(S[S>=new_limit],na.rm = T)
  print(new_m)
  if (new_m>m_p){
    limits[1]=new_limit
  }else{
    limits[2]=new_limit
  }
  print(limits)
}

A_S=S
A_S[A_S<new_limit]=0
diag(A_S)=0
WSaP=graph_from_adjacency_matrix(A_S,mode="upper",weighted = T)
save(WSaP,file=paste(name,"_WSaP.RData",sep=""))


#Networks with matching number of edges

thr=sort(P,decreasing = T)[ecount(WSwS)*2]
A_P=P
A_P[A_P<thr]=0
diag(A_P)=0
WPeS=graph_from_adjacency_matrix(A_P,mode="upper",weighted = T)
save(WPeS,file=paste(name,"_WPeS.RData",sep=""))


thr=sort(S,decreasing = T)[ecount(WPwP)*2]
A_S=S
A_S[A_S<thr]=0
diag(A_S)=0
WSeP=graph_from_adjacency_matrix(A_S,mode="upper",weighted = T)
save(WSeP,file=paste(name,"_WSeP.RData",sep=""))

