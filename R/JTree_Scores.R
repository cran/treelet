JTree_Scores=function(X, K, Zpos,T,PCidx,maxlev,all_nodes){

#Calculates scores for the "best K-basis" at each level of the tree
#input X is the data matrix
#Zpos is the nodes to be merged
#T is a list. Each element is the R matrix at each level
#PCidx is the id of the principal components

#score is the scores for the best K basis
#calls Calc_Energy function


J=dim(Zpos)[1];
if (nargs()<4){
     m=J+1;
     maxlev=m-1;
    
}else{
       m=dim(all_nodes)[2];
}

#initialize

tmpfilts=diag(rep(1,m));   #each row is a component
sums=matrix(rep(0,m*maxlev),ncol=m);
difs=matrix(rep(0,m*maxlev),ncol=m);
levels=seq(0,maxlev,by=1);
scores=rep(0,maxlev);
av_norm=mean(apply(X*X,1,sum));# normalization term, E||x||^2
X_prime=X[,1:K];
scores[1]=mean(apply(X_prime*X_prime,1,sum))/av_norm;  # score at level L=0

for (lev in 1:maxlev){# loop through the levels
   s=tmpfilts[Zpos[lev,],];#s(2,m); components to be merged
   R=T[[lev]];# R is the rotation matrix
   y=t(R)%*%s;
   tmpfilts[Zpos[lev,],]=y;
   y=y[PCidx[lev,],];  #y[1,]=PC1, y[2,]=PC2, sort
   sums[lev,]=y[1,];
   difs[lev,]=y[2,];
    nodes=all_nodes[lev,];
       nodes=nodes[which(nodes>0)];

  tmp=rbind(diag(rep(1,m)),sums);
  basis=rbind(tmp[nodes,],difs[lev:1,]);

#score the best K-basis
  out_E=Calc_Energy(X,basis,K,av_norm);
  scores[lev+1]=out_E$score;
  
}

return(list(scores=scores,levels=levels));
}

