JTree_Basis=function(Zpos,T,PCidx,maxlev,all_nodes){

# input please refer to the output of Build_JTree
# Output:
#   basis(m,m)       Basis functions entered ROW-WISE --- coarse-to-fine 
#                    where m is the dimension of the data
#   sums(J,m)      Part of basis funcs of subspaces V1,...,V_{m-1} 
#                    ---- entered rowwise (fine-to-coarse)
#   difs(J,m)      Basis funcs of subspaces W1,...,W_{m-1}        
#                    ---- entered rowwise (fine-to-coarse)

J=dim(Zpos)[1];

#if (nargin<4###############needs to be fixed)
    #   m=J+1;
#else{
       m=dim(all_nodes)[2];   #dim of the data
       nodes=all_nodes[maxlev,];
       nodes=nodes[which(nodes>0)];          #find the existing nodes at maxlev
      # clear all_nodes;##########needs to be fixed
#}
tmpfilts=diag(rep(1,m));   #each row is a component
ind=list();
sums=matrix(rep(0,m*maxlev),ncol=m);
difs=matrix(rep(0,m*maxlev),ncol=m);

for (lev in 1:maxlev){

#the updating of tmpfilts actually:
# tmpfilts=tmptilts%*%R(j,k,theta)

   s=tmpfilts[Zpos[lev,],];   #S is a 2*m matrix, components to be merged
   R=T[[lev]];
   y=t(R)%*%s;
   tmpfilts[Zpos[lev,],]=y;
   y=y[PCidx[lev,],];  #y[1,]=PC1, y[2,]=PC2, sort
   sums[lev,]=y[1,];
   difs[lev,]=y[2,];
}

#if (nargin<4){
   # basis=rbind(sums[J,],flip(difs[(J-m+2):J,]));
#}
#else {
  tmp=rbind(diag(rep(1,m)),sums);
  basis=rbind(tmp[nodes,],difs[maxlev:1,]);
#}

return(list(basis=basis,sums=sums,difs=difs));
}
