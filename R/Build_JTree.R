
#this is used for building the tree. It combines all covariates until there is only one node. So level=p-1

Build_JTree=function(C,cc,maxlevel){
# input C is covariance matrix, cc is the correlation matrix, by default maxlevel=# covariates-1
# about the output
#  Zpos   positions of the two clusters merged;  number between 1 and d (the length of the comp vector)
#   T         Element T{k} in the list array contains a 2-by-2 rot matrix 
#   PCidx   ID of prinicipal components; [1 2] or [2 1]
#   all_nodes(dim-1,dim)   node labels
#   theta     rotation angles that decorrelates inputs X(:,Z(k,1)) and X(:,Z(k,2))
#   PC_ratio   ratio C_qq/C_pp, where p is the index of the 1st
#                       principal component and q is the index of the 2nd principal comp
#   Z         contain the indices of clusters that are merged; 
#   all_d  labels of the d-components in the comp vector (=0 if s-component, otherwise 1,...,d-1)
#   C      updated covariance matrix 
#   cc     correlation coefficients

dim_C=dim(C)[1];   #dim of C, which is supposed to be a square matrix

# the parts that I commented with "#" seems unnecessary now. But might be useful in the future

#too large matrix
#if (dim_C[1]>65535)    
     #print(ERROR:Too large matrix);

 #must be squared matrix
#if (dim_C[1]!=dim_C[2])   
      #print (Wrong dimensions of input matrix); 

#if (nargs() <3){
#  J=dim_C[1]-1;
#} else {
  # J=maxlevel;
#}

J=maxlevel;

#Create relative matrix

Z=matrix (rep(0,J*2),ncol=2);   # each row contains the rows to be merged
T=list()   # T is going to put the 2-by-2 rotation matrix for each level
theta=rep(0,J); #put the rotation angles
PCidx=matrix(rep(0,J*2),ncol=2) ; #Order of principal components [1,2] or [2,1]

#initialize

L=1;  ##nodes=dim-lev
maskno=matrix();    #mask out the merged elements
nodes=seq(1,dim_C,by=1)      #node indices
dlabels=rep(0,dim_C);           #lables of d-components
PC_ratio=rep(0,dim_C-1);       #ratio of 1st principal components and the 2nd
Zpos=matrix(rep(0,J*2),ncol=2);
all_d=matrix(rep(0,J*dim_C),ncol=dim_C);
all_nodes=matrix(rep(0,J*dim_C),ncol=dim_C);

#loop through scale levels

for (lev in 1:J) {

###newJacobi

#newJacobi=function(C,cc,maskno){
  #if (nargs() <3){
     #maskno=matrix();
     #if (nargs()<2){
     #d=diag(C);
     #cc=C./sqrt(d*t(d));
    # }
###############################################
#finding the maximum correlation of the matrix
###############################################

#dim_C=dim(C)[1];
#C_mask=cc;
mask_C=upper.tri(cc)*cc;

#exclude the merged ones

#mask the merged terms
#C_mask[maskno,]=-1;
#C_mask[,maskno]=-1;

#mask the terms on the diag
#for (i in 1:dim_C){
#C_mask[i,i]=-1;
#}
k=(mask_C==0);
mask_C[k]=-1; #extract upper triangle of correlation matrix
mask_C[maskno,]=-1;  
mask_C[,maskno]=-1; 

compno=which(mask_C==max(mask_C),arr.ind=TRUE);
#maxcc=max(C_mask)
#compno=temp[1,];

############################
#Decorrelate components in pair
#############################

Cred=C[compno,compno]#the 2 by 2 matrix that only contains the two covariates

if(Cred[1,2]==0){
Cnew=C;
ccnew=cc;
R=diag(c(1,1));
theta=0;
idx=c(1,2);
}else{
C11=Cred[1,1];
C22=Cred[2,2];
C12=Cred[1,2];
th=1/2*atan(2*C12/(C11-C22));

cs=cos(th);
sn=sin(th);

R=rbind(c(cs,-sn),c(sn,cs));#rotation matrix R
###################################
#update C. C<-R'CR=(R'C)R
###################################
M=C;
M[compno,]=t(R)%*%C[compno,];
C=M;
C[,compno]=M[,compno]%*%R;

#update Cred now
Cred=C[compno,compno]
#find out which is the first principal component vector,use it as the merged vector
idx=c(Cred[1,1],Cred[2,2]);
idx=sort.list(idx,decreasing = TRUE);

#################
#update cc
#################

#get the diagnal components, (dim,1)
dnew=diag(C);
temp=sqrt(matrix(dnew[compno],ncol=1)%*%dnew);
temp=C[compno,]/temp;
cc[compno,]=temp;
cc[,compno]=t(temp);
}

#return(list(Cnew=Cnew,ccnew=ccnew,R=R,th=th,compno=compno,maxcc=maxcc,idx=idx));
#}

#finishing updating

###############################
#store the updated results
###############################

#dist=(1-newJacobi$maxcc)/2; # not need it right now
PCidx[lev,]=idx;      #[1,2] or [2,1]
theta[lev]=th;           #rotation angle
T[[lev]]=R;                 #rotation angle

#C=newJacobi$Cnew;
#cc=newJacobi$ccnew;

Z[lev,]=nodes[compno];  #nodes that are merged

pind=compno[idx];
p1=pind[1];
p2=pind[2];
nodes[pind]=cbind(dim_C+lev,0);

dlabels[p2]=lev;
if(lev==1){
maskno=p2;
maskno=as.matrix(maskno);
}else{
maskno=cbind(maskno,p2);               #update the maskno
}
PC_ratio[lev]=C[p2,p2]/C[p1,p1];          #see the definition of PC_ratio
Zpos[lev,]=compno;
all_d[lev,]=t(dlabels);
all_nodes[lev,]=nodes;
}

return(list(Z=Z,Zpos=Zpos,T=T,PCidx=PCidx,all_nodes=all_nodes,theta=theta,PC_ratio=PC_ratio,Z=Z,all_d=all_d,C=C,cc=cc))

}




