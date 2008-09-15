# Get_C is a function to Get the (1) Covariance matrix (2) correlation matrix

Get_C=function(X){
N=dim(X)[1];
dim_X=dim(X)[2];
X=X-matrix(rep(1,N),ncol=1)%*%apply(X,2,mean); #center the data
tmp_norm=sqrt(apply(X*X,2,sum));
Xnorm=X/(matrix(rep(1,N),ncol=1)%*%tmp_norm)   #normalize columns
cc=t(Xnorm)%*%Xnorm;
C=t(X)%*%X/N;
remove(tmp_norm);
remove(Xnorm);
return(list(C=C,cc=cc));
}