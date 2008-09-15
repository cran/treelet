Run_JTree=function(X,maxlev){
print("computing the Correlation....")
out_C=Get_C(X);
C=out_C$C;
cc=out_C$cc;

out1=Build_JTree(C,cc,maxlev);
print("building the tree......")
out=JTree_Basis(out1$Zpos,out1$T,out1$PCidx,maxlev,out1$all_nodes);
print("computing the basis for the level......")
return(list(basis=out$basis,sums=out$sums,difs=out$difs,Zpos=out1$Zpos,T=out1$T,PCidx=out1$PCidx,all_nodes=out1$all_nodes))
}