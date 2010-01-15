Run_JTree=function(X,maxlev, drawTree=FALSE){
print("computing the Correlation....")
C=var(X);
cc=cor(X);
print("building the tree......")
maxTree=Build_JTree(C,cc,maxlev);
print("computing the basis for the level......")
out=JTree_Basis(maxTree$Zpos,maxTree$T,maxTree$PCidx,maxlev,maxTree$all_nodes);
if(drawTree=="TRUE" & maxlev==(dim(X)[2]-1)){
	draw_tree(maxTree)
}
return(list(basis=out$basis,sums=out$sums,difs=out$difs,Z=maxTree$Z,Zpos=maxTree$Zpos,T=maxTree$T,PCidx=maxTree$PCidx,all_nodes=maxTree$all_nodes))
}