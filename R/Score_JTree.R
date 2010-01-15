Score_JTree=function(X,K,TrainIdx,TestIdx){

# The function computes the scores for the best K basis, for each level of the tree
# X is the whole matrix we are going to work on
# use TrainIdx and TestIdx to split Training and Testing examples
# Returns the score of the tree on each level. 
# This function calls: Build_JTree and JTree_Scores



############################
#Building the Tree
############################


Xtrain=X[TrainIdx,]
dim_Xtrain=dim(Xtrain)[2];

#covarice and correlation matrix
print("computing correlations for training data......")

C=var(Xtrain);
cc=cor(Xtrain);

#run treelet algorithm
maxlev=dim_Xtrain-1;  # the maximum of levels is # covariates -1
print("Building the tree......")
out_train=Build_JTree(C,cc,maxlev);# Please refer to the function Build_JTree for the returned values

##################################
#computing scores on testing data
##################################
print("Geting the score on testing data......")
Xtest=X[TestIdx,];
out=JTree_Scores(Xtest, K,out_train$Zpos,out_train$T,out_train$PCidx,maxlev,out_train$all_nodes);
scores=out$scores;
levels=out$levels;

return(list(socres=scores,levels=levels,Zpos=out_train$Zpos,T=out_train$T, PCidx=out_train$PCidx,all_nodes=out_train$all_nodes)); #this is what we want!

}






