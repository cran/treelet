Calc_Energy=function(X,basis,K,av_norm){
coeffs=X%*%t(basis);
energy=coeffs*coeffs;
avE=apply(energy,2,mean);  #each row is an observation
avE=avE/av_norm;   #normalize

#compute energy in descending order, get the first K values
avE=sort(avE,decreasing=TRUE);
score=sum(avE[1:K]);
return(list(score=score));
}

