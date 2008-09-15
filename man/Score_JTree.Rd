\name{Score_JTree}
\alias{Score_JTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The function computes the scores for the best K basis, for each level of the tree }
\description{
calls three functions:Build\_JTree,JTree\_Scores,Get\_C.The returned variables are computed on a maximum level tree
}
\usage{
Score_JTree(X, K, TrainIdx, TestIdx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ data matrix (each row is an observation)}
  \item{K}{ number of vectors in final representation}
  \item{TrainIdx}{indices for train data (for building the tree)}
  \item{TestIdx}{indices for test data (for computing a score) }
}
\details{
 Remark: Compare this function with Run\_JTree
}
\value{
 ~Describe the value returned
 
  \item{levels }{at each level}
  \item{scores }{what the score is at the corresponding level}
  \item{Zpos }{Zpos[dim-1,2]   positions of the two clusters merged;  number between 1 and d (the length of the comp vector)}
  \item{ PCidx }{ID of prinicipal components; [1 2] or [2 1]}
  \item{ T}{Element T\{k\} in the list array contains a 2-by-2 rot matrix }
  \item{all\_nodes }{node labels }
  
}
\author{Ann Lee, Di Liu}

\examples{
N=1000;
dim_x=10;
sigma_x=1;


s1_2=290;
s2_2=300;
c1=-0.3;
c2=0.925;
u1=sqrt(s1_2)*matrix(rnorm(N),ncol=1);
u2=sqrt(s2_2)*matrix(rnorm(N),ncol=1);
u3=c1*u1+c2*u2;
s3_2=c1^2*s1_2+c2^2*s2_2;
Xraw=cbind(u1,u1,u1,u1,u2,u2,u2,u2,u3,u3);
#Xraw=cbind(rep(u1,4),rep(u2,4),u3,u3);
X=Xraw+sigma_x*matrix(rnorm(N*dim_x),ncol=dim_x);

###Just create a matrix X, whatever matrix you want


TrainIdx=sample(N,floor(3/4*N));#use 3/4 of the data for training data, randomize
TotalIdx=seq(1,N,by=1);
TestIdx=TotalIdx[-c(TrainIdx)]; # the others are testing data

K=5;# for example,,,,,,


out1=Score_JTree(X,K,TrainIdx,TestIdx);
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree}
