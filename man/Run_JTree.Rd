\name{Run_JTree}
\alias{Run_JTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The function computes the basis representation of the matrix at maxlev}
\description{
 It is the main function of the treelet algorithm.
 calls three functions: Build\_JTree,JTree\_Basis
}
\usage{
Run_JTree(X, maxlev, drawTree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{data matrix (each row is an observation)}
  \item{maxlev}{default: dim-1 (levels L=0,1,2,...,maxlevel)}
  \item{drawTree}{draw the tree structure if we build the whole tree. Use "TRUE"}
}

\value{

  \item{basis }{Basis functions entered ROW-WISE --- coarse-to-fine where m is the dimension of the data}
  \item{sums }{ Part of basis funcs of subspaces V1,...,V\_\{m-1\}---- entered rowwise (fine-to-coarse)}
  \item{difs }{ Basis funcs of subspaces W1,...,W\_\{m-1\}---- entered rowwise (fine-to-coarse)}
  \item{Zpos }{positions of the two clusters merged;  number between 1 and d (the length of the comp vector)}
  \item{PCidx }{ID of principal components; [1 2] or [2 1]}
  \item{T}{Element in the list array contains a 2-by-2 rot matrix }
  \item{all\_nodes}{ node labels }
  \item{Z}{indices of clusters that are merged}
  \item{C\_old}{the covariance matrix}
  \item{cc\_old}{the old correlation matrix}
  \item{basis\_cov}{the same as the basis\_cov as returned by JTree_Basis}
}


\examples{

#Generating the Sample Data
N=1000;
dim_x=10;
sigma_x=1;
set.seed(100);

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
X=Xraw+sigma_x*matrix(rnorm(N*dim_x),ncol=dim_x)

#Run_JTree, Give the full treelet decomposition
k=1;
maxlev=dim(X)[2]-k;
output=Run_JTree(X,maxlev);
# Display the output image
#library(matlab)
#imagesc(output$basis)
#Display the output in tree structure
#draw_tree(output)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }

