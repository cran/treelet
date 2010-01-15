\name{draw_tree}
\alias{draw_tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ draws the tree structure from first level to the highest level }
\description{
  needs library "ape", uses the result of function "Run\_JTree". 
}
\usage{
draw_tree(out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{The output of the function "Build\_JTree"}
}




\examples{

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

out=Run_JTree(X, dim(X)[2]-1)
draw_tree(out)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }

