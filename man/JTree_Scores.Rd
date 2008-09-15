\name{JTree_Scores}
\alias{JTree_Scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates scores for the "best K-basis" at each level of the tree}
\description{
  calls two function flip and Calc\_Energy
}
\usage{
JTree_Scores(X, K, Zpos, T, PCidx, maxlev, all_nodes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{data matrix}
  \item{K}{number of vectors in final representation}
  \item{Zpos}{Zpos[J,2] is the nodes to be merged}
  \item{T}{a list containing T nodes.Element in the cell array contains a 2-by-2 matrix that decorrelates inputs X(:,Z (k,1)) and X(:,Z(k,2));output from LKLB\_Decorrelate.}
  \item{PCidx}{PCidx[J,2]  id of principal components (1st or 2nd)}
  \item{maxlev}{default:dim-1}
  \item{all_nodes}{node labels, output of Build\_JTree}
}

\value{
  ~Discribe the output value
  \item{levels}{at each level}
  \item{scores}{what the score is at the corresponding level}
  
}

\author{Ann Lee, Di Liu}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }
