\name{JTree_Basis}
\alias{JTree_Basis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate the basis for the Jacobi tree algorithm (one or more iterations)}
\description{
  calls function "flip"
}
\usage{
JTree_Basis(Zpos, T, PCidx, maxlev, all_nodes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Zpos}{positions of the two clusters merged;  number between 1 and d (the length of the comp vector)}
  \item{T}{Element in the list array contains a 2-by-2 rot matrix}
  \item{PCidx}{ID of prinicipal components; [1 2] or [2 1]}
  \item{maxlev}{default: dim-1 (levels L=0,1,2,...,maxlevel)}
  \item{all\_nodes}{node labels}
}

\value{
 ~Describe the value returned
  \item{basis}{Basis functions entered ROW-WISE --- coarse-to-fine where m is the dimension of the data}
  \item{sums}{Part of basis funcs of subspaces V1,...,V\_\{m-1\}---- entered rowwise (fine-to-coarse)}
  \item{difs}{ Basis funcs of subspaces W1,...,W\_\{m-1\}---- entered rowwise (fine-to-coarse)}
  
}

\author{Ann Lee, Di Liu}

\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }

