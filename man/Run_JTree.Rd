\name{Run_JTree}
\alias{Run_JTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The function computes the basis representation of the matrix at maxlev}
\description{
 calls three functions: Build\_JTree,Get\_C,JTree\_Basis
}
\usage{
Run_JTree(X, maxlev)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{data matrix (each row is an observation)}
  \item{maxlev}{default: dim-1 (levels L=0,1,2,...,maxlevel)}
}

\value{
~Describe the output values
  \item{basis }{Basis functions entered ROW-WISE --- coarse-to-fine where m is the dimension of the data}
  \item{sums }{ Part of basis funcs of subspaces V1,...,V\_\{m-1\}---- entered rowwise (fine-to-coarse)}
  \item{difs }{ Basis funcs of subspaces W1,...,W\_\{m-1\}---- entered rowwise (fine-to-coarse)}
  \item{Zpos }{positions of the two clusters merged;  number between 1 and d (the length of the comp vector)}
  \item{PCidx }{ID of prinicipal components; [1 2] or [2 1]}
  \item{T}{Element in the list array contains a 2-by-2 rot matrix }
  \item{all\_nodes}{ node labels }
}

\author{Ann Lee, Di Liu}

\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }

