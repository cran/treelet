\name{Build_JTree}
\alias{Build_JTree}
\title{ Build a hierarchical/wavelet tree of 2-by-2 Jacobi rotations}
\description{
  Builds a hierarchical binary tree using the covariance structure of the inputs. Decorrelates groups of two with PCA .
  
}
\usage{
Build_JTree(C, cc, maxlevel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C}{ input covariance matrix}
  \item{cc}{correlation coeffs}
  \item{maxlevel}{ default: dim-1 (levels L=0,1,2,...,maxlevel) }
}

\value{
~Describe the value returned
 
  \item{Zpos }{  positions of the two clusters merged;  number between 1 and d (the length of the comp vector)}
  \item{PCidx }{ID of prinicipal components; [1 2] or [2 1]}
  \item{T}{Element T\{k\} in the list array contains a 2-by-2 rot matrix }
  \item{all_nodes }{node labels }
  \item{theta}{rotation angles that decorrelates inputs X(:,Z(k,1)) and X(:,Z(k,2)) }
  \item{Z }{contain the indices of clusters that are merged}
  \item{all_d }{labels of the d-components in the comp vector (=0 if s-component, otherwise 1,...,d-1)}
  \item{C }{updated covariance matrix}
  \item{cc }{updated coefficients }

}

\author{Ann Lee, Di Liu}

\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }

