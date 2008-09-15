\name{Get_C}
\alias{Get_C}
\title{a function to Get the (1) Covariance matrix (2) correlation matrix}
\description{
  Get both covariance matrix and correlation matrix of a given matrix. The matrix can not have missing data
}
\usage{
Get_C(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{input matrix, no missing value}
}

\value{
~Describe the value returned
   
  \item{C}{covariace matrix}
  \item{cc }{correlation matrix}

}

\author{Ann Lee, Di Liu}

\examples{
   x=matrix(seq(1,9,by=1),ncol=3)
   Get_C(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
