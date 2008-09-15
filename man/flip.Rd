\name{flip}
\alias{flip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{a little function which is equal to flipud in matlab,make the matrix"upside-down"}
\description{
  Use recursive function
}
\usage{
flip(A)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{A}{ The input matrix: matrix to be flipped}
}

\value{
 ~Describe the value returned

  \item{A}{ the flipped matrix}
}
\author{Ann Lee, Di Liu}

\examples{

   x=matrix(seq(1,9,by=1),ncol=3)
   flip(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }

