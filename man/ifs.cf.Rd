% file ifs/man/ifs.cf.Rd
% copyright (C) 2001-2 S.M. Iacus    
%
\name{ifs.cf}
\alias{ifs.cf}
\title{Calculates the main parameters of the IFS estimators}
\description{
  Tool function to construct and find the solution of the minimization
  problem involving the quadratic form \eqn{x'Qx + b'x}. Not an optimal one.
  You can provide one better then this.
}
\usage{
ifs.cf(Q,b)
}
\arguments{
  \item{Q}{the matrix \eqn{Q} of \eqn{x'Qx + b'x}}
  \item{b}{the vector \eqn{b} of \eqn{x'Qx + b'x}}
}
\value{
	\item{p}{the vector of the coefficients to be plugged into the IFS}
}
\seealso{
\code{\link{ifs}}
}
\keyword{misc}


