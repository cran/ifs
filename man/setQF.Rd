% file ifs/man/setQF.Rd
% copyright (C) 2001-2 S.M. Iacus    
%
\name{setQF}
\alias{setQF}
\title{Set up the Quadratic Form}
\description{
  Tool function to construct the quadratic form \eqn{x'Qx + b'x} to be minimized
  to construct the IFS estimator.
}
\usage{
setQF(m, s, a, n = 10)
}
\arguments{
  \item{m}{the vector of the sample or true moments}
  \item{s}{the vector of coefficients \eqn{s_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{a}{the vector of coefficients \eqn{a_i} in: \eqn{w_i = s_i *x + a_i}}
  \item{n}{number of parameter to use in the IFS estimator, default = 10}
}
\value{
	\item{Q}{the matrix of the quadratic form}
	\item{b}{the matrix of the quadratic form}	
}
\details{
  This estimator is intended to estimate the continuous distribution
  function of a random variable on [0,1]. The estimator is a continuous
  function not everywhere differentiable.
 }
\author{S. M. Iacus}
\references{
  S. M. Iacus and D. La Torre (2001)
  \emph{Approximating distribution function by Iterated Function Systems.}
  submitted.
}
\seealso{
\code{\link{ifs}}
}
\keyword{nonparametric}

