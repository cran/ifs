% file ifs/man/ifs.setup.Rd
% copyright (C) 2001-2 S.M. Iacus    
%
\name{ifs.setup}
\alias{ifs.setup}
\title{Set up the parameters for the maps of the IFS estimator}
\description{
 This is called before calling \code{setQF} to prepear the
 parameters to be passed in \code{setQF}.
}
\usage{
ifs.setup(y, maps = "quantile", qtl)
}
\arguments{
  \item{y}{the vector of the sample observations}
  \item{maps}{type of maps: quantile, wl1 or wl2}
  \item{qtl}{instead of passing the data y you can pass a vector
   of quantiles}
}
\value{
	\item{m}{the vector of the empirical moments}
	\item{a}{the vector of the coefficents `a' in the maps}	
	\item{s}{the vector of the coefficents `s' in the maps}	
	\item{n}{the number of maps}	
}
\author{S. M. Iacus}
\seealso{
\code{\link{ifs}}
}
\keyword{misc}


