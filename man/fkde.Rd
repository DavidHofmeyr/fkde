\name{fkde}
\alias{fkde}
\title{
Fast Kernel Density Estimation
}
\description{
Evaluates a kernel based estimate of a density function from a univariate sample. The default is to evaluate the density at all of the sample points. Optional evaluation on a grid is also available by setting argument ngrid to non-zero value.
}
\usage{
fkde(x, h, ord, ngrid, nbin, m, M, hnorm)
}
\arguments{
  \item{x}{vector of (univariate) sample observations.}
  \item{h}{(optional) bandwidth used in density estimate. If omitted then Silverman's rule of thumb is used.}
  \item{ord}{(optional) degree of polynomial part of the kernel function.}
  \item{ngrid}{(optional) number of evaluation points if evaluation on a grid is desired.}
  \item{nbin}{(optional) number of bins to use if approximate density estimate is desired instead of exact.}
  \item{m}{(optional) minimum evaluation point if grid evaluation or binned approximation required. Defaults to min(x)-6*h}
  \item{M}{(optional) maxmimum evaluation point if grid evaluation or binned approximation required. Defaults to max(x)+6*h}
  \item{hnorm}{(optional) to make implementation compatible with existing R functions for bandwidth selection, which usually by default provide the bandwidth value for the Gaussian kernel. If hnorm is provided then fkde will convert this to the corresponding value for the kernel being used, based on the AMISE value relative to that of the Gaussian.}
}
\examples{
# generate a sample of size 1000000 from heavily skewed density
x <- dataGen(1000000, 5)

# obtain exact kernel estimate and a binned approximation on a grid of 1000 points
system.time(f_exact <- fkde(x, ngrid = 1000))
system.time(f_binned <- fkde(x, nbin = 1000))

# estimate the MISE
MISE(5, f_exact$x, f_exact$y)
MISE(5, f_binned$x, f_binned$y)
}
