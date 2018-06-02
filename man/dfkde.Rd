\name{dfkde}
\alias{dfkde}
\title{
Fast Kernel Density Derivative Estimation
}
\description{
Evaluates a kernel based estimate of the first derivative of a density function from a univariate sample. Currently only evaluation at the sample points is implemented. Updates to evaluation on a grid and using a fast binned approximation will be added soon.
}
\usage{
dfkde(x, h, ord, hnorm)
}
\arguments{
  \item{x}{vector of (univariate) sample observations.}
  \item{h}{(optional) bandwidth used in density derivative estimate. If omitted then Silverman's rule of thumb is used.}
  \item{ord}{(optional) degree of polynomial part of the kernel function.}
  \item{hnorm}{(optional) to make implementation compatible with existing R functions for bandwidth selection, which usually by default provide the bandwidth value for the Gaussian kernel. If hnorm is provided then dfkde will convert this to the corresponding value for the kernel being used, based on the AMISE value relative to that of the Gaussian.}
}
\examples{
# generate a sample of size 1000000 from simple bimodal density
x <- dataGen(1000000, 4)

# obtain exact kernel estimate of the density using order 1 and order 7 kernel
system.time(df_7 <- dfkde(x, ord = 7))
system.time(df_1 <- dfkde(x))

# plot exact derivative and kernel estimates (plotting of 1000000 points takes time)
xs = seq(-3.5, 3.5, length = 500)
plot(xs, .5/sqrt(2*pi)*27/8*((1-xs)*exp(-(xs-1)^2*9/8) - (xs+1)*exp(-(xs+1)^2*9/8)), type = 'l')

lines(df_7$x[(1:1000)*990+500], df_7$y[(1:1000)*990+500], col = 2)
lines(df_1$x[(1:1000)*990+500], df_1$y[(1:1000)*990+500], col = 3)
}
