\name{fkde-package}
\alias{fkde-package}
\alias{fkde-package}
\docType{package}
\title{
  Fast and exact kernel density and density derivative estimation
}
\description{
  Fast and exact density and density derivative estimation using kernels of the type poly(|x|)exp(-|x|).
}
\details{
  The main functions are fkde and dfkde which offer fast kernel density and density derivative estimation respectively. Currently dfkde only provides evaluation at the sample points, but will updated soon to include evaluation on a grid and a sped up approximate method using binning. Both functions require only a single argument, a vector containing a univariate sample of observations. Optional arguments can be seen in the documentation for the individual functions.
}
\author{
David P. Hofmeyr
}
