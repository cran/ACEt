\name{ACEt-package}
\alias{ACEt-package}
\alias{ACEt}
\docType{package}
\title{
Estimating Age Modification Effect on Genetic and Environmental Variance Components in Twin Models
}
\description{
The package implements two novel twin models that are able to estimate dynamic behaviour of the variance components in the classical twin models with respect to age or other covariates such as calendar year. 
}
\details{
\tabular{ll}{
Package: \tab ACEt\cr
Type: \tab Package\cr
Version: \tab 1.5.3\cr
Date: \tab 2016-03-27\cr
License: \tab GPL \cr
}
}
\author{
Liang He

Maintainer: Liang He <liang.he@duke.edu>

}
\references{
He, L., Sillanp‰‰, M. J., Silventoinen, K., Kaprio, J., & Pitk‰niemi, J. (2016). Estimating Modifying Effect of Age on Genetic and Environmental Variance Components in Twin Models. Genetics, genetics-115.
}
\keyword{ package }

\examples{

# data(data_ace)

# AtCtEt(data_ace$mz, data_ace$dz, mod=c('d','d','c'), knot_a = 5, knot_c = 4)
# AtCtEtp(data_ace$mz, data_ace$dz, knot_a = 6, mod=c('d','l','c'))

}
