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
Version: \tab 1.3\cr
Date: \tab 2016-02-06\cr
License: \tab GPL \cr
}
}
\author{
Liang He

Maintainer: Liang He <liang.he@duke.edu>

}
\references{
Liang He, Mikko J. Sillanp‰‰, Karri Silventoinen, Jaakko Kaprio, Janne Pitk‰niemi, Estimating modifying effect of age on genetic and environmental variance components in twin models. Genetics, 2016
}
\keyword{ package }

\examples{

# data(data_ace)

# AtCtEt(data_ace$mz, data_ace$dz, model=c('d','d','c'), knot_a = 5, knot_c = 4)
# AtCtEp(data_ace$mz, data_ace$dz, knot_a = 8, knot_c = 7)

}
