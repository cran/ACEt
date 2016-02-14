\name{AtCtEp}
\alias{AtCtEp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The ACE(t)-p model 
}
\description{
The ACE(t)-p model with the A and C variance components as functions with respect to age modelled by P-splines.
}
\usage{
AtCtEp(data_m, data_d, knot_a = 8, knot_c = 8, eps = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_m}{
An \eqn{N_m} x 3 data matrix for MZ twins. \eqn{N_m} is the number of MZ twin pairs. The first two columns are traits and the third column is age (or other covariates). 

}
  \item{data_d}{
An \eqn{N_d} x 3 data matrix for DZ twins. \eqn{N_d} is the number of DZ twin pairs. The first two columns are traits and the third column is age (or other covariates). 
}
  \item{knot_a}{
The number of interior knots for the A component. The default value is 8.
}
  \item{knot_c}{
The number of interior knots for the C component. The default value is 8.
}
  \item{eps}{
Tolerance for convergence of the EM algorithm iterations. The default value is 0.1.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  
  \item{var }{The estimated variance of the E component.}
  \item{var_b_a }{The estimated variance for the penalized coefficient for the A components.}
  \item{var_b_c }{The estimated variance for the penalized coefficient for the C components.}
  \item{beta_c }{The estimated spline coefficients of the A component.}
  \item{beta_a }{The estimated spline coefficients of the C component.}
  \item{con }{The indicator of convergence of the optimization algorithm.}
  \item{lik }{The minus log marginal likelihood.}
  \item{knot_a }{The values of the knots for the A component.}
  \item{knot_c }{The values of the knots for the C component.}
%% ...
}
\references{
Liang He, Mikko J. Sillanp‰‰, Karri Silventoinen, Jaakko Kaprio, Janne Pitk‰niemi, Estimating modifying effect of age on genetic and environmental variance components in twin models. Genetics, 2016
}
\author{
Liang He
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# data(data_ace)

# result <- AtCtEp(data_ace$mz, data_ace$dz, knot_a = 8, knot_c = 8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
