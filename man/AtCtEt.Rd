\name{AtCtEt}
\alias{AtCtEt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The ACE(t) model  
}
\description{
The ACE(t) model with the A, C and E variance components as functions with respect to age modelled by B-splines.
}
\usage{
AtCtEt(data_m, data_d, model = c('d','d','d'), knot_a = 5, knot_c = 5, knot_e = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_m}{
An \eqn{N_m} x 3 data matrix for MZ twins. \eqn{N_m} is the number of MZ twin pairs. The first two columns are traits and the third column is age (or other covariates). 
}
  \item{data_d}{
An \eqn{N_d} x 3 data matrix for DZ twins. \eqn{N_d} is the number of DZ twin pairs. The first two columns are traits and the third column is age (or other covariates). 
}
  \item{model}{
  A character vector of length 3. Each element specifies the function for the A, C or E component respectively. The A and C components can be 'd'(dynamic), 'c'(constant) or 'n'(NA). The E component can only be 'd' or 'c'. Thus, \eqn{model=c('c','c','c')} is corresponding to the classical ACE model.
  }
  \item{knot_a}{
The number of interior knots of the B-spline for the A component, which must be no less than 3. The default value is 5. 
}

  \item{knot_c}{
The number of interior knots of the B-spline for the C component, which must be no less than 3. The default value is 5.
}
  \item{knot_e}{
The number of interior knots of the B-spline for the E component, which must be no less than 3. The default value is 5.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{n_beta_a }{The number of spline coefficients for the A component.}
  \item{n_beta_c }{The number of spline coefficients for the C component.}
  \item{n_beta_e }{The number of spline coefficients for the E component.}
  \item{beta_a }{The estimated spline coefficients of the A component.}
  \item{beta_c }{The estimated spline coefficients of the C component.}
  \item{beta_e }{The estimated spline coefficients of the E component.}
  \item{hessian }{The approximated Hessian matrix from the minus log-likelihood.}
  \item{con }{The indicator of convergence of the optimization algorithm.}
  \item{lik }{The minus log-likelihood.}
  \item{knots_a }{A vector of the knot positions for the A component.}
  \item{knots_c }{A vector of the knot positions for the C component.}
  \item{knots_e }{A vector of the knot positions for the E component.}
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
data(data_ace)

result <- AtCtEt(data_ace$mz, data_ace$dz, model=c('d','d','c'))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
