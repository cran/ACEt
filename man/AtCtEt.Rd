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
AtCtEt(data_m, data_d, model = c('d','d','d'), knot_a = 5, knot_c = 5, knot_e = 5, 
boot=FALSE, num_b = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_m}{
An \eqn{N_m} x 3 data matrix for MZ twins. \eqn{N_m} is the number of MZ twin pairs. The first two columns are centered trait values (i.e. the mean should be zero) and the third column is age (or other covariates). 
}
  \item{data_d}{
An \eqn{N_d} x 3 data matrix for DZ twins. \eqn{N_d} is the number of DZ twin pairs. The first two columns are centered trait values (i.e. the mean should be zero) and the third column is age (or other covariates). 
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

  \item{boot}{
An indicator of whether to use the bootstrap method to calculate the confidence interval. The default is FALSE. 
}

  \item{num_b}{
The number of replicates when the bootstrap method is used (i.e. \eqn{boot=TRUE}). The default is 100.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{n_beta_a }{The number of spline coefficients for the A component.}
  \item{n_beta_c }{The number of spline coefficients for the C component.}
  \item{n_beta_e }{The number of spline coefficients for the E component.}
  \item{beta_a }{The estimated spline coefficients (if the model parameter is 'd') or variance (if the model parameter is 'c') of the A component.}
  \item{beta_c }{The estimated spline coefficients (if the model parameter is 'd') or variance (if the model parameter is 'c') of the C component.}
  \item{beta_e }{The estimated spline coefficients (if the model parameter is 'd') or variance (if the model parameter is 'c') of the E component.}
  \item{hessian }{The approximated Hessian matrix from the minus log-likelihood.}
  \item{con }{An indicator of convergence of the optimization algorithm. An integer code 0 indicates successful completion. See 'optim' for more details.}
  \item{lik }{The minus log-likelihood.}
  \item{knots_a }{A vector of the knot positions for the A component.}
  \item{knots_c }{A vector of the knot positions for the C component.}
  \item{knots_e }{A vector of the knot positions for the E component.}
  \item{boot }{A list containing confidence bands estimated from the bootstrap method when \eqn{boot=TRUE}}
%% ...
}
\references{
He, L., Sillanp‰‰, M. J., Silventoinen, K., Kaprio, J., & Pitk‰niemi, J. (2016). Estimating Modifying Effect of Age on Genetic and Environmental Variance Components in Twin Models. Genetics, genetics-115.
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
