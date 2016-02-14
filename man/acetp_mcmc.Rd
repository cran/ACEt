\name{acetp_mcmc}
\alias{acetp_mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute CIs for the ACE(t)-p model
}
\description{
Compute CIs for the ACE(t)-p model using the MCMC methods and the bootstrap method
}
\usage{
acetp_mcmc(acetp, iter_num = 10000, sd = 0.1, burnin = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acetp}{
An object of the AtCtEtp_model (or AtCtEp_model, AtEtp_model).
}
  \item{iter_num}{
The number of MCMC iteration.
}
  \item{sd}{
The standard error of the normal proposal distribution in the MCMC algorithm. The default value is 0.1.
}
  \item{burnin}{
The number of burn-in, which must be smaller than the number of iteration.
}
}

\value{

  \item{beta_a_mc }{The estimates of the spline coefficients for the A component based on the posterior mean from the MCMC method.}
  \item{beta_c_mc }{The estimates of the spline coefficients for the C component based on the posterior mean from the MCMC method.}
  \item{beta_e_mc }{The estimates of the spline coefficients for the E component based on the posterior mean from the MCMC method.}
  \item{cov_a_mc }{The posterior covariance matrix of the estimates of the spline coefficients for the A component.}
  \item{cov_c_mc }{The posterior covariance matrix of the estimates of the spline coefficients for the C component.}
  \item{cov_e_mc }{The posterior covariance matrix of the estimates of the spline coefficients for the E component.}
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

# result <- AtCtEp(data_ace$mz, data_ace$dz, knot_a = 7, knot_c = 7)
# result_mc <- acetp_mcmc(result, iter_num=10000, burnin = 500)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
