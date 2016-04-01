\name{plot_acet}
\alias{plot_acet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot variance curves for various ACE(t) or ACE(t)-p models with 95% CIs 
}
\description{
Plot variance curves or a heritability curve (with 95% CIs) of the A, C and E components with respect to age modelled by B-splines or P-splines.
}
\usage{
plot_acet(acet, boot = FALSE, heri = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acet}{
An object obtained from the B-splines or P-splines functions. For the P-splines functions, an object from the MCMC method must be used.
}
  \item{boot}{
An logical  indicator of whether the confidence bands estimated from the bootstrap method are plotted. The default is FALSE. Only available for the AtCtEt model.
}
  \item{heri}{
A logical indicator of whether to plot the dynamic heritability curve (currently only available for the B-splines). The default is FALSE.
}
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
# result <- AtCtEtp(data_ace$mz, data_ace$dz, knot_a = 5, knot_c = 4)
# result_mc <- acetp_mcmc(result, iter_num=10000, burnin = 500)
# plot_acet(result_mc)

# result <- AtCtEt(data_ace$mz, data_ace$dz, mod=c('d','c','c'), knot_a = 5)
# plot_acet(result)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
