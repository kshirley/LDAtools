\name{plotLoglik}
\alias{plotLoglik}
\title{Traceplot of log-likelihood.}
\usage{
plotLoglik(data = matrix(), zoom = FALSE, prop = 0.5, start = 0)
}
\arguments{
  \item{data}{a matrix with log-likelihood values. Each
  column of the matrix refers to a different MCMC chain.}

  \item{zoom}{should a second plot of the latter iterations
  be output.}

  \item{prop}{If \code{zoom == TRUE}, this number (between
  0 and 1) sets the proportion of final iterations to
  plot.}

  \item{start}{number of iterations (if any) performed
  before these chains.}
}
\description{
Trace plot of the log-likelihood (the x-axis is MCMC
iteration, and the y-axis is the log-likelihood). The MCMC
algorithm has converged when this plot levels out and the
log-likelihood time series plot looks like white noise. If
the log-likelihood is still increasing at the last
iteration, you should go back and run the algorithm longer
(using \code{fitLDA}.)
}
\examples{
one.chain <- matrix(rnorm(1000, mean=-500, sd=100), ncol=1) #hopefully our logliklihood is essentially random noise!
plotLoglik(one.chain)
par(mfrow=c(1,2))
plotLoglik(one.chain, zoom=TRUE)
three.chains <- matrix(rnorm(3000, mean=-500, sd=100), ncol=3)
plotLoglik(three.chains)
plotLoglik(three.chains, zoom=TRUE, prop=.1, start=5000)
}

