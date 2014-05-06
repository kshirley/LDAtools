\name{fitLDA}
\alias{fitLDA}
\title{Fit LDA model via Gibbs sampler}
\usage{
fitLDA(word.id = integer(), doc.id = integer(), k = 10, n.chains = 1,
  n.iter = 1000, topics.init = NULL, alpha = 0.01, beta = 0.01)
}
\arguments{
  \item{word.id}{Unique token ID. Can be taken directly
  from the output of \code{filter}.}

  \item{doc.id}{Unique document ID. Can be taken directly
  from the output of \code{filter}.}

  \item{k}{number of topics.}

  \item{n.chains}{number of MCMC chains.}

  \item{n.iter}{number of iterations.}

  \item{topics.init}{A vector of topics to initially
  assign. The Markov property of MCMC allows one to input
  the topic assignments from the last iteration of a
  previous model fit. Note that this vector should be the
  same length of the \code{word.id} vector times the number
  of chains.}

  \item{alpha}{Dirichlet hyperparameter}

  \item{beta}{Dirichlet hyperparameter}
}
\value{
A list of length two. The first element is the sampled
latent topic value from the last iteration (for each
token). The second element is a vector with the
log-likelihood values for every iteration of the gibbs
sampler.
}
\description{
This function implements the Gibbs sampling method
described by Griffiths and Steyvers (2004). The Gibbs
sampler portion of the function is a call to C code. Note
that we only return the latent topic assignments (for each
token) from the last iteration. Thus, memory limitations
aren't really an issue. However, the run time is
O(num.chains*n.iter*N*k) where \code{n.chains} is number of
MCMC chains, \code{n.iter} is the number of iterations, N
is the total number of tokens in the data, and k is the
number of topics. It is possible to resume a Gibbs sampler
from a previous fit by using the topics from that fit to
initiate the next set of iterations using
\code{topics.init}.
}
\examples{
data(APinput)
#takes a while
\dontrun{o <- fitLDA(APinput$word.id, APinput$doc.id, k=20)}
}
\references{
Griffiths and Steyvers (2004). Finding Scientific Topics.
Proceedings of the National Academy of Sciences. 101:
5228-5235.
}

