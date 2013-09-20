\name{predictLDA}
\alias{predictLDA}
\title{Estimate topics for new documents using a Gibbs sampler}
\usage{
  predictLDA(word.id = integer(), doc.id = integer(),
    k = 10, n.chains = 1, n.iter = 1000,
    topics.init = NULL, alpha = 0.01, phi)
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

  \item{phi}{The \code{T} x \code{W} matrix containing the
  topic-token probability distributions for each of the
  \code{T} topics in the previously fit topic model.}
}
\value{
  A list of length two. The first element is the sampled
  latent topic value from the last iteration (for each
  token). The second element is a vector with the
  log-likelihood values for every iteration of the Gibbs
  sampler.
}
\description{
  This function estimates topic proportions for a new
  corpus of documents, using the the vocabulary and the
  topic-token probability distributions from a previously
  fit LDA topic model. The function samples the latent
  topics for each token in the new corpus using a Gibbs
  sampler, and returns the latent topics from the last
  iteration.
}
\examples{
data(APinput)
#takes a while
\dontrun{o <- fitLDA(APinput$word.id, APinput$doc.id, k=20)}
}

