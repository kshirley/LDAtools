\name{plotTokens}
\alias{plotTokens}
\title{Plot probable tokens for a given topic}
\usage{
  plotTokens(phi = vector(), vocab = character(),
    n.tokens = 20, lambda = 0.5, p = vector(), ...)
}
\arguments{
  \item{phi}{numeric vector with the probability of each
  token for a given topic.}

  \item{vocab}{character vector of the vocabulary for the
  corpus}

  \item{n.token}{the number of tokens to plot, where the
  default is \code{n.token = 20}.}

  \item{lambda}{a parameter between 0 and 1 to control how
  tokens are ranked within topics}

  \item{p}{the marginal probabilities of the tokens in the
  vocabulary}

  \item{...}{additional arguments to the plot() function}
}
\description{
  Plot probable tokens for a given topic
}
\details{
  The ranking of tokens within topics is based on a
  weighted average of the probability of a token (given the
  topic) and the lift, where the lift of a token is defined
  as the probability of the token (given the topic) divided
  by the marginal probability of the token (i.e. across all
  topics). The ranking that determines the top
  \code{n.token} tokens to plot is simply \code{lambda *
  log(p(token)) + (1 - lambda) * log(p(token |
  topic)/p(token))}.

  Note: the ordering of \code{phi}, \code{vocab}, and
  \code{p} must be the same (i.e. the nth element of each
  vector must correspond to the same token)
}
\examples{
data(APinput)
data(APtopics) #load output instead for demonstration
probs <- getProbs(word.id=APinput$word.id, doc.id=APinput$doc.id, topic.id=APtopics$topics,
               vocab=APinput$vocab)
 #THE ORDERING OF phi, vocab and p MUST MATCH!
tokens <- factor(APinput$vocab[APinput$word.id], levels=colnames(probs$phi.hat))
token.tab <- table(tokens)
p <- token.tab/sum(token.tab)
plotTokens(phi=probs$phi.hat[1,], vocab=names(p), n.tokens=30, lambda=1/3, p)
# plot all the topics!
\dontrun{
 for (i in seq_along(probs$phi.hat[,1])) {
   plotWords(probs$phi.hat[i,], tokens)
 }
}
}

