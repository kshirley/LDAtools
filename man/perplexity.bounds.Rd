\name{perplexity.bounds}
\alias{perplexity.bounds}
\title{Compute the lower bound of the perplexity of a topic model fit
using LDA}
\usage{
perplexity.bounds(term.id = integer(), alpha = double(), beta = double(),
  term.frequency = integer(), doc.frequency = integer(), print = 50)
}
\arguments{
  \item{term.id}{an integer vector containing the term ID
  number of every token in the corpus. Should take values
  between 1 and W, where W is the number of terms in the
  vocabulary.}

  \item{alpha}{the dirichlet prior parameter for the
  document-topic multionomial distributions}

  \item{beta}{the dirichlet prior parameter for the
  topic-term multionomial distributions}

  \item{term.frequency}{an integer vector containing the
  counts of the number of occurences of each term in the
  vocabulary}

  \item{doc.frequency}{an integer vector containing the
  number of tokens per document, whose length is equal to
  the total number of documents in the corpus.}
}
\value{
bounds a numeric vector of length two containing the upper
and lower bound of perplexity. The lower bound is computed
by the C function. The upper bound is simple -- it is just
the perplexity of a 1-topic model. Included as output just
for convenience.
}
\description{
This is an 'experimental' function that computes the lower
bound of the perplexity of the training data in an LDA
topic model. We claim that the perplexity of the training
data is minimized when alpha and beta, the priors for the
document-topic distributions and the topic-term
distributions, respectively, approach zero, and the number
of topics is equal to min(D, W), where D is the number of
training documents and W is the size of the vocabulary.
I'll have to write up the idea some day.
}

