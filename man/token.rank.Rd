\name{token.rank}
\alias{token.rank}
\title{Compute distinctiveness and saliency of the words in the vocabulary for a given topic model}
\usage{
token.rank(word.id = integer(), topic.id = integer(), phi.hat = numeric())
}
\arguments{
  \item{word.id}{an integer vector with the token id of
  each token occurrence in the data.}

  \item{topic.id}{an integer vector with the topic id for
  each token in the data.}

  \item{phi.hat}{a numeric (K x W) matrix with topic
  distributions on each of the K rows}
}
\value{
A list of two length-W vectors, the distinctiveness and
saliency of each of the W tokens in the vocab
}
\description{
This function computes 'distinctiveness' and 'saliency' as
defined in Chuang, et al (2012), 'Termite'.
}

