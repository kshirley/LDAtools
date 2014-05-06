\name{collapse.bigrams}
\alias{collapse.bigrams}
\title{Replace specified bigrams with terms representing the bigrams}
\usage{
collapse.bigrams(bigrams = character(), doc.id = integer(),
  term.id = integer(), vocab = character())
}
\arguments{
  \item{bigrams}{A character vector, each element of which
  is a bigram represented by two terms separated by a
  hyphen, such as 'term1-term2'. Every consecutive
  occurrence of 'term1' and 'term2' in the data will be
  replaced by a single token representing this bigram.}

  \item{doc.id}{an interger vector containing the document
  ID number of every token in the corpus. Should take
  values between 1 and D, where D is the total number of
  documents in the corpus.}

  \item{term.id}{an integer vector containing the term ID
  number of every token in the corpus. Should take values
  between 1 and W, where W is the number of terms in the
  vocabulary.}

  \item{vocab}{a character vector of length W, containing
  the terms in the vocabulary. This vector must align with
  \code{term.id}, such that a term.id of 1 indicates the
  first element of \code{vocab}, a term.id of 2 indicates
  the second element of \code{vocab}, etc.}
}
\value{
Returns a list of length three. The first element,
\code{new.vocab}, is a character vector containing the new
vocabulary. The second element, \code{new.term.id} is the
new vector of term ID numbers for all tokens in the data,
taking integer values from 1 to the length of the new
vocabulary. The third element is \code{new.doc.id}, which
is the new version of the document id vector. If any of the
specified bigrams were present in the data, then
\code{new.term.id} and \code{new.doc.id} will be shorter
vectors than the original \code{term.id} and \code{doc.id}
vectors.
}
\description{
After tokenization, use this function to replace all
occurrences of a given bigram with a single token
representing the bigram, and 'delete' the occurrences of
the two individual tokens that comprised the bigram (so
that it is still a generative model for text).
}

